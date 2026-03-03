# -*- coding: utf-8 -*-
"""
CHIRPS v2.0 global_daily p05 批次下載 + 台灣逐日平均 + 極端指標
Requirements: requests, tqdm, xarray, rioxarray, pandas, numpy
pip install requests tqdm xarray rioxarray pandas numpy
"""
import os
import time
from enum import Enum
from pathlib import Path
import numpy as np
import pandas as pd
import requests
import xarray as xr
from tqdm import tqdm

# Set Working Directory
PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT_ROOT)

# ========== Download Settings ==========
BASE_URL = "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p05"
# annual list
YEARS = list(range(1981, 2025))  # 1981–2024
Base_Data_Dir = "Data_Set"
RAW_DATA_DIR = f"{Base_Data_Dir}/chirps_p05_nc"

# ========== Grid of Taiwan ==========
LON_MIN, LON_MAX = 119, 123
LAT_MIN, LAT_MAX = 21, 26

class RainfallAlert(Enum):
    """Taiwan-style rainfall alert thresholds based on 24h accumulation (mm).
    Values are millimeters in 24 hours.
    """
    HEAVY_200 = 200       # 豪雨
    TORRENTIAL_350 = 350  # 大豪雨
    EXTREME_500 = 500     # 超大豪雨

ALERT_LABELS = {
    RainfallAlert.HEAVY_200: "days_ge_200mm",
    RainfallAlert.TORRENTIAL_350: "days_ge_350mm",
    RainfallAlert.EXTREME_500: "days_ge_500mm",
}

# ========== Download Tools ==========
def download_one(year, out_dir=RAW_DATA_DIR, base_url=BASE_URL, retries=3, timeout=60):
    os.makedirs(out_dir, exist_ok=True)
    fn = f"chirps-v2.0.{year}.days_p05.nc"  #fn = filename
    url = f"{base_url}/{fn}"
    out_path = os.path.join(out_dir, fn)

    # Skip if exists
    if os.path.exists(out_path) and os.path.getsize(out_path) > 10_000:
        print(f"[skip] {fn} 已存在")
        return out_path

    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                total = int(r.headers.get("Content-Length", 0))
                with open(out_path, "wb") as f, tqdm(
                    total=total, unit="B", unit_scale=True, desc=fn
                ) as pbar:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
            # 粗略驗證大小
            if os.path.getsize(out_path) < 10_000:
                raise IOError("下載到的檔案大小異常")
            return out_path
        except Exception as e:
            print(f"[retry {attempt}/{retries}] {fn} 失敗：{e}")
            time.sleep(2 * attempt)
    raise RuntimeError(f"下載多次失敗：{fn}")

# ========== 讀檔、裁切、台灣逐日平均、極端指標 ==========
def build_tw_daily_csv(nc_dir=RAW_DATA_DIR, csv_dir=Base_Data_Dir,
                       lon_min=LON_MIN, lon_max=LON_MAX,
                       lat_min=LAT_MIN, lat_max=LAT_MAX,
                       out_daily_csv=None,
                       out_extreme_csv=None):
    if out_daily_csv is None:
        out_daily_csv = os.path.join(csv_dir, "tw_daily_mean.csv")
    if out_extreme_csv is None:
        out_extreme_csv = os.path.join(csv_dir, "tw_extreme_indices.csv")

    # Skip if both outputs already exist
    if os.path.exists(out_daily_csv) and os.path.getsize(out_daily_csv) > 0 \
       and os.path.exists(out_extreme_csv) and os.path.getsize(out_extreme_csv) > 0:
        print(f"[skip] 已存在：{out_daily_csv} 與 {out_extreme_csv}")
        return

    nc_files = sorted(
        [os.path.join(nc_dir, f) for f in os.listdir(nc_dir) if f.endswith(".nc")]
    )
    if not nc_files:
        raise FileNotFoundError("找不到 .nc 檔，請先執行下載。")

    print(f"讀取 {len(nc_files)} 個年度檔，逐檔處理台灣逐日平均")

    daily_frames = []
    for fn in tqdm(nc_files, desc="處理年度檔"):
        ds = xr.open_dataset(fn)
        try:
            tw = ds.sel(longitude=slice(lon_min, lon_max),
                        latitude=slice(lat_min, lat_max))
            # Area-weighted daily mean (cos(lat))
            weights = np.cos(np.deg2rad(tw["latitude"]))
            series = tw["precip"].weighted(weights).mean(dim=["latitude", "longitude"], skipna=True)
            df = series.to_dataframe().reset_index()[["time", "precip"]]
            df.rename(columns={"time": "Date", "precip": "Precip_mm"}, inplace=True)
            daily_frames.append(df)
        finally:
            ds.close()

    # Combine all daily data
    df_daily = pd.concat(daily_frames, ignore_index=True)
    df_daily.sort_values("Date", inplace=True)
    df_daily.to_csv(out_daily_csv, index=False, float_format="%.3f")
    print(f"[OK] 台灣逐日平均已輸出：{out_daily_csv}  （mm/day）")

    # ====== Extreme index：Rx1day / Rx5day（Yearly） ======
    df_daily["Year"] = pd.to_datetime(df_daily["Date"]).dt.year
    # Rx1day：Most rainfall day annually
    rx1 = (
        df_daily.groupby("Year", as_index=False)["Precip_mm"].max()
        .rename(columns={"Precip_mm": "Rx1day"})
    )
    # Rx5day：Most consecutive 5-day rainfall in the year
    df_daily = df_daily.sort_values(["Year", "Date"]).copy()
    df_daily["Roll5"] = (
        df_daily.groupby("Year")["Precip_mm"].rolling(5, min_periods=1).sum()
        .reset_index(level=0, drop=True)
    )
    rx5 = (
        df_daily.groupby("Year", as_index=False)["Roll5"].max()
        .rename(columns={"Roll5": "Rx5day"})
    )
    extreme = rx1.merge(rx5, on="Year", how="outer").sort_values("Year")
    extreme.to_csv(out_extreme_csv, index=False, float_format="%.3f")
    print(f"[OK] 極端指標已輸出：{out_extreme_csv}  （mm）")

def build_exceedance_counts(nc_dir=RAW_DATA_DIR, csv_dir=Base_Data_Dir,
                            lon_min=LON_MIN, lon_max=LON_MAX,
                            lat_min=LAT_MIN, lat_max=LAT_MAX,
                            levels=None, out_exceed_csv=None):
    """
    Build annual exceedance counts for Taiwan daily mean precipitation.

    Parameters
    ----------
    nc_dir : str
        Directory containing yearly CHIRPS p05 NetCDF files.
    csv_dir : str
        Directory to place output CSV.
    lon_min, lon_max, lat_min, lat_max : float
        Bounding box for Taiwan subset.
    levels : list[RainfallAlert]
        Which alert thresholds to count. Defaults to [350mm, 500mm].
    out_exceed_csv : str
        Output CSV path. Default: Data_Set/tw_exceedance_counts.csv

    Output schema
    -------------
    Year, days_ge_200mm?, days_ge_350mm?, days_ge_500mm?
    Columns appear only for requested levels.
    """
    if levels is None:
        levels = [RainfallAlert.TORRENTIAL_350, RainfallAlert.EXTREME_500]

    if out_exceed_csv is None:
        out_exceed_csv = os.path.join(csv_dir, "tw_exceedance_counts.csv")

    nc_files = sorted(
        [os.path.join(nc_dir, f) for f in os.listdir(nc_dir) if f.endswith(".nc")]
    )
    if not nc_files:
        raise FileNotFoundError("找不到 .nc 檔，請先執行下載。")

    print(f"讀取 {len(nc_files)} 個年度檔，建立多檔資料集以計算超閾值日數…")
    ds = xr.open_mfdataset(nc_files, combine="by_coords")

    # Subset Taiwan grid
    tw = ds.sel(longitude=slice(lon_min, lon_max),
                latitude=slice(lat_min, lat_max))

    print("採用『全台格點日最大值』計算超閾值日數（非面平均）。")
    # Daily maximum over Taiwan grid (closer to alert definitions than area-mean)
    # 這裡取島內任一格點的當日最大降雨，對應「任一地點達到警戒」
    daily_max = tw["precip"].max(dim=["latitude", "longitude"], skipna=True)

    # To pandas and add Year
    df_daily = daily_max.to_dataframe().reset_index()[["time", "precip"]]
    df_daily.rename(columns={"time": "Date", "precip": "Precip_mm"}, inplace=True)
    df_daily.sort_values("Date", inplace=True)
    df_daily["Year"] = pd.to_datetime(df_daily["Date"]).dt.year

    # Prepare output frame with all years present in the daily series
    years = (
        df_daily["Year"].dropna().astype(int).sort_values().unique().tolist()
    )
    out = pd.DataFrame({"Year": years})

    # Count exceedance days for each requested alert level
    print("計算各等級超閾值日數…")
    for lvl in tqdm(levels, desc="Exceedance levels"):
        label = ALERT_LABELS.get(lvl, f"days_ge_{int(lvl.value)}mm")
        mask = df_daily["Precip_mm"] >= float(lvl.value)
        cnt = (df_daily.assign(__hit=mask)
               .groupby("Year", as_index=False)["__hit"].sum()
               .rename(columns={"__hit": label}))
        out = out.merge(cnt, on="Year", how="left")

    # Sort and write
    out = out.sort_values("Year")
    out.to_csv(out_exceed_csv, index=False)
    print(f"[OK] 年度超閾值日數已輸出：{out_exceed_csv}")

    return out

def main():
    # 1) Downloading CHIRPS p05 annual files
    print("== 下載 CHIRPS p05 年度檔 ==")
    for y in YEARS:
        download_one(y)

    # 2) Generating Taiwan daily average & extreme indices (csv)
    print("== 建立台灣逐日平均 & 極端指標 ==")
    build_tw_daily_csv()

    # 3) Annual exceedance counts for heavy rain alert thresholds
    print("== 建立各等級（大豪雨/超大豪雨）年度超閾值日數表 ==")
    build_exceedance_counts(
        levels=[
            RainfallAlert.HEAVY_200,  # 豪雨 ≥200mm/24h
            RainfallAlert.TORRENTIAL_350,  # 大豪雨 ≥350mm/24h
            RainfallAlert.EXTREME_500      # 超大豪雨 ≥500mm/24h
        ]
    )

if __name__ == "__main__":
    main()
