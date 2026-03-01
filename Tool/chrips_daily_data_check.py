# -*- coding: utf-8 -*-
"""
檢查 CHIRPS p05 NetCDF 檔案完整性
1. 比對 Content-Length
2. 嘗試用 xarray 開啟並讀取一小部分數據
"""

import os
import requests
import xarray as xr

BASE_URL = "https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p05"
DATA_DIR = "../Data_Set/chirps_p05_nc"   # 你的檔案路徑

def check_file(fn, base_url=BASE_URL, data_dir=DATA_DIR, timeout=30):
    path = os.path.join(data_dir, fn)
    url = f"{base_url}/{fn}"

    if not os.path.exists(path):
        return (fn, "MISSING", "檔案不存在")

    local_size = os.path.getsize(path)

    # --- HEAD 檢查 Content-Length ---
    try:
        h = requests.head(url, timeout=timeout)
        h.raise_for_status()
        remote_size = int(h.headers.get("Content-Length", 0))
    except Exception as e:
        return (fn, "WARN", f"無法取得 Content-Length ({e})")

    if remote_size and local_size != remote_size:
        return (fn, "CORRUPT", f"大小不符，本地 {local_size}, 遠端 {remote_size}")

    # --- 嘗試開檔 ---
    try:
        ds = xr.open_dataset(path)
        # 嘗試讀取一個時間點的數值
        sample = ds["precip"].isel(time=0).values
        ds.close()
    except Exception as e:
        return (fn, "CORRUPT", f"無法讀取 NetCDF ({e})")

    return (fn, "OK", f"大小 {local_size} bytes")

def main():
    all_files = sorted([f for f in os.listdir(DATA_DIR) if f.endswith(".nc")])
    if not all_files:
        print("❌ 資料夾裡沒有 .nc 檔")
        return

    results = []
    for f in all_files:
        results.append(check_file(f))

    print("\n=== 檢查結果 ===")
    for fn, status, msg in results:
        print(f"{status:8} {fn} : {msg}")

if __name__ == "__main__":
    main()