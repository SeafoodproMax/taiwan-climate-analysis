# GlobalWarming

## Overview
Class final project on global warming trends using open datasets. The main focus
is NASA GISTEMP temperature anomalies and Taiwan-region analyses
(119–123E, 21–26N), with rainfall trend and extreme-event indicators.

## Data Sources
- **NASA GISTEMP (Land+Ocean, 1200 km)**  
  `Data_Set/gistemp1200_GHCNv4_ERSSTv5.nc`
- **CHIRPS v2.0 Daily Precipitation (p05)**  
  Downloaded to `Data_Set/chirps_p05_nc/` via `Tool/chrips_tw_daily.py`
- **Taiwan gridded monthly rainfall CSVs**  
  `Data_Set/RainFall/ObsRain_*.csv`

## Methods (Concise)
- **Temperature anomaly (GISTEMP)**: area-mean time series for Taiwan box with
  5-year (60-month) moving average.
- **Extreme rainfall (CHIRPS)**: Rx1day / Rx5day and exceedance counts for
  200/350/500 mm; compared with annual temperature.
- **Monthly rainfall**: Taiwan-wide monthly mean and 10-year moving average.

## Project Structure
- `R/TaiwanLandTemp.R` temperature anomaly analysis and plot
- `R/RainFall_daily.R` extreme rainfall analysis and plots
- `R/RainFall_monthly.R` monthly rainfall aggregation and plot
- `Tool/chrips_tw_daily.py` CHIRPS download + CSV preprocessing
- `Plot/` example figures

## How to Run
### 1) Download and build CHIRPS rainfall datasets
```bash
python Tool/chrips_tw_daily.py
```

### 2) Run R analyses
```bash
Rscript R/TaiwanLandTemp.R
Rscript R/RainFall_daily.R
Rscript R/RainFall_monthly.R
```

## Dependencies
### R
- ncdf4, zoo, ggplot2
- terra, dplyr, readr, lubridate, stringr, purrr, tidyr

### Python
- requests, tqdm, xarray, rioxarray, pandas, numpy

## Notes
- Scripts contain hard-coded `setwd()` / `os.chdir()` paths; update to your
  local project path or remove them.
- To analyze a different region, adjust the lat/lon bounds in the scripts.

