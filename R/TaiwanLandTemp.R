# --- Packages ---
library(ncdf4)
library(zoo)
library(ggplot2)
library(here)

# --- PATH ---
setwd(here())

# --- Read NetCDF ---
nc <- nc_open("Data_Set/gistemp1200_GHCNv4_ERSSTv5.nc")
lon      <- ncvar_get(nc, "lon")
lat      <- ncvar_get(nc, "lat")
time_bnd <- ncvar_get(nc, "time_bnds")        # [2, time]
anom_raw <- ncvar_get(nc, "tempanomaly")      # [lon, lat, time]

# 屬性：縮放與缺值
sf   <- ncatt_get(nc, "tempanomaly", "scale_factor")
ao   <- ncatt_get(nc, "tempanomaly", "add_offset")
fill <- ncatt_get(nc, "tempanomaly", "_FillValue")$value
nc_close(nc)

# --- TimeLine ---
dates <- as.Date("1800-01-01") + colMeans(time_bnd)  # Timeline: 1880-01-16 ~ 2025-07-16
cat("Date range:", format(min(dates)), "to", format(max(dates)), "\n")

# --- 缺值處理 + 縮放 ---
anom_raw[anom_raw == fill] <- NA
anom <- anom_raw * ifelse(sf$hasatt, sf$value, 1)
if (ao$hasatt) anom <- anom + ao$value   # 單位 K == °C

# --- Grid of Taiwan (3*3) ---
lon360  <- ifelse(lon < 0, lon + 360, lon)
idx_lon <- which(lon360 >= 119 & lon360 <= 123)  # 119–123E
idx_lat <- which(lat    >=  21 & lat    <=  26)  # 21–26N

sub <- anom[idx_lon, idx_lat, , drop = FALSE]    # [x, y, t]

# 區域平均：只平均有值的格點
taiwan_mean <- apply(sub, 3, function(v) {
  m <- mean(v, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
})

# --- Mean of 5 yrs (60 months)  ---
df <- data.frame(Date = dates, Anomaly = taiwan_mean)
df$MA10 <- rollmean(df$Anomaly, k = 60, fill = NA, align = "center")

# --- Plot ---
p <- ggplot(df, aes(Date)) +
  geom_line(aes(y = Anomaly), color = "grey70", alpha = 0.7) +  # grey：mean of month
  geom_line(aes(y = MA10), color = "red", linewidth = 1) +      # red：mean of prev 10 yrs
  labs(
    title = "Taiwan temperature anomaly (baseline 1951–1980)",
    subtitle = paste(format(min(dates)), "to", format(max(dates)),
                     " | Box: 119–123E, 21–26N | GISTEMP Land+Ocean, 1200 km"),
    y = "°C anomaly", x = NULL
  ) +
  theme_minimal()
print(p)
