# === Library ===
library(terra)
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)

# === Working Dir ===
setwd("/Users/bill/Documents/NCKU/大二/R語言&黑客松/GlobalWarming")

# === Raw Data Dir ===
tw_extreme_rainfalldata_dir <- "Data_Set/tw_extreme_indices.csv"
temp_data_dir <- "Data_Set/gistemp1200_GHCNv4_ERSSTv5.nc"
exceed_csv  <- "Data_Set/tw_exceedance_counts.csv"
gistemp_nc  <- "Data_Set/gistemp1200_GHCNv4_ERSSTv5.nc"

# === Get Raw data ===
ext <- read_csv(tw_extreme_rainfalldata_dir) %>%
  arrange(Year) %>%
  mutate(
    Rx1day_roll10 = zoo::rollmean(Rx1day, 10, fill = NA, align = "right"),
    Rx5day_roll10 = zoo::rollmean(Rx5day, 10, fill = NA, align = "right")
  )

exc <- read_csv(exceed_csv) %>%
  arrange(Year) %>%
  mutate(
    days200_roll10 = zoo::rollmean(days_ge_200mm, 10, fill = NA, align = "right")
  )

# === The Time Range is from 1981 to 2024 ===
# Plot Picture 1 -> Rx1 day rainfall
ggplot(ext, aes(x = Year, y = Rx1day)) +
  geom_line(color = "grey50") +
  geom_point(size = 1, color = "grey50") +
  geom_line(aes(y = Rx1day_roll10), color = "red", linewidth = 1) +
  stat_smooth(method = "lm", se = FALSE, linetype = 2) +
  labs(x = "Year", y = "Most rainfall day annually : Rx1day (mm)",
       title = "Taiwan's annual maximum daily rainfall trend (including 10-year moving average)") +
  theme_minimal()

# Plot Picture 2 -> Rx5 day rainfall
ggplot(ext, aes(x = Year, y = Rx5day)) +
  geom_line(color = "grey50") +
  geom_point(size = 1, color = "grey50") +
  geom_line(aes(y = Rx5day_roll10), color = "red", linewidth = 1) +
  stat_smooth(method = "lm", se = FALSE, linetype = 2) +
  labs(x = "Year", y = "Taiwan's annual maximum 5-day rainfall : Rx5day (mm)",
       title = "Taiwan's annual maximum 5-day rainfall trend (including 10-year moving average)") +
  theme_minimal()

# === Plot Picture 3 -> days >=200mm ===
ggplot(exc, aes(x = Year, y = days_ge_200mm)) +
  geom_line(color = "grey50") +
  geom_point(size = 1, color = "grey50") +
  geom_line(aes(y = days200_roll10), color = "red", linewidth = 1) +
  stat_smooth(method = "lm", se = FALSE, linetype = 2) +
  labs(x = "Year",
       y = "Days with rainfall ≥200mm",
       title = "Taiwan's annual count of days ≥200mm (including 10-year moving average)") +
  theme_minimal()

# === Plot Picture4 -> extreme count - year temp ===

# Taiwan Box
lon_min <- 119; lon_max <- 123
lat_min <- 21;  lat_max <- 26

r <- rast(gistemp_nc)                                # 月層 NetCDF
r_tw <- crop(r, ext(lon_min, lon_max, lat_min, lat_max))
tvec <- time(r_tw); if (is.null(tvec)) stop("GISTEMP nc 缺少 time 軸。")

# Weight cos(lat)
lat_rast <- init(r_tw, "y")          # 每個 cell 的緯度
w_rast   <- cos(pi * lat_rast / 180) # 權重 raster

# >>> A) 修正：逐層計算加權平均，取純數值 <<<
monthly_vals <- sapply(1:nlyr(r_tw), function(i){
  m   <- r_tw[[i]]
  num <- global(m * w_rast, "sum", na.rm = TRUE)[1,1]
  den <- global(w_rast * !is.na(m), "sum", na.rm = TRUE)[1,1]
  num / den
})

# 月 -> 年
temp_annual <- tibble(
  Date = as.Date(time(r_tw)),
  Temp = as.numeric(monthly_vals)
) |>
  arrange(Date) |>
  mutate(Year = year(Date)) |>
  group_by(Year) |>
  summarise(Temp_Annual = mean(Temp, na.rm = TRUE), .groups = "drop")

# Read days >= 200 (mm)
exc <- read_csv(exceed_csv, show_col_types = FALSE)
names_lc <- tolower(names(exc))

if (!("year" %in% names_lc)) {
  dc <- names(exc)[grepl("date", names(exc), ignore.case = TRUE)][1]
  if (is.na(dc)) stop("exceedance 檔缺少 Year/Date 欄位。")
  exc <- exc |>
    mutate(Year = year(as.Date(.data[[dc]])))
} else {
  names(exc)[which(names_lc == "year")] <- "Year"
}

cand <- names(exc)[tolower(names(exc)) %in% c("days_ge_200mm","days_ge_200","ge_200mm_days")]
if (!length(cand)) stop("找不到 days_ge_200mm 欄位。現有欄位：", paste(names(exc), collapse=", "))

# >>> B) 修正：用 sym() 讓 rename 支援字串欄名 <<<
exc <- exc |>
  rename(RainVar = !!rlang::sym(cand[1])) |>
  mutate(RainVar = suppressWarnings(as.numeric(RainVar))) |>
  group_by(Year) |>
  summarise(RainVar = sum(RainVar, na.rm = TRUE), .groups = "drop")

# 合併
df_join <- exc |>
  inner_join(temp_annual, by = "Year") |>
  arrange(Year)

cat("對應年份：", min(df_join$Year), "~", max(df_join$Year), "\n")
cat("使用降雨欄位：days_ge_200mm\n")

# 標準化年序列
df_plot <- df_join |>
  mutate(
    Z_Rain = scale(RainVar)[,1],
    Z_Temp = scale(Temp_Annual)[,1]
  )

p_ts <- ggplot(df_plot, aes(x = Year)) +
  geom_line(aes(y = Z_Rain), color = "grey40") +
  geom_line(aes(y = Z_Temp), linetype = "dashed", color = "red") +
  labs(title = "Taiwan: Days ≥200mm vs Temperature (standardized)",
       subtitle = "solid = days ≥200mm (z), dashed = temperature (z)",
       x = "Year", y = "z-score") +
  theme_minimal(base_size = 12)
print(p_ts)

# 散點 + 回歸線 + r/p
lm_fit <- lm(RainVar ~ Temp_Annual, data = df_join)
r_val  <- cor(df_join$Temp_Annual, df_join$RainVar, use = "complete.obs")
p_val  <- summary(lm_fit)$coefficients[2,4]

p_scatter <- ggplot(df_join, aes(x = Temp_Annual, y = RainVar)) +
  geom_point(alpha = 0.8, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Days ≥200mm per Year vs Annual Temperature (GISTEMP, Taiwan)",
       subtitle = sprintf("Pearson r = %.2f,  p = %.3g", r_val, p_val),
       x = "Temperature (annual mean)",
       y = "Days ≥200mm") +
  theme_minimal(base_size = 12)
print(p_scatter)




