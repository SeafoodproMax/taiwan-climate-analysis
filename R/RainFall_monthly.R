library(readr)
library(dplyr)
library(stringr)
library(lubridate)
library(purrr)
library(tidyr)
library(ggplot2)
suppressWarnings({
  if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
})
library(zoo)  # mean of 10 years
library(here)

# The Dir of Data_set
data_dir <- "../Data_Set/RainFall"
# The Dir of wd
setwd(here("R"))

files <- list.files(data_dir, pattern = "^ObsRain_.*\\.csv$", full.names = TRUE)

# county name getter
extract_county <- function(path) {
  fn <- basename(path)
  name <- fn %>%
    str_remove("^ObsRain_") %>%
    str_remove("\\.csv$")
  name
}
# Function to Read one file
read_one <- function(path, encoding = "UTF-8", to_mm = FALSE) {
  df <- read_csv(path, locale = locale(encoding = encoding), show_col_types = FALSE)
  
  need <- c("CityName", "YY", "MM", "WGS84_Lon", "WGS84_Lat", "RainValue")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop(sprintf("缺少欄位：%s\n檔案：%s", paste(miss, collapse = ", "), path))
  
  date <- make_date(year = as.integer(df$YY), month = as.integer(df$MM), day = 1)
  precip <- suppressWarnings(as.numeric(df$RainValue))
  if (isTRUE(to_mm)) precip <- precip * 1000
  
  tibble(
    county = df$CityName,
    date   = date,
    year   = year(date),
    month  = month(date),
    lon    = as.numeric(df$WGS84_Lon),
    lat    = as.numeric(df$WGS84_Lat),
    precip = precip
  )
}
# Map 20 of the Data.csv files to a map
all_cells <- map_dfr(files, read_one)

time_head <- min(all_cells$date, na.rm = TRUE)
time_tail <- max(all_cells$date, na.rm = TRUE)
message(sprintf("時間範圍：%s ~ %s", format(time_head, "%Y-%m"), format(time_tail, "%Y-%m")))
cat(sprintf("時間範圍：%s ~ %s\n", format(time_head, "%Y-%m"), format(time_tail, "%Y-%m")))
cat("縣市數：", dplyr::n_distinct(all_cells$county), 
    "；總格網列數：", nrow(all_cells), "\n")

# Do average with 1871 grids in Taiwan for each month, from 1960-1~2022-12 756 months and 1.4 million data
taiwan_monthly_grid_eq <- all_cells |>
  dplyr::group_by(date) |>
  dplyr::summarise(
    n_cells       = dplyr::n(),
    mean_precip   = mean(precip, na.rm = TRUE),
    median_precip = median(precip, na.rm = TRUE),
    sum_precip    = sum(precip, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(date) |>
  dplyr::mutate(
    ma_10yr = zoo::rollmean(mean_precip, k = 120, fill = NA, align = "right")
  )
# Data Test
cat("格網等權：月份筆數 =", nrow(taiwan_monthly_grid_eq), "\n")

# Plot
p1 <- ggplot(taiwan_monthly_grid_eq, aes(x = date)) +
  geom_line(aes(y = mean_precip)) +
  geom_line(aes(y = ma_10yr), linewidth = 1) +
  labs(title = "台灣逐月降雨（格網等權空間平均）",
       subtitle = "深色線：10年移動平均",
       x = "時間", y = "降雨（原始單位）") +
  theme_minimal(base_size = 12)
print(p1)

