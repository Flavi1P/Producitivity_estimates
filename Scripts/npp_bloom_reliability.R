# Diagnostic map of per-cell bloom-season NPP variability.
#
# For each year (2015-2023), integrates monthly npp_int over Apr-Jul
# (sum * 30.4 days -> g C m-2 over the season).  Plots, side by side:
#   - mean bloom-season integrated NPP (climatological intensity)
#   - CV across years (sd / mean) -- the "reliability" metric
#
# High CV = a cell whose bloom magnitude swings widely year to year;
# may be a "doesn't always bloom" cell.
#
# Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/npp_bloom_reliability.R

suppressPackageStartupMessages({
  library(ncdf4)
  library(tidyverse)
  library(RColorBrewer)
  library(patchwork)
})

NC_FILE <- "Output/cmems_cbpm_monthly_0p25deg.nc"
OUT_DIR <- "Output"
BLOOM_MONTHS <- 4:7
DAYS_PER_MONTH <- 30.4

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading ", NC_FILE)
nc <- nc_open(NC_FILE)
lons <- as.numeric(ncvar_get(nc, "longitude"))
lats <- as.numeric(ncvar_get(nc, "latitude"))
t_units <- ncatt_get(nc, "time", "units")$value
t_raw   <- ncvar_get(nc, "time")
unit_word  <- sub(" since.*$", "", t_units)
origin_str <- sub("^[A-Za-z]+ since ", "", t_units)
origin     <- as.Date(substr(origin_str, 1, 10))
times <- switch(unit_word,
                "days"    = origin + t_raw,
                "hours"   = origin + t_raw / 24,
                "seconds" = origin + t_raw / 86400,
                stop("Unsupported time unit: ", unit_word))
npp_int <- ncvar_get(nc, "npp_int")    # (lon, lat, time)
nc_close(nc)

npp_int <- aperm(npp_int, c(3, 2, 1))  # (time, lat, lon)
months  <- as.numeric(format(times, "%m"))
years   <- as.numeric(format(times, "%Y"))
unique_years <- sort(unique(years))

# Per-year, per-cell bloom-season integral (g C m-2 over Apr-Jul)
bloom_mat <- array(NA_real_, dim = c(length(unique_years),
                                     length(lats), length(lons)))
for (yi in seq_along(unique_years)) {
  yr <- unique_years[yi]
  t_mask <- (years == yr) & (months %in% BLOOM_MONTHS)
  if (sum(t_mask) == 0) next
  monthly <- npp_int[t_mask, , , drop = FALSE]   # (n_months, lat, lon)
  # integral = sum over months * days_per_month, in mg C m-2 -> /1000 = g C m-2
  bloom_mat[yi, , ] <- apply(monthly, c(2, 3),
                             function(v) sum(v, na.rm = TRUE) *
                                         DAYS_PER_MONTH / 1000)
}

# Mean and sd across years
bloom_mean <- apply(bloom_mat, c(2, 3), mean, na.rm = TRUE)
bloom_sd   <- apply(bloom_mat, c(2, 3), sd,   na.rm = TRUE)
bloom_cv   <- bloom_sd / bloom_mean
bloom_cv[!is.finite(bloom_cv) | bloom_mean < 1] <- NA_real_  # mask near-zero mean

# Flatten to long form for ggplot
df <- expand.grid(lat = lats, lon = lons) |>
  as_tibble() |>
  mutate(mean = as.vector(bloom_mean),
         cv   = as.vector(bloom_cv)) |>
  filter(is.finite(mean))

message("Cells with finite bloom mean: ", nrow(df))
message("Bloom-mean: median = ", sprintf("%.1f g C m-2", median(df$mean)),
        "  range = [", sprintf("%.1f", min(df$mean)), ", ",
        sprintf("%.1f", max(df$mean)), "]")
message("Bloom-CV  : median = ", sprintf("%.2f", median(df$cv, na.rm = TRUE)),
        "  P10-P90 = [", sprintf("%.2f", quantile(df$cv, 0.10, na.rm = TRUE)),
        ", ", sprintf("%.2f", quantile(df$cv, 0.90, na.rm = TRUE)), "]")

p_mean <- ggplot(df, aes(lon, lat, fill = mean)) +
  geom_tile(width = 0.25, height = 0.25) +
  borders("world", xlim = c(-45, -5), ylim = c(56, 67),
          fill = "grey90", colour = "black", linewidth = 0.3) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1,
                       name = expression("g C m"^-2)) +
  coord_quickmap(xlim = c(-41, -9), ylim = c(57, 66)) +
  labs(x = "Longitude", y = "Latitude",
       title = "Bloom-season mean intensity",
       subtitle = "Apr-Jul integrated NPP, averaged over 2015-2023") +
  theme_bw()

p_cv <- ggplot(df, aes(lon, lat, fill = cv)) +
  geom_tile(width = 0.25, height = 0.25) +
  borders("world", xlim = c(-45, -5), ylim = c(56, 67),
          fill = "grey90", colour = "black", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "CV", na.value = "grey80") +
  coord_quickmap(xlim = c(-41, -9), ylim = c(57, 66)) +
  labs(x = "Longitude", y = "Latitude",
       title = "Bloom-season interannual variability",
       subtitle = "CV (sd / mean) of Apr-Jul integrated NPP across 9 years",
       caption = "High CV = unreliable bloom (large year-to-year swings)") +
  theme_bw()

out_png <- file.path(OUT_DIR, "npp_bloom_reliability.png")
ggsave(out_png, p_mean / p_cv, width = 9, height = 9, dpi = 140)
message("Wrote ", out_png)

# Also persist the per-cell metrics as CSV for downstream use
out_csv <- file.path(OUT_DIR, "npp_bloom_reliability.csv")
write_csv(df |> dplyr::select(lon, lat, bloom_mean = mean, bloom_cv = cv),
          out_csv)
message("Wrote ", out_csv)
