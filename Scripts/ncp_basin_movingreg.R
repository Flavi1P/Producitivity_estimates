# ============================================================
#  Basin-wide NCP for the Iceland Basin from nitrate drawdown
#
#  Strategy requested:
#   1) aggregate the entire BGC-Argo dataset to the SMALLEST time bin
#      that still keeps >= 3 profiles per bin (bins below that are dropped);
#   2) apply a MOVING linear regression of depth-integrated nitrate over a
#      30-day window to get the drawdown slope -> NCP (C units, Redfield).
#
#  Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/ncp_basin_movingreg.R
# ============================================================

suppressMessages({
  library(tidyverse)
  library(lubridate)
  library(castr)   # integrate() over depth
  library(zoo)     # rollapply
})

# ---- params ----------------------------------------------------------------
LON <- c(-40, -10); LAT <- c(58, 65)        # Iceland Basin box
MIN_PROF <- 3                                # >= 3 profiles per kept bin
REG_HALFWIN <- 15                            # +/- 15 d  -> 30-day moving window
CN <- 6.625                                  # Redfield C:N
ZEU <- 40                                    # fallback euphotic depth (m)
BIN_CANDIDATES <- c(2, 3, 5, 7, 10, 12, 15)  # days
COVER_THRESH <- 0.90                         # frac of bins that must reach MIN_PROF

# ---- 1. load + spatial filter ----------------------------------------------
files <- c("Data/Processed/doxy_db_2015_2019_float_with_canyon.csv",
           "Data/Processed/doxy_db_2020_2023_float_with_canyon.csv",
           "Data/Processed/doxy_db_2024_2025_float_with_canyon.csv")
dat <- map_dfr(files, ~read_csv2(.x, show_col_types = FALSE)) |>
  filter(lon >= LON[1], lon <= LON[2], lat >= LAT[1], lat <= LAT[2],
         !is.na(date)) |>
  mutate(date = as.Date(date), zeu = ZEU)

prof <- distinct(dat, float_wmo, prof_number, date)
message(sprintf("ICB: %d profiles, %d floats, %s -> %s",
                nrow(prof), n_distinct(dat$float_wmo),
                min(prof$date), max(prof$date)))

# ---- 2. pick the smallest bin width that keeps >= MIN_PROF per bin ----------
pick_bin <- function(B) {
  prof |> mutate(bin = as.Date(cut(date, breaks = paste(B, "days")))) |>
    count(bin) |> summarise(frac = mean(n >= MIN_PROF)) |> pull(frac)
}
cover <- map_dbl(BIN_CANDIDATES, pick_bin)
BIN <- BIN_CANDIDATES[which(cover >= COVER_THRESH)[1]]
if (is.na(BIN)) BIN <- max(BIN_CANDIDATES)
message(sprintf("chosen bin width = %d days (%.0f%% of bins reach >=%d profiles)",
                BIN, 100 * cover[match(BIN, BIN_CANDIDATES)], MIN_PROF))

binit <- function(d) as.Date(cut(d, breaks = paste(BIN, "days")))

# ---- 3. bins to keep (>= MIN_PROF distinct profiles) ------------------------
bin_keep <- prof |>
  mutate(bin = binit(date)) |>
  group_by(bin) |>
  summarise(n_prof = n_distinct(paste(float_wmo, prof_number)),
            date_bin = mean(date), .groups = "drop") |>
  filter(n_prof >= MIN_PROF)

# ---- 4. per-bin MLD and integration depth ----------------------------------
mld_bin <- dat |>
  distinct(float_wmo, prof_number, date, MLD) |>
  filter(!is.na(MLD)) |>
  mutate(bin = binit(date)) |>
  group_by(bin) |>
  summarise(mld = mean(MLD), .groups = "drop")

depth_bin <- bin_keep |>
  left_join(mld_bin, by = "bin") |>
  arrange(date_bin) |>
  mutate(
    zeu = ZEU,
    prev_mld = lag(mld), prev_zeu = lag(zeu),
    int_depth      = pmax(mld, prev_mld, zeu, prev_zeu, na.rm = TRUE)
  )

# ---- 5. per-bin mean nitrate profile, integrate over depth ------------------
nit_bin <- dat |>
  filter(!is.na(canyon_nitrate), depth <= 250) |>
  mutate(bin = binit(date)) |>
  semi_join(bin_keep, by = "bin") |>
  group_by(bin, depth) |>
  summarise(nitrate = mean(canyon_nitrate), .groups = "drop")

intN <- nit_bin |>
  left_join(select(depth_bin, bin, int_depth), by = "bin") |>
  group_by(bin, int_depth) |>
  summarise(
    int_N        = integrate(nitrate, depth, from = 0, to = unique(int_depth)),
    mld_conc     = integrate(nitrate, depth,
                             from = pmax(unique(int_depth) - 20, 0),
                             to   = pmax(unique(int_depth) - 10, 1)) / 10,
    sub_mld_conc = integrate(nitrate, depth,
                             from = unique(int_depth) + 20,
                             to   = unique(int_depth) + 30) / 10,
    .groups = "drop"
  )

# ---- 6. assemble time series -----------------------------------------------
ts <- depth_bin |>
  left_join(intN, by = c("bin", "int_depth")) |>
  arrange(date_bin) |>
  filter(!is.na(int_N)) |>
  mutate(
    int_N_smooth = rollapply(int_N, width = 3, FUN = mean,
                             align = "center", fill = NA),
    dt = as.numeric(date_bin - lag(date_bin))
  )

# ---- 7. moving 30-day regression of integrated nitrate ----------------------
# for each bin, slope of int_N_smooth ~ time over bins within +/- REG_HALFWIN d
ts$nitrate_slope <- map_dbl(seq_len(nrow(ts)), function(i) {
  d0 <- ts$date_bin[i]
  w <- ts |> filter(date_bin >= d0 - REG_HALFWIN, date_bin <= d0 + REG_HALFWIN,
                    !is.na(int_N_smooth))
  if (nrow(w) < MIN_PROF) return(NA_real_)
  coef(lm(int_N_smooth ~ as.numeric(date_bin), data = w))[2]   # mmol N m-2 d-1
})

# ---- 8. NCP (C units) + entrainment ----------------------------------------
ncp <- ts |>
  mutate(
    ncp_drawdown = -nitrate_slope * CN,                 # mmol C m-2 d-1
    we           = pmax(0, mld - lag(mld)),             # MLD deepening (m)
    delta        = sub_mld_conc - mld_conc,             # N jump across MLD
    entrain      = ifelse(dt > 0, delta * we * CN / dt, 0),
    NCP          = ncp_drawdown + entrain
  ) |>
  mutate(
    ncp_smooth = rollapply(NCP, width = 5, FUN = mean, align = "center", fill = NA),
    ncp_sd     = rollapply(NCP, width = 5, FUN = sd,   align = "center", fill = NA)
  ) |>
  filter(!is.na(NCP))

write_csv(select(ncp, date_bin, n_prof, int_depth, int_N, nitrate_slope,
                 ncp_drawdown, entrain, NCP, ncp_smooth, ncp_sd),
          "Data/Processed/ncp_basin_icb_movingreg.csv")

# ---- 9. NPP overlay (basin product, where available) -----------------------
npp <- read_csv("Data/Processed/icb_npp_ncp.csv", show_col_types = FALSE) |>
  transmute(date = as.Date(date_10day), npp = npp / 12) |>
  filter(!is.na(npp))

# ---- 10. plot --------------------------------------------------------------
ggplot(ncp) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(date_bin, ymin = ncp_smooth - ncp_sd,
                  ymax = ncp_smooth + ncp_sd), fill = "#08519c", alpha = 0.15) +
  geom_line(aes(date_bin, NCP, color = "NCP (per bin)"), alpha = 0.4) +
  geom_line(aes(date_bin, ncp_smooth, color = "NCP (smoothed)"), linewidth = 1) +
  geom_line(aes(date, npp, color = "NPP (basin)"), data = npp, linewidth = 0.7) +
  scale_color_manual(name = "",
    values = c("NCP (per bin)" = "#6baed6", "NCP (smoothed)" = "#08519c",
               "NPP (basin)" = "#252525")) +
  scale_x_date(date_labels = "%Y", breaks = "1 year") +
  labs(x = "Date", y = expression(mmol~C~m^-2~d^-1),
       title = sprintf("Iceland Basin NCP from nitrate drawdown (%d-day bins, 30-day moving regression)", BIN),
       subtitle = sprintf("%d floats, %d profiles, %d bins kept (>=%d profiles each)",
                          n_distinct(dat$float_wmo), nrow(prof), nrow(ncp), MIN_PROF)) +
  theme_bw()

ggsave("Output/ncp_basin_icb_movingreg.png", width = 13, height = 6, dpi = 200)
message("saved -> Output/ncp_basin_icb_movingreg.png")
