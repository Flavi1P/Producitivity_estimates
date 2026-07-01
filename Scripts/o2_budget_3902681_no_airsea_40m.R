# o2_budget_3902681_no_airsea_40m.R
# -----------------------------------------------------------------------------
# Sanity check: can the raw (NOT flux-corrected) 0-40 m diel O2 budget for float
# 3902681 reproduce the reference "net change" / "night loss" series in
# Data/Processed/O2_float_net_change.xlsx and O2_float_night_loss.xlsx?
#
# Same profile ingest and phase-pairing as Scripts/o2_budget_3902681.R, but
# NO air-sea flux term at all (no ERA5 needed): rate = raw dO2/dt of the fixed
# 0-40 m inventory. This isolates whether the reference series are themselves
# uncorrected 0-40 m budgets.
#
# Output:
#   Output/o2_budget_3902681_no_airsea_40m_comparison.png
#
# Run from repo root with R 4.5.2:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681_no_airsea_40m.R
# -----------------------------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(lubridate)
library(gsw)
library(zoo)
library(readxl)

# ---- paths ------------------------------------------------------------------
nc_path        <- "Data/Raw/Floats/3902681_Sprof.nc"
ref_net_path   <- "Data/Processed/O2_float_net_change.xlsx"
ref_night_path <- "Data/Processed/O2_float_night_loss.xlsx"
out_png        <- "Output/o2_budget_3902681_no_airsea_40m_comparison.png"

# =============================================================================
# Step 1 - per-profile O2 inventory (0-40 m), no surface/flux state needed
# =============================================================================
nc <- nc_open(nc_path)

lon  <- ncvar_get(nc, "LONGITUDE")
lat  <- ncvar_get(nc, "LATITUDE")
juld <- ncvar_get(nc, "JULD")
pres <- ncvar_get(nc, "PRES")
temp <- ncvar_get(nc, "TEMP")
psal <- ncvar_get(nc, "PSAL")
doxy <- ncvar_get(nc, "DOXY_ADJUSTED")
doxy_qc <- ncvar_get(nc, "DOXY_ADJUSTED_QC")

nc_close(nc)

n_prof <- length(juld)
time_prof <- as.POSIXct(juld * 86400, origin = "1950-01-01", tz = "UTC")

good_qc <- c(1L, 2L, 5L, 8L)
grid <- 0:40

parse_qc_col <- function(qc_entry, n_levels) {
  chars <- strsplit(qc_entry, "")[[1]]
  flags <- suppressWarnings(as.integer(chars))
  length(flags) <- n_levels
  flags
}

inventory_vec <- rep(NA_real_, n_prof)

for (p in seq_len(n_prof)) {
  pr <- pres[, p]
  tt <- temp[, p]
  ss <- psal[, p]
  oo <- doxy[, p]
  qc <- parse_qc_col(doxy_qc[p], length(pr))

  keep <- is.finite(pr) & is.finite(tt) & is.finite(ss) & is.finite(oo) &
    (qc %in% good_qc)
  if (sum(keep) < 5) next

  pr <- pr[keep]; tt <- tt[keep]; ss <- ss[keep]; oo <- oo[keep]

  SA  <- gsw_SA_from_SP(ss, pr, lon[p], lat[p])
  CT  <- gsw_CT_from_t(SA, tt, pr)
  rho <- gsw_rho(SA, CT, pr)            # kg m-3
  o2_vol <- oo * rho / 1000            # mmol O2 m-3

  depth <- gsw_z_from_p(pr, lat[p]) * -1
  ord <- order(depth)
  depth <- depth[ord]; o2_vol <- o2_vol[ord]

  if (max(depth, na.rm = TRUE) < 40) next
  if (sum(depth <= 60) < 5) next

  o2_grid <- approx(x = depth, y = o2_vol, xout = grid, rule = 2)$y
  if (any(is.na(o2_grid))) next

  inventory_vec[p] <- sum(diff(grid) * (head(o2_grid, -1) + tail(o2_grid, -1)) / 2)
}

prof <- tibble(
  prof_index = seq_len(n_prof),
  time       = time_prof,
  lat        = lat,
  lon        = lon,
  inventory  = inventory_vec
) |>
  mutate(
    lst   = (hour(time) + minute(time) / 60 + lon / 15) %% 24,
    phase = ifelse(lst < 12, "dawn", "dusk")
  ) |>
  filter(!is.na(inventory)) |>
  arrange(time)

# same edge-trimming as the parent script (irregular near-midday casts)
n_drop <- 5
stopifnot(nrow(prof) > 2 * n_drop)
prof <- prof |> slice((n_drop + 1):(n() - n_drop))

message(sprintf("Profiles with valid 0-40 m inventory: %d / %d (after trimming %d at each end)",
                nrow(prof), n_prof, n_drop))

# =============================================================================
# Step 2 - build net/night pairs, RAW rate only (no air-sea flux term)
# =============================================================================
make_pair <- function(d, i, j, type) {
  dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))
  rate    <- (d$inventory[j] - d$inventory[i]) / dt_days   # signed dO2/dt
  tibble(
    type    = type,
    mtime   = d$time[i] + (d$time[j] - d$time[i]) / 2,
    dt_days = dt_days,
    phase_start = d$phase[i],
    phase_end   = d$phase[j],
    rate_mmol_o2_m2_d = rate
  )
}

n <- nrow(prof)

prev_tbl <- map_dfr(seq_len(n - 1), ~ make_pair(prof, .x, .x + 1, "prev"))

night_tbl <- prev_tbl |>
  filter(phase_start == "dusk", phase_end == "dawn", dt_days < 0.7) |>
  mutate(type = "night")

net_tbl <- map_dfr(c("dusk", "dawn"), function(ph) {
  sub <- prof |> filter(phase == ph)
  if (nrow(sub) < 2) return(tibble())
  map_dfr(seq_len(nrow(sub) - 1), ~ make_pair(sub, .x, .x + 1, "net"))
}) |>
  arrange(mtime)

# =============================================================================
# Step 3 - reference series (xlsx, Matlab datenum -> POSIXct)
# =============================================================================
matlab_to_posix <- function(mtime) {
  as.POSIXct((mtime - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
}

# readxl reads these columns as character on this file (mixed formatting) --
# coerce to numeric before converting the Matlab datenum to a timestamp.
ref_net <- read_excel(ref_net_path) |>
  mutate(mtime = as.numeric(mtime),
         Int_dO2dt_40m_mmol_m2_d = as.numeric(Int_dO2dt_40m_mmol_m2_d),
         time = matlab_to_posix(mtime)) |>
  transmute(time, value = Int_dO2dt_40m_mmol_m2_d, series = "net")

# Reference night loss is a positive loss magnitude; negate our dO2/dt (which is
# negative overnight) to match that convention wherever we compare/plot night.
ref_night <- read_excel(ref_night_path) |>
  mutate(mtime = as.numeric(mtime),
         Int_O2_loss_40m_mmol_m2_d = as.numeric(Int_O2_loss_40m_mmol_m2_d),
         time = matlab_to_posix(mtime)) |>
  transmute(time, value = Int_O2_loss_40m_mmol_m2_d, series = "night")

# =============================================================================
# Step 4 - verification: correlation + bias vs reference (nearest within 1 day)
# =============================================================================
nearest_join <- function(mine, ref, tol_days = 1) {
  ref <- ref |> arrange(time)
  out <- mine |>
    arrange(time) |>
    rowwise() |>
    mutate(
      .di      = which.min(abs(as.numeric(difftime(time, ref$time, units = "days")))),
      ref_time = ref$time[.di],
      ref_val  = ref$value[.di],
      .gap     = abs(as.numeric(difftime(time, ref_time, units = "days")))
    ) |>
    ungroup() |>
    filter(.gap <= tol_days)
  out
}

report_fit <- function(label, mine, ref) {
  m <- mine |> filter(is.finite(rate_mmol_o2_m2_d)) |>
    transmute(time = mtime, value = rate_mmol_o2_m2_d)
  j <- nearest_join(m, ref)
  if (nrow(j) < 3) {
    cat(sprintf("%s: too few matches (%d)\n", label, nrow(j)))
    return(invisible())
  }
  r     <- cor(j$value, j$ref_val, use = "complete.obs")
  bias  <- median(j$value - j$ref_val, na.rm = TRUE)
  slope <- coef(lm(value ~ 0 + ref_val, data = j))[1]
  cat(sprintf("%s: n=%d  Pearson r=%.3f  median(mine-ref)=%.3f  slope(mine~0+ref)=%.3f\n",
              label, nrow(j), r, bias, slope))
}

cat("\n--- Verification (0-40 m, NO air-sea flux correction) ---\n")
cat(sprintf("night pairs: %d (ref rows: %d)\n", nrow(night_tbl), nrow(ref_night)))
cat(sprintf("net pairs:   %d (ref rows: %d)\n", nrow(net_tbl), nrow(ref_net)))
report_fit("NET   (mine vs ref net change)     ", net_tbl, ref_net)
report_fit("NIGHT (mine-loss vs ref loss)      ",
           night_tbl |> mutate(rate_mmol_o2_m2_d = -rate_mmol_o2_m2_d),
           ref_night)

# =============================================================================
# Step 5 - comparison plot
# =============================================================================
roll_smooth <- function(x, k = 8) {
  if (length(x) < k) return(x)
  zoo::rollapply(x, width = k, FUN = mean, na.rm = TRUE, fill = NA, align = "center")
}

mine_net <- net_tbl |>
  transmute(time = mtime, value = rate_mmol_o2_m2_d,
            quantity = "net change", source = "computed (no flux-corr)")
mine_night <- night_tbl |>
  transmute(time = mtime, value = -rate_mmol_o2_m2_d,
            quantity = "night loss", source = "computed (no flux-corr)")
ref_net_pl <- ref_net |>
  transmute(time, value, quantity = "net change", source = "reference")
ref_night_pl <- ref_night |>
  transmute(time, value, quantity = "night loss", source = "reference")

plot_df <- bind_rows(mine_net, mine_night, ref_net_pl, ref_night_pl) |>
  filter(is.finite(value)) |>
  arrange(quantity, source, time) |>
  group_by(quantity, source) |>
  mutate(value_smooth = roll_smooth(value, 8)) |>
  ungroup() |>
  mutate(source = factor(source, levels = c("computed (no flux-corr)", "reference")))

p <- ggplot(plot_df, aes(x = time, colour = quantity, linetype = source)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(aes(y = value_smooth), linewidth = 0.9) +
  scale_colour_manual(values = c("net change" = "#3690c0",
                                 "night loss" = "#cb181d")) +
  scale_linetype_manual(values = c("computed (no flux-corr)" = "solid",
                                   "reference" = "dashed")) +
  coord_cartesian(ylim = c(-200, 400)) +
  labs(
    title = "Float 3902681: O2 budget (0-40 m) - computed vs reference (NO air-sea flux correction)",
    subtitle = "Raw dO2/dt of the fixed 0-40 m inventory (no gas-exchange term) vs reference. Reference sign convention (night = loss, positive). Lines: 8-pt rolling mean.",
    x = "Date", y = expression("mmol O"[2] ~ "m"^-2 ~ "d"^-1),
    colour = "Quantity", linetype = "Source"
  ) +
  theme_bw()

ggsave(out_png, p, width = 11, height = 6, dpi = 200)
cat(sprintf("Wrote %s\n", out_png))
