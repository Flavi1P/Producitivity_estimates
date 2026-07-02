# o2_budget_3902681_movingmld.R
# -----------------------------------------------------------------------------
# Diel O2 budget for BGC-Argo float 3902681, integrated to a MOVING bottom depth
# that tracks the mixed-layer depth (MLD), clamped to [50, 200] m, instead of a
# single fixed slab (parent: 0-40 m) or the fixed 40/100/150/200 m candidates in
# o2_budget_3902681_multidepth.R.
#
#   z_bot(pair) = min(200, max(50, max(MLD_i, MLD_j)))
#
# Rationale: the fixed-200 m figure (Scripts/gop_vs_ncp_figure.py) was chosen
# because 200 m is deep enough to capture the O2 inventory change in phase with
# the nitrate-budget NCP during the deep-mixing bloom onset, at the cost of a
# much noisier smoothed line in summer (integrating far below the productive
# layer). A moving MLD bottom keeps the integration matched to the layer that is
# actually ventilating/producing: deep in winter/spring (up to the 200 m cap),
# shallow in the stratified summer (down to the 50 m floor). The 50 m floor
# avoids integrating over an unrealistically thin summer ML where the diel O2
# signal is small and the areal budget would collapse; the 200 m cap matches the
# fixed-depth figure's deepest candidate (and the depth most casts reliably
# reach).
#
# The integration bottom is COMMON to the two casts in a pair (max of the two
# MLDs, clamped) so the inventory difference inv_j - inv_i is a true change of a
# fixed control volume, not a differencing of inventories to two different
# depths. This mirrors the NCP convention (integration depth = the deeper MLD
# across the differenced pair) in Scripts/ncp_function.R and the parent MLD
# budget in Scripts/o2_budget_3902681.R (which floors at 40 m, no upper cap).
#
# Everything else -- profile ingest, QC, surface state, diffusive-only
# Wanninkhof 2014 air-sea flux, phase-aware dusk/dawn pairing, edge trimming --
# is identical to o2_budget_3902681_multidepth.R. MLD is computed per profile
# with castr::mld() (sigma0 0.03 kg m-3 criterion, ref 0-2 m), matching the
# parent script and format_for_ncp.R.
#
# Output: Data/Processed/o2_budget_3902681_movingmld.csv
#   columns: type (net/night), mtime, dt_days, z_bot_m, mld_i_m, mld_j_m,
#            rate_mmol_o2_m2_d, Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d
#            (rate_corr = rate - Fas, diffusive-only correction)
#
# Run from repo root with R:
#   & "C:/Users/flapet/AppData/Local/Programs/R/R-4.3.3/bin/x64/Rscript.exe" Scripts/o2_budget_3902681_movingmld.R
# -----------------------------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(lubridate)
library(gsw)
library(castr)     # mld() - same MLD convention as the parent script / format_for_ncp.R

# ---- paths ------------------------------------------------------------------
nc_path   <- "Data/Raw/Floats/3902681_Sprof.nc"
era5_path <- "Data/Raw/ERA5/era5_wind_slp_3902681.nc"
out_csv   <- "Data/Processed/o2_budget_3902681_movingmld.csv"

MLD_FLOOR <- 50    # m, shallowest integration bottom (avoid a too-thin summer ML)
MLD_CAP   <- 200   # m, deepest integration bottom (matches the fixed-depth figure)

# =============================================================================
# Step 1 - per-profile full O2 profile + surface state + MLD
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

e_nc  <- nc_open(era5_path)
e_lon <- ncvar_get(e_nc, "longitude")
e_lat <- ncvar_get(e_nc, "latitude")
e_t   <- as.POSIXct(ncvar_get(e_nc, "time") * 3600, origin = "1900-01-01", tz = "UTC")
e_u10 <- ncvar_get(e_nc, "u10")
e_v10 <- ncvar_get(e_nc, "v10")
e_msl <- ncvar_get(e_nc, "msl")
nc_close(e_nc)
stopifnot(identical(dim(e_u10), c(length(e_lon), length(e_lat), length(e_t))))
message(sprintf("ERA5: %d hourly steps, %s -> %s",
                length(e_t), format(min(e_t)), format(max(e_t))))

era5_at <- function(t0, t1, lon0, lat0) {
  i  <- which.min(abs(e_lon - lon0))
  j  <- which.min(abs(e_lat - lat0))
  kt <- which(e_t >= t0 & e_t <= t1)
  if (length(kt) == 0) kt <- which.min(abs(as.numeric(e_t - (t0 + (t1 - t0) / 2))))
  u <- e_u10[i, j, kt]; v <- e_v10[i, j, kt]; m <- e_msl[i, j, kt]
  list(u2 = mean(u^2 + v^2, na.rm = TRUE), msl = mean(m, na.rm = TRUE))
}

n_prof <- length(juld)
time_prof <- as.POSIXct(juld * 86400, origin = "1950-01-01", tz = "UTC")
good_qc <- c(1L, 2L, 5L, 8L)

parse_qc_col <- function(qc_entry, n_levels) {
  chars <- strsplit(qc_entry, "")[[1]]
  flags <- suppressWarnings(as.integer(chars))
  length(flags) <- n_levels
  flags
}

o2_surf_vec <- rep(NA_real_, n_prof)
o2_eq_vec   <- rep(NA_real_, n_prof)
sst_vec     <- rep(NA_real_, n_prof)
sss_vec     <- rep(NA_real_, n_prof)
rho_surf_vec <- rep(NA_real_, n_prof)
mld_vec     <- rep(NA_real_, n_prof)   # mixed-layer depth, m (castr, sigma0 0.03)
o2_prof_list <- vector("list", n_prof) # depth-ordered full O2 profile, per prof_index
max_depth_vec <- rep(NA_real_, n_prof)

for (p in seq_len(n_prof)) {
  pr <- pres[, p]; tt <- temp[, p]; ss <- psal[, p]; oo <- doxy[, p]
  qc <- parse_qc_col(doxy_qc[p], length(pr))
  keep <- is.finite(pr) & is.finite(tt) & is.finite(ss) & is.finite(oo) &
    (qc %in% good_qc)
  if (sum(keep) < 5) next
  pr <- pr[keep]; tt <- tt[keep]; ss <- ss[keep]; oo <- oo[keep]

  SA  <- gsw_SA_from_SP(ss, pr, lon[p], lat[p])
  CT  <- gsw_CT_from_t(SA, tt, pr)
  rho <- gsw_rho(SA, CT, pr)
  o2_vol <- oo * rho / 1000
  depth  <- gsw_z_from_p(pr, lat[p]) * -1

  ord <- order(depth)
  depth <- depth[ord]; o2_vol <- o2_vol[ord]
  tt <- tt[ord]; ss <- ss[ord]; SA <- SA[ord]; CT <- CT[ord]; rho <- rho[ord]

  if (max(depth, na.rm = TRUE) < 40) next
  if (sum(depth <= 60) < 5) next

  max_depth_vec[p] <- max(depth, na.rm = TRUE)
  o2_prof_list[[p]] <- list(depth = depth, o2_vol = o2_vol)

  # ---- mixed-layer depth (sigma0 0.03 criterion, castr) ----------------------
  sigma0 <- gsw_sigma0(SA, CT)
  mld_vec[p] <- castr::mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03,
                           default.depth = max(depth, na.rm = TRUE))

  top <- depth <= 10
  if (sum(top) == 0) top <- seq_len(min(3, length(depth)))
  SA_s  <- mean(SA[top]); CT_s <- mean(CT[top]); rho_s <- mean(rho[top])
  sst_vec[p]      <- mean(tt[top])
  sss_vec[p]      <- mean(ss[top])
  rho_surf_vec[p] <- rho_s
  o2_grid10 <- approx(x = depth, y = o2_vol, xout = 0:10, rule = 2)$y
  o2_surf_vec[p] <- mean(o2_grid10)
  o2_eq_umolkg   <- gsw_O2sol(SA_s, CT_s, 0, lon[p], lat[p])
  o2_eq_vec[p]   <- o2_eq_umolkg * rho_s / 1000
}

prof <- tibble(
  prof_index = seq_len(n_prof),
  time = time_prof, lat = lat, lon = lon,
  o2_surf = o2_surf_vec, o2_eq = o2_eq_vec,
  sst = sst_vec, sss = sss_vec, rho_surf = rho_surf_vec,
  mld = mld_vec, max_depth = max_depth_vec
) |>
  mutate(
    lst   = (hour(time) + minute(time) / 60 + lon / 15) %% 24,
    phase = ifelse(lst < 12, "dawn", "dusk")
  ) |>
  filter(!is.na(o2_surf)) |>
  arrange(time)

# same edge-trimming as the parent / multidepth script (irregular near-midday casts)
n_drop <- 5
stopifnot(nrow(prof) > 2 * n_drop)
prof <- prof |> slice((n_drop + 1):(n() - n_drop))

message(sprintf("Profiles retained: %d / %d (after trimming %d at each end)",
                nrow(prof), n_prof, n_drop))
message(sprintf("MLD: median=%.0f m, IQR=[%.0f, %.0f] m; clamped bottom range [%d, %d] m",
                median(prof$mld, na.rm = TRUE),
                quantile(prof$mld, 0.25, na.rm = TRUE),
                quantile(prof$mld, 0.75, na.rm = TRUE),
                MLD_FLOOR, MLD_CAP))

# =============================================================================
# Step 2 - helpers: fixed-depth integral, diffusive-only air-sea flux, pairing
# =============================================================================
integrate_o2 <- function(pi, z_bot) {
  pr <- o2_prof_list[[pi]]
  if (is.null(pr)) return(NA_real_)
  if (max(pr$depth, na.rm = TRUE) < z_bot) return(NA_real_)
  g <- seq(0, z_bot, by = 1)
  if (g[length(g)] < z_bot) g <- c(g, z_bot)
  yv <- approx(x = pr$depth, y = pr$o2_vol, xout = g, rule = 2)$y
  if (any(is.na(yv))) return(NA_real_)
  sum(diff(g) * (head(yv, -1) + tail(yv, -1)) / 2)
}

# diffusive-only air-sea O2 flux (Wanninkhof 2014); surface boundary term,
# independent of the column integration depth. Identical to the parent /
# multidepth script.
airsea_flux_diffusive <- function(d, i, j) {
  T_s <- mean(c(d$sst[i], d$sst[j]))
  Sc  <- 1920.4 - 135.6 * T_s + 5.2122 * T_s^2 - 0.10939 * T_s^3 +
         0.00093777 * T_s^4
  lon_mid <- mean(c(d$lon[i], d$lon[j]))
  lat_mid <- mean(c(d$lat[i], d$lat[j]))
  w   <- era5_at(d$time[i], d$time[j], lon_mid, lat_mid)
  k   <- 0.251 * w$u2 * (Sc / 660)^(-0.5) * 0.24            # m d-1

  o2_eq_p_i <- d$o2_eq[i] * (w$msl / 101325)
  o2_eq_p_j <- d$o2_eq[j] * (w$msl / 101325)
  dO2 <- mean(c(d$o2_surf[i], d$o2_surf[j])) - mean(c(o2_eq_p_i, o2_eq_p_j))
  Fas <- -k * dO2                                           # mmol O2 m-2 d-1, +into ocean
  Fas
}

# Moving common bottom for the pair = max of the two MLDs, clamped to [floor, cap].
pair_z_bot <- function(d, i, j) {
  z <- max(d$mld[i], d$mld[j], na.rm = TRUE)
  min(MLD_CAP, max(MLD_FLOOR, z))
}

make_pair_moving <- function(d, i, j, type) {
  dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))
  z_bot <- pair_z_bot(d, i, j)
  inv_i <- integrate_o2(d$prof_index[i], z_bot)
  inv_j <- integrate_o2(d$prof_index[j], z_bot)
  rate  <- (inv_j - inv_i) / dt_days
  Fas   <- airsea_flux_diffusive(d, i, j)
  tibble(
    type = type,
    mtime = d$time[i] + (d$time[j] - d$time[i]) / 2,
    dt_days = dt_days,
    z_bot_m = z_bot,
    mld_i_m = d$mld[i], mld_j_m = d$mld[j],
    rate_mmol_o2_m2_d = rate,
    Fas_mmol_o2_m2_d = Fas,
    rate_corr_mmol_o2_m2_d = rate - Fas
  )
}

n <- nrow(prof)

# --- night pairs: dusk -> dawn (< 0.7 d apart), phase-aware (mirrors parent) ---
night_tbl <- map_dfr(seq_len(n - 1), function(k) {
  if (prof$phase[k] == "dusk" && prof$phase[k + 1] == "dawn") {
    dt <- as.numeric(difftime(prof$time[k + 1], prof$time[k], units = "days"))
    if (dt < 0.7) return(make_pair_moving(prof, k, k + 1, "night"))
  }
  tibble()
})

# --- net pairs: consecutive within the same phase (dusk-dusk, dawn-dawn) -------
net_tbl <- map_dfr(c("dusk", "dawn"), function(ph) {
  sub <- prof |> filter(phase == ph)
  if (nrow(sub) < 2) return(tibble())
  map_dfr(seq_len(nrow(sub) - 1), ~ make_pair_moving(sub, .x, .x + 1, "net"))
}) |>
  arrange(mtime)

budget_all <- bind_rows(night_tbl, net_tbl)

n_ok_net   <- sum(is.finite(net_tbl$rate_corr_mmol_o2_m2_d))
n_ok_night <- sum(is.finite(night_tbl$rate_corr_mmol_o2_m2_d))
message(sprintf(
  "moving MLD: net pairs %d/%d finite (%.0f%%), night pairs %d/%d finite (%.0f%%)",
  n_ok_net, nrow(net_tbl), 100 * n_ok_net / nrow(net_tbl),
  n_ok_night, nrow(night_tbl), 100 * n_ok_night / nrow(night_tbl)))

budget_out <- budget_all |>
  filter(is.finite(rate_corr_mmol_o2_m2_d)) |>
  select(type, mtime, dt_days, z_bot_m, mld_i_m, mld_j_m,
         rate_mmol_o2_m2_d, Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d) |>
  arrange(type, mtime)

message(sprintf("z_bot used: median=%.0f m, at floor(%d) %.0f%%, at cap(%d) %.0f%%",
                median(budget_out$z_bot_m),
                MLD_FLOOR, 100 * mean(budget_out$z_bot_m <= MLD_FLOOR),
                MLD_CAP,   100 * mean(budget_out$z_bot_m >= MLD_CAP)))

write_csv(budget_out, out_csv)
cat(sprintf("Wrote %s (%d rows)\n", out_csv, nrow(budget_out)))
