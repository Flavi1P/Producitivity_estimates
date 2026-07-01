# o2_budget_3902681_multidepth.R
# -----------------------------------------------------------------------------
# Diel O2 budget for BGC-Argo float 3902681, integrated to several FIXED
# candidate bottom depths (40, 100, 150, 200 m) instead of the parent script's
# single fixed 0-40 m slab.
#
# Reuses the exact same profile ingest, surface state, and air-sea flux
# (Wanninkhof 2014 diffusive correction -- no bubble term, matching the paper
# convention in mat_and_meth.md S2.2.1) as Scripts/o2_budget_3902681.R. Only
# the integration bottom differs: each profile's full depth-ordered O2 profile
# is re-integrated (trapezoidal, 1 m grid) from 0 to z_bot for z_bot in
# CANDIDATE_DEPTHS. A pair is dropped for a given depth if either cast in the
# pair does not physically reach z_bot (no extrapolation).
#
# Purpose: pick the integration depth that makes GOP (net diel change +
# flux-corrected night loss, both 18-day smoothed) easiest to interpret next
# to the float's own nitrate-drawdown NCP -- see Scripts/gop_vs_ncp_figure.py.
#
# Output: Data/Processed/o2_budget_3902681_multidepth.csv
#   columns: depth_m, type (net/night), mtime, dt_days,
#            rate_mmol_o2_m2_d, Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d
#            (rate_corr = rate - Fas, diffusive-only correction)
#
# Run from repo root with R:
#   & "C:/Users/flapet/AppData/Local/Programs/R/R-4.3.3/bin/x64/Rscript.exe" Scripts/o2_budget_3902681_multidepth.R
# -----------------------------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(lubridate)
library(gsw)

# ---- paths ------------------------------------------------------------------
nc_path   <- "Data/Raw/Floats/3902681_Sprof.nc"
era5_path <- "Data/Raw/ERA5/era5_wind_slp_3902681.nc"
out_csv   <- "Data/Processed/o2_budget_3902681_multidepth.csv"

CANDIDATE_DEPTHS <- c(40, 100, 150, 200)

# =============================================================================
# Step 1 - per-profile full O2 profile + surface state + MLD (unused here but
# kept for parity with the parent script's profile table)
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
z_max_need <- max(CANDIDATE_DEPTHS)

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
o2_prof_list <- vector("list", n_prof)   # depth-ordered full O2 profile, per prof_index
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
  max_depth = max_depth_vec
) |>
  mutate(
    lst   = (hour(time) + minute(time) / 60 + lon / 15) %% 24,
    phase = ifelse(lst < 12, "dawn", "dusk")
  ) |>
  filter(!is.na(o2_surf)) |>
  arrange(time)

# same edge-trimming as the parent script (irregular near-midday casts)
n_drop <- 5
stopifnot(nrow(prof) > 2 * n_drop)
prof <- prof |> slice((n_drop + 1):(n() - n_drop))

message(sprintf("Profiles retained: %d / %d (after trimming %d at each end)",
                nrow(prof), n_prof, n_drop))
message(sprintf("Profiles reaching >= %d m: %d / %d (%.0f%%)",
                z_max_need,
                sum(prof$max_depth >= z_max_need, na.rm = TRUE), nrow(prof),
                100 * mean(prof$max_depth >= z_max_need, na.rm = TRUE)))

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

# diffusive-only air-sea O2 flux (Wanninkhof 2014), identical formulation to
# airsea_flux() in the parent script (surface boundary term, independent of
# the column integration depth). No Liang 2013 bubble term -- matches the
# paper's stated GOP convention.
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

make_pair_depth <- function(d, i, j, type, z_bot) {
  dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))
  inv_i <- integrate_o2(d$prof_index[i], z_bot)
  inv_j <- integrate_o2(d$prof_index[j], z_bot)
  rate  <- (inv_j - inv_i) / dt_days
  Fas   <- airsea_flux_diffusive(d, i, j)
  tibble(
    depth_m = z_bot, type = type,
    mtime = d$time[i] + (d$time[j] - d$time[i]) / 2,
    dt_days = dt_days,
    rate_mmol_o2_m2_d = rate,
    Fas_mmol_o2_m2_d = Fas,
    rate_corr_mmol_o2_m2_d = rate - Fas
  )
}

n <- nrow(prof)

# --- prev: every consecutive pair (needed to build "net" and "night" subsets)
budget_all <- map_dfr(CANDIDATE_DEPTHS, function(z_bot) {
  # phase-aware pairing (mirrors the parent script exactly)
  night_tbl <- map_dfr(seq_len(n - 1), function(k) {
    if (prof$phase[k] == "dusk" && prof$phase[k + 1] == "dawn") {
      dt <- as.numeric(difftime(prof$time[k + 1], prof$time[k], units = "days"))
      if (dt < 0.7) return(make_pair_depth(prof, k, k + 1, "night", z_bot))
    }
    tibble()
  })

  net_tbl <- map_dfr(c("dusk", "dawn"), function(ph) {
    sub <- prof |> filter(phase == ph)
    if (nrow(sub) < 2) return(tibble())
    map_dfr(seq_len(nrow(sub) - 1), ~ make_pair_depth(sub, .x, .x + 1, "net", z_bot))
  }) |>
    arrange(mtime)

  out <- bind_rows(night_tbl, net_tbl)

  n_ok_net   <- sum(is.finite(net_tbl$rate_corr_mmol_o2_m2_d))
  n_ok_night <- sum(is.finite(night_tbl$rate_corr_mmol_o2_m2_d))
  message(sprintf(
    "depth=%3d m: net pairs %d/%d finite (%.0f%%), night pairs %d/%d finite (%.0f%%)",
    z_bot, n_ok_net, nrow(net_tbl), 100 * n_ok_net / nrow(net_tbl),
    n_ok_night, nrow(night_tbl), 100 * n_ok_night / nrow(night_tbl)))

  out
})

budget_out <- budget_all |>
  filter(is.finite(rate_corr_mmol_o2_m2_d)) |>
  select(depth_m, type, mtime, dt_days, rate_mmol_o2_m2_d,
         Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d) |>
  arrange(depth_m, type, mtime)

write_csv(budget_out, out_csv)
cat(sprintf("Wrote %s (%d rows, depths: %s)\n", out_csv, nrow(budget_out),
            paste(CANDIDATE_DEPTHS, collapse = ", ")))
