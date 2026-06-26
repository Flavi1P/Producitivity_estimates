# o2_budget_3902681.R
# -----------------------------------------------------------------------------
# Profile-by-profile dissolved-oxygen budget for BGC-Argo float 3902681.
#
# Computes the rate of change of the 0-40 m O2 inventory between profiles, in
# mmol O2 m-2 d-1 (signed; negative = O2 loss). Three budget variants:
#   prev  - every consecutive profile pair
#   night - dusk -> dawn overnight pairs (respiration / night loss)
#   net   - consecutive same-phase pairs (net diel change)
#
# No air-sea flux, no O2->C /2 conversion. Outputs:
#   Data/Processed/o2_budget_3902681.csv
#   Output/o2_budget_3902681_comparison.png
#
# Run from repo root with R 4.5.2:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681.R
# -----------------------------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(lubridate)
library(gsw)
library(zoo)
library(castr)   # mld() - same MLD convention as Scripts/format_for_ncp.R
library(patchwork)
library(geosphere)   # great-circle inter-profile displacement

# ---- paths ------------------------------------------------------------------
nc_path        <- "Data/Raw/Floats/3902681_Sprof.nc"
ref_net_path   <- "Data/Processed/O2_float_net_change.csv"
ref_night_path <- "Data/Processed/O2_float_night_loss.csv"
out_csv        <- "Data/Processed/o2_budget_3902681.csv"
out_png        <- "Output/o2_budget_3902681_comparison.png"
out_csv_mld    <- "Data/Processed/o2_budget_3902681_mld.csv"
out_png_mld    <- "Output/o2_budget_3902681_mld_comparison.png"
out_png_vol    <- "Output/o2_budget_3902681_mld_volumetric.png"
out_csv_resp   <- "Data/Processed/o2_budget_3902681_respiration.csv"
out_png_resp   <- "Output/o2_budget_3902681_respiration.png"

# =============================================================================
# Step 1 - per-profile O2 inventory (0-40 m)
# =============================================================================
nc <- nc_open(nc_path)

lon  <- ncvar_get(nc, "LONGITUDE")              # length N_PROF
lat  <- ncvar_get(nc, "LATITUDE")               # length N_PROF
juld <- ncvar_get(nc, "JULD")                   # days since 1950-01-01
pres <- ncvar_get(nc, "PRES")                   # [N_LEVELS x N_PROF]
temp <- ncvar_get(nc, "TEMP")
psal <- ncvar_get(nc, "PSAL")
doxy <- ncvar_get(nc, "DOXY_ADJUSTED")          # umol/kg
doxy_qc <- ncvar_get(nc, "DOXY_ADJUSTED_QC")    # char matrix, one char per level

nc_close(nc)

# ---- ERA5 hourly wind + SLP (air-sea flux forcing) --------------------------
# Cleaned, classic-NetCDF file written by Scripts/download_era5_wind.py.
# Dims confirmed as [longitude x latitude x time]; time is "hours since
# 1900-01-01". Variables u10, v10 (m s-1), msl (Pa).
era5_path <- "Data/Raw/ERA5/era5_wind_slp_3902681.nc"
e_nc  <- nc_open(era5_path)
e_lon <- ncvar_get(e_nc, "longitude")
e_lat <- ncvar_get(e_nc, "latitude")
e_t   <- as.POSIXct(ncvar_get(e_nc, "time") * 3600, origin = "1900-01-01", tz = "UTC")
e_u10 <- ncvar_get(e_nc, "u10")      # [lon x lat x time]
e_v10 <- ncvar_get(e_nc, "v10")
e_msl <- ncvar_get(e_nc, "msl")      # Pa
nc_close(e_nc)
stopifnot(identical(dim(e_u10), c(length(e_lon), length(e_lat), length(e_t))))
message(sprintf("ERA5: %d hourly steps, %s -> %s",
                length(e_t), format(min(e_t)), format(max(e_t))))

# Interval-mean forcing at the grid cell nearest the pair midpoint position.
# Returns <U2> = mean of (u10^2 + v10^2) over the interval (preserves gustiness
# that drives quadratic gas exchange) and the interval-mean SLP.
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
grid <- 0:40

# DOXY_ADJUSTED precision used for per-pair uncertainty propagation (Fix 3, the
# MLD budget below). The adjusted optode is good to ~1 umol/kg; at ~1027 kg m-3
# that is ~1 mmol m-3. We use it as the uncertainty on the LAYER-MEAN O2, i.e.
# sigma_<O2>. The z_int / MLD factor cancels in volumetric units (so the
# volumetric rate uncertainty is depth-independent) and scales the areal sigma,
# which is exactly why deep winter pairs are the least certain.
SIG_O2 <- 1.0   # mmol O2 m-3

# DOXY_ADJUSTED_QC comes back as a character vector of length N_PROF, each
# element being an N_LEVELS-long string (one flag char per level).
parse_qc_col <- function(qc_entry, n_levels) {
  chars <- strsplit(qc_entry, "")[[1]]
  flags <- suppressWarnings(as.integer(chars))
  length(flags) <- n_levels   # pad with NA if short
  flags
}

inventory_vec <- rep(NA_real_, n_prof)
mld_vec       <- rep(NA_real_, n_prof)   # mixed-layer depth, m (castr, sigma0 0.03)
# Surface state for the air-sea flux (means over the top 0-10 m, mmol m-3 / degC)
o2_surf_vec <- rep(NA_real_, n_prof)   # surface O2 inventory mean, mmol m-3
o2_eq_vec   <- rep(NA_real_, n_prof)   # equilibrium O2 at 1 atm, mmol m-3
sst_vec     <- rep(NA_real_, n_prof)   # surface in-situ temperature, degC

# Full depth-ordered O2 profile per cast, kept so the MLD-based budget can
# re-integrate each profile to an arbitrary depth later (the 0-40 m inventory
# above is precomputed; the MLD budget needs a variable bottom). Indexed by the
# ORIGINAL profile index (prof_index), so it survives the prof-table filtering.
o2_prof_list <- vector("list", n_prof)

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

  # in-situ density -> O2 in mmol m-3
  SA  <- gsw_SA_from_SP(ss, pr, lon[p], lat[p])
  CT  <- gsw_CT_from_t(SA, tt, pr)
  rho <- gsw_rho(SA, CT, pr)            # kg m-3
  o2_vol <- oo * rho / 1000            # mmol O2 m-3

  # depth from pressure (positive metres down)
  depth <- gsw_z_from_p(pr, lat[p]) * -1

  # reorder ALL level vectors by depth so surface T/S/rho are taken from the
  # genuinely shallowest levels (not a deep level) for the gas-exchange term.
  ord <- order(depth)
  depth <- depth[ord]; o2_vol <- o2_vol[ord]
  tt <- tt[ord]; SA <- SA[ord]; CT <- CT[ord]; rho <- rho[ord]

  # require coverage to 40 m and >= 5 valid levels in 0-60 m
  if (max(depth, na.rm = TRUE) < 40) next
  if (sum(depth <= 60) < 5) next

  o2_grid <- approx(x = depth, y = o2_vol, xout = grid, rule = 2)$y
  if (any(is.na(o2_grid))) next

  inventory_vec[p] <- sum(diff(grid) * (head(o2_grid, -1) + tail(o2_grid, -1)) / 2)

  # ---- mixed-layer depth (sigma0 0.03 criterion, castr) ---------------------
  # depth/SA/CT are already sorted shallow->deep above. Same call as
  # Scripts/format_for_ncp.R; default.depth = 300 when the criterion is never
  # reached (deep winter mixing beyond the cast or near-homogeneous column).
  sigma0 <- gsw_sigma0(SA, CT)
  mld_vec[p] <- castr::mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03,
                           default.depth = 300)

  # ---- store full depth-ordered O2 profile for the MLD-based budget ----------
  o2_prof_list[[p]] <- list(depth = depth, o2_vol = o2_vol)

  # ---- surface state for air-sea O2 flux (top 0-10 m means) -----------------
  top <- depth <= 10
  if (sum(top) == 0) top <- seq_len(min(3, length(depth)))  # fallback: shallowest
  SA_s  <- mean(SA[top]); CT_s <- mean(CT[top]); rho_s <- mean(rho[top])
  sst_vec[p]    <- mean(tt[top])                              # degC, in-situ
  o2_surf_vec[p] <- mean(o2_grid[grid <= 10])                 # mmol m-3
  # Garcia & Gordon (1992) solubility via gsw (umol/kg) -> mmol m-3 at 1 atm.
  o2_eq_umolkg  <- gsw_O2sol(SA_s, CT_s, 0, lon[p], lat[p])   # umol/kg
  o2_eq_vec[p]  <- o2_eq_umolkg * rho_s / 1000                # mmol m-3 (1 atm)
}

prof <- tibble(
  prof_index = seq_len(n_prof),
  time       = time_prof,
  lat        = lat,
  lon        = lon,
  inventory  = inventory_vec,
  mld        = mld_vec,       # mixed-layer depth, m
  o2_surf    = o2_surf_vec,   # surface O2, mmol m-3
  o2_eq      = o2_eq_vec,     # equilibrium O2 at 1 atm, mmol m-3
  sst        = sst_vec        # surface temperature, degC
) |>
  mutate(
    lst   = (hour(time) + minute(time) / 60 + lon / 15) %% 24,
    phase = ifelse(lst < 12, "dawn", "dusk")
  ) |>
  filter(!is.na(inventory)) |>
  arrange(time)

# Drop the first and last 5 profiles of the record: the early ones are the
# irregular near-midday casts (huge spurious diel rates) and we trim the tail
# symmetrically.
n_drop <- 5
stopifnot(nrow(prof) > 2 * n_drop)
prof <- prof |> slice((n_drop + 1):(n() - n_drop))

message(sprintf("Profiles with valid 0-40 m inventory: %d / %d (after trimming %d at each end)",
                nrow(prof), n_prof, n_drop))

# =============================================================================
# Step 2 - build the three budgets
# =============================================================================

# helper: diffusive air-sea O2 flux (Wanninkhof 2014) for an ordered pair.
# Surface boundary term -> independent of the integration depth, so it is shared
# by the 0-40 m budget and the MLD-based budget. Returns the pieces needed for
# both the CSV columns and rate_corr = rate - Fas.
airsea_flux <- function(d, i, j) {
  # Interval-mean surface temperature for the O2 Schmidt number (W2014 Table 1,
  # 35 psu seawater; salinity term not needed).
  T_s <- mean(c(d$sst[i], d$sst[j]))
  Sc  <- 1920.4 - 135.6 * T_s + 5.2122 * T_s^2 - 0.10939 * T_s^3 +
         0.00093777 * T_s^4
  lon_mid <- mean(c(d$lon[i], d$lon[j]))
  lat_mid <- mean(c(d$lat[i], d$lat[j]))
  w   <- era5_at(d$time[i], d$time[j], lon_mid, lat_mid)
  # Piston velocity, cm h-1 -> m d-1 (factor 0.24 = 24/100). <U2> = mean of
  # squared wind over the interval.
  k   <- 0.251 * w$u2 * (Sc / 660)^(-0.5) * 0.24            # m d-1

  # SLP correction of the equilibrium concentration to actual surface pressure
  # (first-order; water-vapour term omitted). msl in Pa, 101325 Pa = 1 atm.
  o2_eq_p_i <- d$o2_eq[i] * (w$msl / 101325)
  o2_eq_p_j <- d$o2_eq[j] * (w$msl / 101325)
  # dO2 > 0 = supersaturated. Both surface O2 and equilibrium are the average of
  # the pair's two profiles, representing the interval.
  dO2 <- mean(c(d$o2_surf[i], d$o2_surf[j])) - mean(c(o2_eq_p_i, o2_eq_p_j))
  Fas <- -k * dO2                                           # mmol O2 m-2 d-1, +into ocean
  list(k = k, u2 = w$u2, o2_eq = mean(c(o2_eq_p_i, o2_eq_p_j)),
       dO2 = dO2, Fas = Fas)
}

# helper: budget for an ordered pair of prof-table rows (later = j)
make_pair <- function(d, i, j, type) {
  dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))
  rate    <- (d$inventory[j] - d$inventory[i]) / dt_days   # signed dO2/dt

  fx <- airsea_flux(d, i, j)

  # NOTE: F_as is a surface boundary flux applied to the whole 0-40 m inventory.
  # This is exact only when MLD <= 40 m. In Iceland Basin / Irminger winter the
  # MLD is far deeper than 40 m, so gas exchange ventilates water below 40 m and
  # this simple correction over-attributes flux to the layer in winter. The
  # MLD-based budget below addresses this by integrating to the deeper of the two
  # MLDs (so F_as acts on the whole mixed layer); here we report the
  # unpartitioned 0-40 m correction.
  tibble(
    type        = type,
    mtime       = d$time[i] + (d$time[j] - d$time[i]) / 2,
    t_start     = d$time[i],
    t_end       = d$time[j],
    dt_days     = dt_days,
    phase_start = d$phase[i],
    phase_end   = d$phase[j],
    inv_start_mmol_m2 = d$inventory[i],
    inv_end_mmol_m2   = d$inventory[j],
    rate_mmol_o2_m2_d = rate,
    k_m_d             = fx$k,
    u2_m2_s2          = fx$u2,
    o2_eq_mmol_m3     = fx$o2_eq,
    delta_o2_mmol_m3  = fx$dO2,
    Fas_mmol_o2_m2_d  = fx$Fas,
    rate_corr_mmol_o2_m2_d = rate - fx$Fas
  )
}

# helper: integrate the stored O2 profile (prof_index pi) from 0 to z_bot (m),
# trapezoidal on a 1-m grid that ends exactly at z_bot. Returns NA if the cast
# does not reach z_bot (rule=2 would otherwise fabricate values below the cast).
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

# helper: layer-MEAN O2 (mmol m-3) of stored profile pi over [z_top, z_bot].
# Trapezoidal integral / thickness. Returns NA if the cast does not reach z_bot
# (same no-extrapolation guard as integrate_o2). Used by the moving mixed-layer
# budget for both the ML-mean O2 (0..MLD) and the just-below-MLD O2 (entrainment).
layer_mean_o2 <- function(pi, z_top, z_bot) {
  pr <- o2_prof_list[[pi]]
  if (is.null(pr)) return(NA_real_)
  if (z_bot <= z_top) return(NA_real_)
  if (max(pr$depth, na.rm = TRUE) < z_bot) return(NA_real_)
  g <- seq(z_top, z_bot, by = 1)
  if (g[length(g)] < z_bot) g <- c(g, z_bot)
  yv <- approx(x = pr$depth, y = pr$o2_vol, xout = g, rule = 2)$y
  if (any(is.na(yv))) return(NA_real_)
  sum(diff(g) * (head(yv, -1) + tail(yv, -1)) / 2) / (z_bot - z_top)
}

# helper: MLD-based budget for an ordered pair. Both profiles are integrated to a
# COMMON bottom = the deeper of the two MLDs (so the inventory difference is a
# like-for-like comparison), and never shallower than 40 m. Same air-sea flux as
# make_pair; here F_as acts on the whole mixed layer, which is the physically
# correct depth for a surface gas-exchange term.
make_pair_mld <- function(d, i, j, type) {
  z_int <- max(d$mld[i], d$mld[j], 40)
  inv_i <- integrate_o2(d$prof_index[i], z_int)
  inv_j <- integrate_o2(d$prof_index[j], z_int)
  dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))
  rate    <- (inv_j - inv_i) / dt_days                     # signed dO2/dt

  # great-circle horizontal displacement between the two casts (advection proxy)
  disp_km <- distHaversine(c(d$lon[i], d$lat[i]),
                           c(d$lon[j], d$lat[j])) / 1000
  dmld_m  <- d$mld[j] - d$mld[i]                            # MLD change (entrainment proxy)

  fx <- airsea_flux(d, i, j)

  # ---- Fix 1: volumetric (concentration) tendency on the common z_int --------
  # Within a pair both casts share z_int, so rate_areal = z_int * d<O2>/dt: the
  # z_int factor mechanically amplifies a few-mmol-m-3 signal into hundreds of
  # mmol m-2 d-1 in deep winter. Dividing out z_int reports the genuine
  # concentration tendency (mmol O2 m-3 d-1), comparable across seasons. The
  # surface flux diluted over the layer is Fas / z_int -- Cornec & Fassbender
  # (2025) Eq. 5, the standard mixed-layer-budget convention.
  rate_vol      <- rate / z_int
  Fas_vol       <- fx$Fas / z_int
  rate_corr_vol <- (rate - fx$Fas) / z_int

  # ---- Fix 2: moving mixed-layer budget with explicit entrainment ------------
  # Instead of differencing two casts over a COMMON z_int (a fixed control volume
  # whose depth changes between consecutive pairs, so the rates do not string into
  # one consistent inventory budget), track the mixed-layer-MEAN O2 over each
  # cast's OWN MLD and model the deepening explicitly. MLD floored at 40 m so
  # summer (MLD < 40) coincides with the 0-40 m budget. Mirrors the nitrate
  # entrainment in compute_ncp() (Scripts/ncp_function.R) and Cornec & Fassbender
  # (2025) Eq. 6, in concentration units.
  mld_i   <- max(d$mld[i], 40)
  mld_j   <- max(d$mld[j], 40)
  mld_bar <- mean(c(mld_i, mld_j))
  o2ml_i  <- layer_mean_o2(d$prof_index[i], 0, mld_i)      # ML-mean O2 start (mmol m-3)
  o2ml_j  <- layer_mean_o2(d$prof_index[j], 0, mld_j)      # ML-mean O2 end
  rate_ml_vol <- (o2ml_j - o2ml_i) / dt_days               # moving-ML tendency (mmol m-3 d-1)
  Fas_ml_vol  <- fx$Fas / mld_bar                          # surface flux diluted over the ML

  # Entrainment: deepening rate we (m d-1) draws water from just below the earlier
  # MLD into the layer. O2_below < O2ml (lower-O2 water below) => negative
  # tendency (ML O2 decreases) -- the sign check in the verification block.
  we       <- max(0, mld_j - mld_i) / dt_days              # deepening rate, m d-1 (0 if shoaling)
  o2_below <- layer_mean_o2(d$prof_index[i], mld_i, mld_i + 20)  # mean O2 in 20 m below earlier MLD
  entrain_vol <- if (we > 0 && is.finite(o2_below) && is.finite(o2ml_i)) {
    (we / mld_bar) * (o2_below - o2ml_i)                   # mmol m-3 d-1, Cornec Eq. 6
  } else 0
  # biological / uncorrected residual (concentration units)
  rate_resid_vol <- rate_ml_vol - Fas_ml_vol - entrain_vol

  # areal counterparts of the moving-ML budget (the headline units: x mixed-layer
  # depth). entrain_areal = we * (o2_below - o2ml_i) matches the compute_ncp()
  # diff = delta * we form.
  rate_ml_areal    <- rate_ml_vol    * mld_bar
  entrain_areal    <- entrain_vol    * mld_bar
  rate_resid_areal <- rate_resid_vol * mld_bar

  # ---- Fix 3: per-pair uncertainty (propagated layer-mean O2 precision) -------
  # sigma_rate_vol = sigma_<O2> * sqrt(2) / dt  (two independent ML means
  # differenced). Depth cancels in volumetric; the areal sigma scales with the ML
  # depth, so deep winter pairs are correctly the least certain and get
  # down-weighted in the error-weighted smoother (roll_smooth_w).
  sigma_rate_vol         <- SIG_O2 * sqrt(2) / dt_days     # mmol m-3 d-1
  sigma_rate_resid_areal <- sigma_rate_vol * mld_bar       # mmol m-2 d-1

  tibble(
    type        = type,
    mtime       = d$time[i] + (d$time[j] - d$time[i]) / 2,
    t_start     = d$time[i],
    t_end       = d$time[j],
    dt_days     = dt_days,
    phase_start = d$phase[i],
    phase_end   = d$phase[j],
    z_int_m     = z_int,
    mld_start_m = d$mld[i],
    mld_end_m   = d$mld[j],
    dmld_m      = dmld_m,
    disp_km     = disp_km,
    inv_start_mmol_m2 = inv_i,
    inv_end_mmol_m2   = inv_j,
    rate_mmol_o2_m2_d = rate,
    k_m_d             = fx$k,
    u2_m2_s2          = fx$u2,
    o2_eq_mmol_m3     = fx$o2_eq,
    delta_o2_mmol_m3  = fx$dO2,
    Fas_mmol_o2_m2_d  = fx$Fas,
    rate_corr_mmol_o2_m2_d = rate - fx$Fas,
    # Fix 1 - volumetric (concentration) versions of the common-z_int rate
    rate_vol_mmol_o2_m3_d      = rate_vol,
    Fas_vol_mmol_o2_m3_d       = Fas_vol,
    rate_corr_vol_mmol_o2_m3_d = rate_corr_vol,
    # Fix 2 - moving mixed-layer budget (own-MLD, floored at 40 m) + entrainment
    mld_bar_m              = mld_bar,
    o2ml_start_mmol_m3     = o2ml_i,
    o2ml_end_mmol_m3       = o2ml_j,
    o2_below_mmol_m3       = o2_below,
    we_m_d                 = we,
    rate_ml_vol_mmol_o2_m3_d       = rate_ml_vol,
    Fas_ml_vol_mmol_o2_m3_d        = Fas_ml_vol,
    entrain_vol_mmol_o2_m3_d       = entrain_vol,
    rate_resid_vol_mmol_o2_m3_d    = rate_resid_vol,
    rate_ml_mmol_o2_m2_d           = rate_ml_areal,
    entrain_mmol_o2_m2_d           = entrain_areal,
    rate_resid_mmol_o2_m2_d        = rate_resid_areal,
    # Fix 3 - propagated per-pair uncertainty
    sigma_rate_vol_mmol_o2_m3_d    = sigma_rate_vol,
    sigma_rate_resid_mmol_o2_m2_d  = sigma_rate_resid_areal
  )
}

n <- nrow(prof)

# --- prev: every consecutive pair --------------------------------------------
prev_tbl <- map_dfr(seq_len(n - 1), ~ make_pair(prof, .x, .x + 1, "prev"))

# --- night: consecutive dusk -> dawn within ~11 h ----------------------------
night_tbl <- prev_tbl |>
  filter(phase_start == "dusk", phase_end == "dawn", dt_days < 0.7) |>
  mutate(type = "night")

# --- net: consecutive same-phase pairs (dusk->dusk, dawn->dawn) --------------
net_tbl <- map_dfr(c("dusk", "dawn"), function(ph) {
  sub <- prof |> filter(phase == ph)
  if (nrow(sub) < 2) return(tibble())
  map_dfr(seq_len(nrow(sub) - 1), ~ make_pair(sub, .x, .x + 1, "net"))
}) |>
  arrange(mtime)

budget <- bind_rows(prev_tbl, night_tbl, net_tbl)

# =============================================================================
# Step 2c - MLD-based budgets (integrate to the deeper of the two MLDs, >= 40 m)
# =============================================================================
# Same three pairings as above, but each pair integrates BOTH profiles to a
# common bottom = max(mld_start, mld_end, 40). This deepens the inventory into
# the full mixed layer instead of the fixed 0-40 m slab, so winter ventilation
# below 40 m is captured and the air-sea flux acts on the correct depth.
prev_tbl_mld <- map_dfr(seq_len(n - 1), ~ make_pair_mld(prof, .x, .x + 1, "prev"))

night_tbl_mld <- prev_tbl_mld |>
  filter(phase_start == "dusk", phase_end == "dawn", dt_days < 0.7) |>
  mutate(type = "night")

net_tbl_mld <- map_dfr(c("dusk", "dawn"), function(ph) {
  sub <- prof |> filter(phase == ph)
  if (nrow(sub) < 2) return(tibble())
  map_dfr(seq_len(nrow(sub) - 1), ~ make_pair_mld(sub, .x, .x + 1, "net"))
}) |>
  arrange(mtime)

budget_mld <- bind_rows(prev_tbl_mld, night_tbl_mld, net_tbl_mld)

n_drop_mld <- sum(is.na(net_tbl_mld$rate_mmol_o2_m2_d))
message(sprintf(
  "MLD budget: integration depth median=%.0f m, max=%.0f m; %d/%d net pairs dropped (cast shallower than common MLD)",
  median(net_tbl_mld$z_int_m, na.rm = TRUE), max(net_tbl_mld$z_int_m, na.rm = TRUE),
  n_drop_mld, nrow(net_tbl_mld)))

# =============================================================================
# Step 3 - write CSV
# =============================================================================
budget_out <- budget |>
  select(type, mtime, t_start, t_end, dt_days, phase_start, phase_end,
         inv_start_mmol_m2, inv_end_mmol_m2, rate_mmol_o2_m2_d,
         k_m_d, u2_m2_s2, o2_eq_mmol_m3, delta_o2_mmol_m3,
         Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d)

write_csv(budget_out, out_csv)
cat(sprintf("Wrote %s (%d rows)\n", out_csv, nrow(budget_out)))

budget_mld_out <- budget_mld |>
  select(type, mtime, t_start, t_end, dt_days, phase_start, phase_end,
         z_int_m, mld_start_m, mld_end_m, dmld_m, disp_km,
         inv_start_mmol_m2, inv_end_mmol_m2, rate_mmol_o2_m2_d,
         k_m_d, u2_m2_s2, o2_eq_mmol_m3, delta_o2_mmol_m3,
         Fas_mmol_o2_m2_d, rate_corr_mmol_o2_m2_d,
         # Fix 1 - volumetric versions of the common-z_int rate
         rate_vol_mmol_o2_m3_d, Fas_vol_mmol_o2_m3_d, rate_corr_vol_mmol_o2_m3_d,
         # Fix 2 - moving mixed-layer budget + entrainment (vol + areal)
         mld_bar_m, o2ml_start_mmol_m3, o2ml_end_mmol_m3, o2_below_mmol_m3, we_m_d,
         rate_ml_vol_mmol_o2_m3_d, Fas_ml_vol_mmol_o2_m3_d,
         entrain_vol_mmol_o2_m3_d, rate_resid_vol_mmol_o2_m3_d,
         rate_ml_mmol_o2_m2_d, entrain_mmol_o2_m2_d, rate_resid_mmol_o2_m2_d,
         # Fix 3 - per-pair uncertainty
         sigma_rate_vol_mmol_o2_m3_d, sigma_rate_resid_mmol_o2_m2_d)

write_csv(budget_mld_out, out_csv_mld)
cat(sprintf("Wrote %s (%d rows)\n", out_csv_mld, nrow(budget_mld_out)))

# =============================================================================
# Reference series
# =============================================================================
matlab_to_posix <- function(mtime) {
  as.POSIXct((mtime - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
}

ref_net <- read_csv(ref_net_path, show_col_types = FALSE) |>
  mutate(time = matlab_to_posix(mtime)) |>
  transmute(time, value = Int_dO2dt_40m_mmol_m2_d, series = "net")

# Reference "night loss" is stored as a positive loss magnitude. We follow the
# comparison's sign convention: night is reported as a LOSS (positive = O2
# consumed overnight). Our computed night rate is dO2/dt (negative when O2
# drops), so we negate it to the loss convention wherever we compare/plot night.
ref_night <- read_csv(ref_night_path, show_col_types = FALSE) |>
  mutate(time = matlab_to_posix(mtime)) |>
  transmute(time, value = Int_O2_loss_40m_mmol_m2_d, series = "night")

# =============================================================================
# Step 5 - verification: correlation + bias vs reference (nearest within 1 day)
# =============================================================================
nearest_join <- function(mine, ref, tol_days = 1) {
  # mine, ref: tibbles with columns time (POSIXct) and value
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
  # slope through the origin: is mine a constant multiple of the reference?
  # slope ~1 => 1:1 agreement; slope ~2 => mine is twice the reference (i.e. a
  # 0.5 respiration coefficient applied on the reference side).
  slope <- coef(lm(value ~ 0 + ref_val, data = j))[1]
  cat(sprintf("%s: n=%d  Pearson r=%.3f  median(mine-ref)=%.3f  slope(mine~0+ref)=%.3f\n",
              label, nrow(j), r, bias, slope))
}

cat("\n--- Verification ---\n")
cat(sprintf("prev pairs:  %d\n", nrow(prev_tbl)))
cat(sprintf("night pairs: %d (ref rows: %d)\n", nrow(night_tbl), nrow(ref_night)))
cat(sprintf("net pairs:   %d (ref rows: %d)\n", nrow(net_tbl), nrow(ref_net)))
report_fit("NET   (mine vs ref net change)     ", net_tbl, ref_net)
# Night compared on the reference's LOSS convention: negate our dO2/dt to a loss.
report_fit("NIGHT (mine-loss vs ref loss)      ",
           night_tbl |> mutate(rate_mmol_o2_m2_d = -rate_mmol_o2_m2_d),
           ref_night)
# /2 check: does halving our estimate improve agreement with the reference?
report_fit("NIGHT/2 (mine-loss/2 vs ref loss)  ",
           night_tbl |> mutate(rate_mmol_o2_m2_d = -rate_mmol_o2_m2_d / 2),
           ref_night)
# MLD-based net change vs the 0-40 m reference. Expect agreement in summer (MLD
# <= 40 m, so the integration depth floors at 40 m and the two series coincide)
# and divergence in winter (MLD >> 40 m brings in sub-40 m ventilation).
cat(sprintf("MLD net: pairs at 40 m floor = %d / %d (%.0f%%)\n",
            sum(net_tbl_mld$z_int_m <= 40, na.rm = TRUE),
            sum(is.finite(net_tbl_mld$z_int_m)),
            100 * mean(net_tbl_mld$z_int_m <= 40, na.rm = TRUE)))
report_fit("MLD-NET (mine vs ref net change)   ", net_tbl_mld, ref_net)

# =============================================================================
# Step 4 - comparison plot
# =============================================================================
# Smoother helper: centered rolling mean (robust to irregular spacing here).
roll_smooth <- function(x, k = 8) {
  if (length(x) < k) return(x)
  zoo::rollapply(x, width = k, FUN = mean, na.rm = TRUE,
                 fill = NA, align = "center")
}

# Inverse-variance-weighted centered rolling mean (Fix 3). Weights w = 1/sigma^2,
# so noisy deep-winter pairs (large sigma_rate_resid, which scales with the ML
# depth) do not dominate the seasonal line. Falls back to NA windows with no
# weight. x and w must be the same length and pre-sorted by time.
roll_smooth_w <- function(x, w, k = 8) {
  if (length(x) < k) return(x)
  w[!is.finite(w) | !is.finite(x)] <- 0
  xx <- ifelse(is.finite(x), x, 0)
  num <- zoo::rollapply(xx * w, width = k, FUN = sum, na.rm = TRUE,
                        fill = NA, align = "center")
  den <- zoo::rollapply(w,      width = k, FUN = sum, na.rm = TRUE,
                        fill = NA, align = "center")
  out <- num / den
  out[!is.finite(out)] <- NA
  out
}

# Assemble long plotting table: source (mine/reference) x quantity (net/night)
mine_net <- net_tbl |>
  transmute(time = mtime, value = rate_mmol_o2_m2_d,
            quantity = "net change", source = "computed")
# Night on the reference's LOSS convention (positive = O2 consumed overnight):
# negate our dO2/dt.
mine_night <- night_tbl |>
  transmute(time = mtime, value = -rate_mmol_o2_m2_d,
            quantity = "night loss", source = "computed")
# /2 check: halve our estimate to see whether a 0.5 respiration coefficient
# would be needed to match the reference (it should NOT, slope is already ~1).
mine_night_half <- night_tbl |>
  transmute(time = mtime, value = -rate_mmol_o2_m2_d / 2,
            quantity = "night loss", source = "computed / 2")
# Flux-corrected series (rate_corr = rate - Fas). Same loss-convention flip for
# night as the raw computed series.
mine_net_corr <- net_tbl |>
  transmute(time = mtime, value = rate_corr_mmol_o2_m2_d,
            quantity = "net change", source = "computed (flux-corr)")
mine_night_corr <- night_tbl |>
  transmute(time = mtime, value = -rate_corr_mmol_o2_m2_d,
            quantity = "night loss", source = "computed (flux-corr)")
mine_prev <- prev_tbl |>
  transmute(time = mtime, value = rate_mmol_o2_m2_d,
            quantity = "prev (all pairs)", source = "computed")

ref_net_pl <- ref_net |>
  transmute(time, value, quantity = "net change", source = "reference")
# ref_night kept on its native loss convention (positive).
ref_night_pl <- ref_night |>
  transmute(time, value, quantity = "night loss", source = "reference")

plot_df <- bind_rows(mine_net, mine_night, mine_night_half,
                     mine_net_corr, mine_night_corr,
                     ref_net_pl, ref_night_pl, mine_prev) |>
  filter(is.finite(value)) |>
  arrange(quantity, source, time) |>
  group_by(quantity, source) |>
  mutate(value_smooth = roll_smooth(value, 8)) |>
  ungroup()

# Only the flux-corrected estimate vs the reference series.
src_levels <- c("computed (flux-corr)", "reference")
main_df  <- plot_df |>
  filter(quantity %in% c("net change", "night loss"), source %in% src_levels) |>
  mutate(source = factor(source, levels = src_levels))

# Smoothed seasonal lines only (no per-pair scatter) for a clean comparison.
p <- ggplot(main_df, aes(x = time, colour = quantity, linetype = source)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(aes(y = value_smooth), linewidth = 0.9) +
  scale_colour_manual(values = c("net change" = "#3690c0",
                                 "night loss" = "#cb181d")) +
  scale_linetype_manual(values = c("computed (flux-corr)" = "solid",
                                   "reference" = "dashed")) +
  # Clip residual diel noise so the seasonal signal stays legible.
  coord_cartesian(ylim = c(-200, 400)) +
  labs(
    title = "Float 3902681: O2 budget (0-40 m) - computed vs reference (+ air-sea flux correction)",
    subtitle = "Flux-corrected estimate (rate - F_as, Wanninkhof 2014) vs reference. Reference sign convention (night = loss, positive). Lines: 8-pt rolling mean.",
    x = "Date", y = expression("mmol O"[2] ~ "m"^-2 ~ "d"^-1),
    colour = "Quantity", linetype = "Source"
  ) +
  theme_bw()

ggsave(out_png, p, width = 11, height = 6, dpi = 200)
cat(sprintf("Wrote %s\n", out_png))

# =============================================================================
# Step 5b - moving mixed-layer budget vs 0-40 m comparison plot (net change)
# =============================================================================
# Top panel (AREAL headline, mmol O2 m-2 d-1): three computed net-change series
# against the reference -
#   * 0-40 m flux-corrected (the validated fixed-slab budget, unchanged);
#   * MLD common-z flux-corrected (the old deeper-of-two-MLDs differencing -
#     carries the winter depth amplification, shown as the "before");
#   * MLD moving + flux + entrainment corrected residual (rate_resid_mmol_o2_m2_d,
#     Fix 2) - the deepening of the mixed layer is now an explicit term instead of
#     sitting in the residual. This is the headline series.
# The MLD moving residual is error-weighted-smoothed (Fix 3): each pair's weight
# is 1/sigma_rate_resid^2, so the noisy deep-winter pairs do not dominate the
# line. The other series keep the plain 8-pt mean.
# Bottom panel: the common integration depth per MLD pair, so the divergence can
# be read against how deep the mixed layer was.
rate_40m <- net_tbl |>
  transmute(time = mtime, value = rate_corr_mmol_o2_m2_d,
            source = "0-40 m (flux-corr)", w = 1)
rate_mld_old <- net_tbl_mld |>
  transmute(time = mtime, value = rate_corr_mmol_o2_m2_d,
            source = "MLD common-z (flux-corr)", w = 1)
rate_mld_resid <- net_tbl_mld |>
  transmute(time = mtime, value = rate_resid_mmol_o2_m2_d,
            source = "MLD moving (flux+entrain-corr)",
            w = 1 / sigma_rate_resid_mmol_o2_m2_d^2)
rate_ref <- ref_net |>
  transmute(time, value, source = "reference", w = 1)

src_lv2 <- c("0-40 m (flux-corr)", "MLD common-z (flux-corr)",
             "MLD moving (flux+entrain-corr)", "reference")
rate_df <- bind_rows(rate_40m, rate_mld_old, rate_mld_resid, rate_ref) |>
  filter(is.finite(value)) |>
  mutate(source = factor(source, levels = src_lv2)) |>
  arrange(source, time) |>
  group_by(source) |>
  # error-weighted smoothing for the moving-ML residual; plain mean otherwise
  mutate(value_smooth = if (first(source) == "MLD moving (flux+entrain-corr)")
                          roll_smooth_w(value, w, 8) else roll_smooth(value, 8)) |>
  ungroup()

p_rate <- ggplot(rate_df, aes(time, colour = source, linetype = source)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(aes(y = value_smooth), linewidth = 0.9) +
  scale_colour_manual(values = c("0-40 m (flux-corr)" = "#3690c0",
                                 "MLD common-z (flux-corr)" = "#fb9a29",
                                 "MLD moving (flux+entrain-corr)" = "#238b45",
                                 "reference" = "grey30")) +
  scale_linetype_manual(values = c("0-40 m (flux-corr)" = "solid",
                                   "MLD common-z (flux-corr)" = "dotted",
                                   "MLD moving (flux+entrain-corr)" = "solid",
                                   "reference" = "dashed")) +
  coord_cartesian(ylim = c(-300, 600)) +
  labs(
    title = "Float 3902681: O2 net-change budget - 0-40 m vs moving mixed-layer (entrainment-resolved)",
    subtitle = paste("Areal (mmol O2 m-2 d-1). Green = moving-ML residual after air-sea AND entrainment removal,",
                     "error-weighted 8-pt mean (Fix 2/3);\norange dotted = old common-z flux-corr (before entrainment term)."),
    x = NULL, y = expression("mmol O"[2] ~ "m"^-2 ~ "d"^-1),
    colour = "Source", linetype = "Source"
  ) +
  theme_bw() +
  theme(legend.position = "top")

depth_df <- net_tbl_mld |>
  filter(is.finite(z_int_m)) |>
  transmute(time = mtime, z_int_m, mld_start_m, mld_end_m)

p_depth <- ggplot(depth_df, aes(time)) +
  geom_hline(yintercept = 40, colour = "grey60", linetype = "dashed") +
  geom_point(aes(y = pmax(mld_start_m, mld_end_m)), colour = "grey75",
             size = 0.5, alpha = 0.4) +
  geom_line(aes(y = roll_smooth(z_int_m, 8)), colour = "#238b45",
            linewidth = 0.8) +
  scale_y_reverse() +
  labs(
    subtitle = "Common integration depth per pair = max(MLD_start, MLD_end, 40 m). Grey = raw per-pair value; line = 8-pt mean; dashed = 40 m floor.",
    x = "Date", y = "Integration depth (m)"
  ) +
  theme_bw()

p_mld <- p_rate / p_depth + plot_layout(heights = c(2, 1))
ggsave(out_png_mld, p_mld, width = 11, height = 8, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_mld))

# =============================================================================
# Step 5d - moving mixed-layer budget decomposition (Fix 2 diagnostic)
# =============================================================================
# Show how the moving-ML areal tendency splits into the three modelled terms so
# the entrainment contribution is visible:
#   rate_ml  =  F_as  +  entrainment  +  residual (biological)
# All areal (mmol O2 m-2 d-1), 8-pt rolling means of the net pairs. Confirms the
# entrainment term is non-trivial in winter (deepening) and ~0 in summer.
out_png_decomp <- "Output/o2_budget_3902681_mld_decomposition.png"

decomp_df <- net_tbl_mld |>
  filter(is.finite(rate_ml_mmol_o2_m2_d)) |>
  arrange(mtime) |>
  transmute(
    time = mtime,
    `total ML tendency`     = rate_ml_mmol_o2_m2_d,
    `air-sea flux F_as`     = Fas_mmol_o2_m2_d,
    `entrainment`           = entrain_mmol_o2_m2_d,
    `residual (biological)` = rate_resid_mmol_o2_m2_d
  ) |>
  pivot_longer(-time, names_to = "term", values_to = "value") |>
  filter(is.finite(value)) |>
  mutate(term = factor(term, levels = c("total ML tendency", "air-sea flux F_as",
                                        "entrainment", "residual (biological)"))) |>
  arrange(term, time) |>
  group_by(term) |>
  mutate(value_smooth = roll_smooth(value, 8)) |>
  ungroup()

p_decomp <- ggplot(decomp_df, aes(time, colour = term)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(aes(y = value_smooth), linewidth = 0.9) +
  scale_colour_manual(values = c("total ML tendency"     = "grey30",
                                 "air-sea flux F_as"     = "#2166ac",
                                 "entrainment"           = "#d94801",
                                 "residual (biological)" = "#238b45")) +
  coord_cartesian(ylim = c(-400, 400)) +
  labs(
    title = "Float 3902681: moving mixed-layer O2 budget decomposition (areal)",
    subtitle = "total ML tendency = air-sea flux + entrainment + residual. Each: 8-pt rolling mean of net pairs. Entrainment < 0 when the ML deepens into lower-O2 water.",
    x = "Date", y = expression("mmol O"[2] ~ "m"^-2 ~ "d"^-1),
    colour = "Term"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(out_png_decomp, p_decomp, width = 11, height = 6, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_decomp))

# =============================================================================
# Step 5e - amplification diagnostic: areal vs volumetric spread (Fix 1)
# =============================================================================
# The depth amplification is a units artifact: rate_areal = z_int * d<O2>/dt, so
# deep-winter pairs carry a huge z_int factor. Dividing by z_int (volumetric)
# removes it. Compare the per-pair spread of the flux-corrected rate, areal vs
# volumetric, by season - the winter blow-up in areal collapses in volumetric.
out_png_amp <- "Output/o2_budget_3902681_amplification.png"

amp_df <- net_tbl_mld |>
  filter(is.finite(rate_corr_mmol_o2_m2_d), is.finite(rate_corr_vol_mmol_o2_m3_d)) |>
  mutate(season = ifelse(month(mtime) %in% c(11, 12, 1, 2, 3),
                         "winter (Nov-Mar)", "summer (Apr-Oct)")) |>
  transmute(season,
            `areal  (mmol O2 m-2 d-1)`     = rate_corr_mmol_o2_m2_d,
            `volumetric  (mmol O2 m-3 d-1)` = rate_corr_vol_mmol_o2_m3_d) |>
  pivot_longer(-season, names_to = "units", values_to = "value")

p_amp <- ggplot(amp_df, aes(season, value, fill = season)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
  facet_wrap(~units, scales = "free_y") +
  scale_fill_manual(values = c("winter (Nov-Mar)" = "#4292c6",
                               "summer (Apr-Oct)" = "#fdae6b"), guide = "none") +
  labs(
    title = "Float 3902681: depth amplification is a units artifact (areal vs volumetric)",
    subtitle = "Flux-corrected MLD net-change rate per pair. The winter spread that dominates the areal panel collapses in the volumetric (per-z_int) panel.",
    x = NULL, y = NULL
  ) +
  theme_bw()

ggsave(out_png_amp, p_amp, width = 10, height = 5, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_amp))

# =============================================================================
# Step 5f - volumetric (concentration) net-change timeseries (Fix 1 headline)
# =============================================================================
# The amplification-free view: the same net-change budget in mmol O2 m-3 d-1, so
# winter and summer are directly comparable. Three computed series plus the
# reference (the 0-40 m reference is areal -> divided by 40 to compare in
# concentration units):
#   * MLD moving residual (flux + entrainment corrected, rate_resid_vol) - the
#     headline volumetric quantity, error-weighted 8-pt mean (Fix 3);
#   * MLD common-z flux-corrected (rate_corr_vol) - same pairs, before Fix 2;
#   * 0-40 m flux-corrected, expressed volumetrically as rate_corr / 40;
#   * reference / 40.
# Per-pair points (grey) are shown under the residual line so the genuine
# concentration-level scatter is visible (not hidden) - cf. the huge areal spread.
rate_mld_resid_vol <- net_tbl_mld |>
  transmute(time = mtime, value = rate_resid_vol_mmol_o2_m3_d,
            source = "MLD moving (flux+entrain-corr)",
            w = 1 / sigma_rate_vol_mmol_o2_m3_d^2)
rate_mld_corr_vol <- net_tbl_mld |>
  transmute(time = mtime, value = rate_corr_vol_mmol_o2_m3_d,
            source = "MLD common-z (flux-corr)", w = 1)
rate_40m_vol <- net_tbl |>
  transmute(time = mtime, value = rate_corr_mmol_o2_m2_d / 40,
            source = "0-40 m (flux-corr) / 40", w = 1)
rate_ref_vol <- ref_net |>
  transmute(time, value = value / 40, source = "reference / 40", w = 1)

src_lv_vol <- c("0-40 m (flux-corr) / 40", "MLD common-z (flux-corr)",
                "MLD moving (flux+entrain-corr)", "reference / 40")
vol_df <- bind_rows(rate_mld_resid_vol, rate_mld_corr_vol,
                    rate_40m_vol, rate_ref_vol) |>
  filter(is.finite(value)) |>
  mutate(source = factor(source, levels = src_lv_vol)) |>
  arrange(source, time) |>
  group_by(source) |>
  mutate(value_smooth = if (first(source) == "MLD moving (flux+entrain-corr)")
                          roll_smooth_w(value, w, 8) else roll_smooth(value, 8)) |>
  ungroup()

# raw per-pair scatter of the headline volumetric residual (honesty, Fix 3)
vol_pts <- net_tbl_mld |>
  filter(is.finite(rate_resid_vol_mmol_o2_m3_d)) |>
  transmute(time = mtime, value = rate_resid_vol_mmol_o2_m3_d)

p_vol <- ggplot(vol_df, aes(time, colour = source, linetype = source)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(data = vol_pts, aes(time, value), inherit.aes = FALSE,
             colour = "grey75", size = 0.5, alpha = 0.4) +
  geom_line(aes(y = value_smooth), linewidth = 0.9) +
  scale_colour_manual(values = c("0-40 m (flux-corr) / 40" = "#3690c0",
                                 "MLD common-z (flux-corr)" = "#fb9a29",
                                 "MLD moving (flux+entrain-corr)" = "#238b45",
                                 "reference / 40" = "grey30")) +
  scale_linetype_manual(values = c("0-40 m (flux-corr) / 40" = "solid",
                                   "MLD common-z (flux-corr)" = "dotted",
                                   "MLD moving (flux+entrain-corr)" = "solid",
                                   "reference / 40" = "dashed")) +
  coord_cartesian(ylim = c(-6, 6)) +
  labs(
    title = "Float 3902681: O2 net-change budget in concentration units (volumetric)",
    subtitle = paste("mmol O2 m-3 d-1 (= areal / layer depth, Fix 1). Winter and summer are now directly comparable.",
                     "Grey = raw per-pair residual;\ngreen line = moving-ML residual, error-weighted 8-pt mean. Reference and 0-40 m divided by 40 m to compare."),
    x = "Date", y = expression("mmol O"[2] ~ "m"^-3 ~ "d"^-1),
    colour = "Source", linetype = "Source"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(out_png_vol, p_vol, width = 11, height = 6, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_vol))

# =============================================================================
# Step 5g - nighttime respiration index R (clean dark dusk->dawn budget)
# =============================================================================
# Respiration is isolated from the CLEAN OVERNIGHT window only. This float's
# sampling gives a clean ~11 h dark interval (dusk -> dawn, dt ~ 0.47 d) but NO
# clean daytime window (dawn -> dusk is ~37 h and spans day+night+day), so the
# dark period is the only place photosynthesis is guaranteed zero. In the dark,
# the moving mixed-layer O2 budget reduces to
#       d<O2>/dt = -R + F_as/MLD + entrainment
# so the biological residual (rate_resid = rate - F_as - entrainment, already
# computed by make_pair_mld) equals -R. The respiration index is therefore
#       R = -rate_resid       (positive = O2 consumed by the community)
# reported both areal (mmol O2 m-2 d-1) and volumetric (mmol O2 m-3 d-1, the
# amplification-free version, Fix 1). Per-pair uncertainty (Fix 3) carries
# straight over from sigma_rate_resid / sigma_rate_vol; the seasonal line is the
# inverse-variance-weighted 8-pt mean so deep-winter pairs do not dominate.
resp_tbl <- night_tbl_mld |>
  mutate(
    R_mmol_o2_m2_d           = -rate_resid_mmol_o2_m2_d,
    R_vol_mmol_o2_m3_d       = -rate_resid_vol_mmol_o2_m3_d,
    sigma_R_mmol_o2_m2_d     = sigma_rate_resid_mmol_o2_m2_d,
    sigma_R_vol_mmol_o2_m3_d = sigma_rate_vol_mmol_o2_m3_d
  )

resp_out <- resp_tbl |>
  select(mtime, t_start, t_end, dt_days,
         mld_bar_m, we_m_d, o2ml_start_mmol_m3, o2ml_end_mmol_m3, o2_below_mmol_m3,
         Fas_mmol_o2_m2_d, entrain_mmol_o2_m2_d,
         R_mmol_o2_m2_d, R_vol_mmol_o2_m3_d,
         sigma_R_mmol_o2_m2_d, sigma_R_vol_mmol_o2_m3_d)
write_csv(resp_out, out_csv_resp)
cat(sprintf("Wrote %s (%d rows)\n", out_csv_resp, nrow(resp_out)))

# ---- timeseries plot: areal R (top) + volumetric R (bottom) -----------------
resp_pl <- resp_tbl |>
  filter(is.finite(R_mmol_o2_m2_d), is.finite(R_vol_mmol_o2_m3_d)) |>
  arrange(mtime) |>
  mutate(
    w_areal = 1 / sigma_R_mmol_o2_m2_d^2,
    w_vol   = 1 / sigma_R_vol_mmol_o2_m3_d^2,
    R_areal_smooth = roll_smooth_w(R_mmol_o2_m2_d,     w_areal, 8),
    R_vol_smooth   = roll_smooth_w(R_vol_mmol_o2_m3_d, w_vol,   8)
  )

p_resp_areal <- ggplot(resp_pl, aes(mtime)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(aes(y = R_mmol_o2_m2_d), colour = "grey75", size = 0.5, alpha = 0.4) +
  geom_line(aes(y = R_areal_smooth), colour = "#cb181d", linewidth = 0.9) +
  coord_cartesian(ylim = c(-200, 400)) +
  labs(
    title = "Float 3902681: nighttime community respiration index R (dark dusk->dawn budget)",
    subtitle = paste("R = -(rate - F_as - entrainment) from the clean ~11 h dark window.",
                     "Positive = O2 consumed. Grey = raw per-pair; line = inverse-variance-weighted 8-pt mean."),
    x = NULL, y = expression("R  (mmol O"[2] ~ "m"^-2 ~ "d"^-1 * ")")
  ) +
  theme_bw()

p_resp_vol <- ggplot(resp_pl, aes(mtime)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(aes(y = R_vol_mmol_o2_m3_d), colour = "grey75", size = 0.5, alpha = 0.4) +
  geom_line(aes(y = R_vol_smooth), colour = "#cb181d", linewidth = 0.9) +
  coord_cartesian(ylim = c(-4, 6)) +
  labs(
    subtitle = "Volumetric (amplification-free, Fix 1): R per unit volume of the mixed layer.",
    x = "Date", y = expression("R  (mmol O"[2] ~ "m"^-3 ~ "d"^-1 * ")")
  ) +
  theme_bw()

p_resp <- p_resp_areal / p_resp_vol + plot_layout(heights = c(1, 1))
ggsave(out_png_resp, p_resp, width = 11, height = 8, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_resp))

# ---- printed respiration sanity ---------------------------------------------
# (sumr is also defined later for the air-sea block; define locally so this
# section can run standalone in source order.)
sumr <- function(x) sprintf("min=%.1f  median=%.1f  mean=%.1f  max=%.1f",
                            min(x), median(x), mean(x), max(x))
cat("\n--- Nighttime respiration index R (dark dusk->dawn, corrected) ---\n")
rf <- resp_tbl |> filter(is.finite(R_mmol_o2_m2_d))
cat(sprintf("n = %d night pairs; median entrainment correction = %.2f mmol m-2 d-1 (|.| median %.2f)\n",
            nrow(rf), median(rf$entrain_mmol_o2_m2_d, na.rm = TRUE),
            median(abs(rf$entrain_mmol_o2_m2_d), na.rm = TRUE)))
cat(sprintf("R areal : %s  | frac R>0 (net O2 consumed) = %.0f%%\n",
            sumr(rf$R_mmol_o2_m2_d), 100 * mean(rf$R_mmol_o2_m2_d > 0)))
cat(sprintf("R vol   : %s\n", sumr(rf$R_vol_mmol_o2_m3_d)))

# ---- printed amplification + entrainment sanity (Fix 1/2/3) -----------------
cat("\n--- MLD budget amplification & entrainment sanity (net pairs) ---\n")
amp_stat <- net_tbl_mld |>
  filter(is.finite(rate_corr_mmol_o2_m2_d), is.finite(rate_corr_vol_mmol_o2_m3_d)) |>
  mutate(season = ifelse(month(mtime) %in% c(11, 12, 1, 2, 3), "winter", "summer"))
sd_by <- amp_stat |>
  group_by(season) |>
  summarise(sd_areal = sd(rate_corr_mmol_o2_m2_d, na.rm = TRUE),
            sd_vol   = sd(rate_corr_vol_mmol_o2_m3_d, na.rm = TRUE), .groups = "drop")
print(sd_by)
wa <- sd_by$sd_areal[sd_by$season == "winter"]; sa <- sd_by$sd_areal[sd_by$season == "summer"]
wv <- sd_by$sd_vol[sd_by$season == "winter"];   sv <- sd_by$sd_vol[sd_by$season == "summer"]
cat(sprintf("winter/summer SD ratio: areal = %.2f  ->  volumetric = %.2f (amplification removed)\n",
            wa / sa, wv / sv))

# entrainment sign check: deepening pairs (we>0) with lower O2 below the MLD
# (o2_below < o2ml_start) must give a NEGATIVE entrainment tendency.
ent_chk <- net_tbl_mld |>
  filter(we_m_d > 0, is.finite(o2_below_mmol_m3), is.finite(o2ml_start_mmol_m3),
         is.finite(entrain_vol_mmol_o2_m3_d), o2_below_mmol_m3 < o2ml_start_mmol_m3)
cat(sprintf("entrainment sign check (deepening into lower-O2 => entrain<0): %s (%d/%d pairs)\n",
            if (nrow(ent_chk) == 0) "no qualifying pairs"
            else if (all(ent_chk$entrain_vol_mmol_o2_m3_d <= 0)) "OK" else "VIOLATED",
            sum(ent_chk$entrain_vol_mmol_o2_m3_d <= 0), nrow(ent_chk)))

# =============================================================================
# Step 5h - vertically-resolved (10 m) nighttime respiration field R(z, t)
# =============================================================================
# Step 5g collapses the clean dark window into a single mixed-layer-integrated
# respiration index per night pair. Here we keep the depth axis: in a FIXED 10 m
# bin during the dark window photosynthesis is zero, so the bin's O2 tendency is
# respiration plus physical exchange (air-sea gas exchange + vertical turbulent
# diffusion; entrainment and advection are negligible at the 11 h timescale - see
# Step 5g entrainment ~0 and Step 5c displacement-vs-rate ~0):
#       dO2/dt(z) = -R(z) + S_airsea(z) + S_diff(z)
#   ->  R(z)      = S_airsea(z) + S_diff(z) - dO2/dt(z)   [mmol O2 m-3 d-1]
# Positive R = O2 consumed. Volumetric only (no z_int / MLD multiplier - that is
# the depth-amplification artifact the volumetric form avoids, Fix 1) and no
# O2->C /2. Reuses layer_mean_o2 (bin-mean O2, NA below the cast), airsea_flux,
# and the same dusk->dawn dt<0.7 d night filter as night_tbl.

# ---- parameters (tunable; documented in the plot subtitle) -------------------
Z_EDGES   <- seq(0, 200, by = 10)        # 21 edges -> 20 bins of 10 m
N_BIN     <- length(Z_EDGES) - 1L
KZ_M2_S   <- 1e-5                         # background vertical diffusivity, m2 s-1
KZ_M2_D   <- KZ_M2_S * 86400             # -> m2 d-1 (~0.864)
out_csv_respz <- "Data/Processed/o2_budget_3902681_respiration_profile.csv"
out_png_respz <- "Output/o2_budget_3902681_respiration_profile.png"
out_png_respz_comp <- "Output/o2_budget_3902681_respiration_profile_components.png"

# ---- one night pair -> long tibble (20 bins) --------------------------------
bin_respiration <- function(d, i, j) {
  pi_i <- d$prof_index[i]; pi_j <- d$prof_index[j]
  dt_days <- as.numeric(difftime(d$time[j], d$time[i], units = "days"))

  z_top <- Z_EDGES[-length(Z_EDGES)]
  z_bot <- Z_EDGES[-1]
  z_mid <- (z_top + z_bot) / 2

  # bin-mean O2 (mmol m-3) for each cast; NA where the cast is too shallow.
  o2_i <- vapply(seq_len(N_BIN), function(k) layer_mean_o2(pi_i, z_top[k], z_bot[k]), numeric(1))
  o2_j <- vapply(seq_len(N_BIN), function(k) layer_mean_o2(pi_j, z_top[k], z_bot[k]), numeric(1))

  # concentration tendency per bin (NA where either bin-mean is NA - no extrapolation)
  dO2dt <- (o2_j - o2_i) / dt_days

  # vertical-diffusion flux convergence: centered 2nd difference per cast
  # (dz^2 = 100 m2), ends (no neighbour) -> 0, then averaged over the two casts.
  curv_one <- function(o2) {
    cv <- rep(NA_real_, length(o2))
    if (length(o2) > 2L) {
      k <- 2:(length(o2) - 1L)
      cv[k] <- (o2[k - 1L] - 2 * o2[k] + o2[k + 1L]) / 100
    }
    cv[1L] <- 0; cv[length(o2)] <- 0
    cv
  }
  curv   <- (curv_one(o2_i) + curv_one(o2_j)) / 2
  S_diff <- KZ_M2_D * curv
  S_diff[!is.finite(S_diff)] <- 0          # missing neighbour -> drop the diffusion term

  # air-sea exchange enters only the mixed layer, diluted over it (Fas/MLD).
  mld_bar  <- mean(c(max(d$mld[i], 40), max(d$mld[j], 40)))
  fx       <- airsea_flux(d, i, j)
  S_airsea <- ifelse(z_mid < mld_bar, fx$Fas / mld_bar, 0)

  R_vol       <- S_airsea + S_diff - dO2dt
  sigma_R_vol <- SIG_O2 * sqrt(2) / dt_days   # depth-independent (Fix 3 convention)

  tibble(
    mtime    = d$time[i] + (d$time[j] - d$time[i]) / 2,
    t_start  = d$time[i],
    t_end    = d$time[j],
    dt_days  = dt_days,
    z_top    = z_top,
    z_bot    = z_bot,
    z_mid    = z_mid,
    mld_bar  = mld_bar,
    o2_start = o2_i,
    o2_end   = o2_j,
    dO2dt    = dO2dt,
    S_airsea = S_airsea,
    S_diff   = S_diff,
    R_vol    = R_vol,
    sigma_R_vol = sigma_R_vol
  )
}

# ---- build the night set (same filter as night_tbl) and map over it ----------
night_idx <- which(prof$phase[-n] == "dusk" & prof$phase[-1] == "dawn" &
                   as.numeric(difftime(prof$time[-1], prof$time[-n], units = "days")) < 0.7)
respz_tbl <- map_dfr(night_idx, ~ bin_respiration(prof, .x, .x + 1))

write_csv(respz_tbl, out_csv_respz)
cat(sprintf("Wrote %s (%d rows = %d night pairs x %d bins)\n",
            out_csv_respz, nrow(respz_tbl), length(night_idx), N_BIN))

# ---- plot: R(z, t) heatmap (full + confidence-masked) ------------------------
mld_line <- respz_tbl |> distinct(mtime, mld_bar) |> arrange(mtime)

base_heat <- function(df, subt) {
  ggplot(df, aes(x = mtime, y = z_mid)) +
    geom_tile(aes(fill = R_vol), height = 10) +
    geom_line(data = mld_line, aes(x = mtime, y = mld_bar), inherit.aes = FALSE,
              colour = "grey15", linewidth = 0.5) +
    geom_hline(yintercept = 40, linetype = "dashed", colour = "grey30") +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-6, 6), oob = scales::squish,
                         name = expression("R (mmol O"[2] ~ "m"^-3 ~ "d"^-1 * ")")) +
    coord_cartesian(ylim = c(200, 0)) +
    labs(subtitle = subt, x = "Date", y = "Depth (m)") +
    theme_bw()
}

# confidence-masked copy: blank bins at the per-pair detection limit (|R|<sigma).
respz_masked <- respz_tbl |>
  mutate(R_vol = ifelse(is.finite(R_vol) & abs(R_vol) >= sigma_R_vol, R_vol, NA_real_))

p_full <- base_heat(respz_tbl,
  "All bins. Black line = mean MLD; dashed = 40 m floor. White = no respiration / O2 source.") +
  labs(title = paste0(
    "Float 3902681: dark-window respiration R(z, t) - volumetric, air-sea + vertical-diffusion corrected\n",
    sprintf("Kz = %.0e m2 s-1, dusk->dawn pairs only, 10 m bins 0-200 m", KZ_M2_S)))
p_masked <- base_heat(respz_masked,
  "Confidence-masked: bins with |R| < sigma_R (per-pair detection limit) blanked - trust summer/interior, not winter.")

p_respz <- p_full / p_masked + plot_layout(heights = c(1, 1))
ggsave(out_png_respz, p_respz, width = 11, height = 9, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_respz))

# ---- visualisation check: per-bin time-smoothed R(z, t) ----------------------
# Each column above is a single night pair, so the field is high-frequency. Smooth
# each depth bin in time with the inverse-variance-weighted rolling mean (w=1/s^2,
# Fix 3 convention) to expose the coherent seasonal/interior signal underneath the
# per-pair noise - this is the headline view of the respiration FIELD.
out_png_respz_smooth <- "Output/o2_budget_3902681_respiration_profile_smooth.png"
respz_smooth <- respz_tbl |>
  filter(is.finite(R_vol)) |>
  arrange(z_mid, mtime) |>
  group_by(z_mid) |>
  mutate(R_vol_smooth = roll_smooth_w(R_vol, 1 / sigma_R_vol^2, 8)) |>
  ungroup()

p_smooth <- ggplot(respz_smooth, aes(mtime, z_mid, fill = R_vol_smooth)) +
  geom_tile(height = 10) +
  geom_line(data = mld_line, aes(x = mtime, y = mld_bar), inherit.aes = FALSE,
            colour = "grey15", linewidth = 0.5) +
  geom_hline(yintercept = 40, linetype = "dashed", colour = "grey30") +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-6, 6), oob = scales::squish,
                       name = expression("R (mmol O"[2] ~ "m"^-3 ~ "d"^-1 * ")")) +
  coord_cartesian(ylim = c(200, 0)) +
  labs(title = "Float 3902681: dark-window respiration R(z, t) - per-bin 8-pt error-weighted time mean",
       subtitle = paste("Volumetric, air-sea + vertical-diffusion corrected. Black = mean MLD; dashed = 40 m floor.",
                        "Expect a summer subsurface respiration maximum beneath a near-zero/noisy winter column."),
       x = "Date", y = "Depth (m)") +
  theme_bw()
ggsave(out_png_respz_smooth, p_smooth, width = 11, height = 5.5, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_respz_smooth))

# ---- visualisation check: the three terms of R(z, t) side by side ------------
# R = S_airsea + S_diff - dO2dt. Faceting the terms confirms S_airsea acts only
# in-ML, S_diff is small in the stratified interior, and dO2dt carries the signal.
comp_df <- respz_tbl |>
  transmute(mtime, z_mid,
            `S_airsea  (Fas/MLD, in-ML only)` = S_airsea,
            `S_diff  (Kz d2O2/dz2)`           = S_diff,
            `-dO2/dt  (tendency)`             = -dO2dt,
            `R = S_airsea + S_diff - dO2/dt`  = R_vol) |>
  pivot_longer(-c(mtime, z_mid), names_to = "term", values_to = "value") |>
  mutate(term = factor(term, levels = c("S_airsea  (Fas/MLD, in-ML only)",
                                        "S_diff  (Kz d2O2/dz2)",
                                        "-dO2/dt  (tendency)",
                                        "R = S_airsea + S_diff - dO2/dt")))

p_comp <- ggplot(comp_df, aes(mtime, z_mid, fill = value)) +
  geom_tile(height = 10) +
  facet_wrap(~term, ncol = 2) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-6, 6), oob = scales::squish,
                       name = expression("mmol O"[2] ~ "m"^-3 ~ "d"^-1)) +
  coord_cartesian(ylim = c(200, 0)) +
  labs(title = "Float 3902681: respiration-budget terms by depth and time (volumetric)",
       subtitle = "R(z,t) decomposition. S_airsea is non-zero only inside the mixed layer; S_diff activates at the pycnocline.",
       x = "Date", y = "Depth (m)") +
  theme_bw()

ggsave(out_png_respz_comp, p_comp, width = 12, height = 8, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_respz_comp))

# ---- printed sanity (mirror Step 5g style) ----------------------------------
cat("\n--- Vertically-resolved nighttime respiration R(z, t) (Step 5h) ---\n")
cat(sprintf("n = %d night pairs x %d bins = %d rows (%d with finite R)\n",
            length(night_idx), N_BIN, nrow(respz_tbl), sum(is.finite(respz_tbl$R_vol))))

rz <- respz_tbl |>
  filter(is.finite(R_vol)) |>
  mutate(season  = ifelse(month(mtime) %in% c(11, 12, 1, 2, 3), "winter", "summer"),
         in_ml   = z_mid < mld_bar)
seas <- rz |>
  group_by(season) |>
  summarise(median_R_vol = median(R_vol),
            frac_R_pos   = mean(R_vol > 0),
            .groups = "drop")
print(seas)

cat(sprintf("median |S_airsea| in-ML = %.3f | below-ML = %.3f (should be ~0 below ML)\n",
            median(abs(rz$S_airsea[rz$in_ml])),
            median(abs(rz$S_airsea[!rz$in_ml]))))

# stratified interior = summer, below the mixed layer: S_diff should be << |R|.
strat <- rz |> filter(season == "summer", !in_ml)
cat(sprintf("stratified interior (summer, below ML): median |S_diff| = %.3f vs median |R_vol| = %.3f\n",
            median(abs(strat$S_diff)), median(abs(strat$R_vol))))

# Kz sensitivity: S_diff scales linearly with Kz, so x10 Kz -> x10 |S_diff|.
cat(sprintf("Kz sensitivity: median |S_diff| = %.3f at Kz=%.0e -> %.3f at Kz=1e-4 (x10)\n",
            median(abs(rz$S_diff)), KZ_M2_S, median(abs(rz$S_diff)) * 10))

# =============================================================================
# Step 5c - advection vs entrainment diagnostic: displacement vs delta-MLD
# =============================================================================
# Is the strong winter budget variability driven by the float drifting across a
# horizontal O2 gradient (-> displacement) or by the mixed layer deepening and
# entraining low-O2 water (-> |delta MLD|)? Plot displacement vs |delta MLD| per
# consecutive (prev) MLD pair, coloured by the budget-rate magnitude (the
# variability we are trying to explain). Whichever axis the high-|rate| points
# line up with is the more likely driver. Season (winter = Nov-Mar) is faceted.
out_png_diag2 <- "Output/o2_budget_3902681_disp_vs_dmld.png"

diag2 <- prev_tbl_mld |>
  filter(is.finite(disp_km), is.finite(dmld_m), is.finite(rate_corr_mmol_o2_m2_d)) |>
  mutate(
    abs_dmld   = abs(dmld_m),
    abs_rate   = abs(rate_corr_mmol_o2_m2_d),
    month      = month(mtime),
    season     = ifelse(month %in% c(11, 12, 1, 2, 3), "winter (Nov-Mar)",
                        "summer (Apr-Oct)")
  )

p_diag2 <- ggplot(diag2, aes(disp_km, abs_dmld, colour = abs_rate)) +
  geom_point(size = 1.6, alpha = 0.8) +
  scale_colour_viridis_c(option = "magma", direction = -1, trans = "sqrt",
                         name = expression("|rate| (mmol O"[2] ~ "m"^-2 ~ "d"^-1 * ")")) +
  facet_wrap(~season) +
  labs(
    title = "Float 3902681: horizontal displacement vs |MLD change| per consecutive pair",
    subtitle = "Colour = |flux-corrected O2 budget rate| (the variability to explain). Advection driver -> hot points spread along x; entrainment driver -> hot points spread along y.",
    x = "Inter-profile displacement (km)",
    y = expression("|" * Delta * "MLD| between profiles (m)")
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(out_png_diag2, p_diag2, width = 11, height = 5.5, dpi = 200)
cat(sprintf("Wrote %s\n", out_png_diag2))

# Quantitative read: rank-correlate |rate| against each candidate driver,
# overall and for winter only (where the variability is strongest).
corr_block <- function(label, d) {
  if (nrow(d) < 5) { cat(sprintf("%s: too few rows (%d)\n", label, nrow(d))); return(invisible()) }
  cs_disp <- cor(d$abs_rate, d$disp_km,  method = "spearman", use = "complete.obs")
  cs_mld  <- cor(d$abs_rate, d$abs_dmld, method = "spearman", use = "complete.obs")
  cat(sprintf("%s (n=%d): Spearman |rate| vs displacement = %+.2f | vs |dMLD| = %+.2f\n",
              label, nrow(d), cs_disp, cs_mld))
}
cat("\n--- Advection vs entrainment diagnostic (prev MLD pairs) ---\n")
corr_block("ALL    ", diag2)
corr_block("WINTER ", diag2 |> filter(season == "winter (Nov-Mar)"))
corr_block("SUMMER ", diag2 |> filter(season == "summer (Apr-Oct)"))

# =============================================================================
# Step 6 - air-sea flux diagnostics (the new term, all along the way)
# =============================================================================
# Built from the dense `prev` series (every consecutive pair). Each facet is one
# intermediate quantity of the Wanninkhof-2014 diffusive flux, so the chain
# wind -> k -> disequilibrium -> F_as is visible end to end.
out_diag <- "Output/o2_budget_3902681_airsea_diagnostics.png"

diag <- prev_tbl |>
  arrange(mtime) |>
  mutate(
    `wind speed  sqrt<U2>  (m s-1)`      = sqrt(u2_m2_s2),
    `piston velocity  k  (m d-1)`        = k_m_d,
    `O2 saturation  (%)`                 = 100 * (o2_eq_mmol_m3 + delta_o2_mmol_m3) /
                                            o2_eq_mmol_m3,
    `surface O2 disequilibrium  dO2  (mmol m-3)` = delta_o2_mmol_m3,
    `air-sea O2 flux  F_as  (mmol m-2 d-1)`      = Fas_mmol_o2_m2_d
  ) |>
  select(mtime,
         `wind speed  sqrt<U2>  (m s-1)`,
         `piston velocity  k  (m d-1)`,
         `O2 saturation  (%)`,
         `surface O2 disequilibrium  dO2  (mmol m-3)`,
         `air-sea O2 flux  F_as  (mmol m-2 d-1)`) |>
  pivot_longer(-mtime, names_to = "panel", values_to = "value") |>
  filter(is.finite(value)) |>
  mutate(panel = factor(panel, levels = c(
    "wind speed  sqrt<U2>  (m s-1)",
    "piston velocity  k  (m d-1)",
    "O2 saturation  (%)",
    "surface O2 disequilibrium  dO2  (mmol m-3)",
    "air-sea O2 flux  F_as  (mmol m-2 d-1)"))) |>
  arrange(panel, mtime) |>
  group_by(panel) |>
  mutate(value_smooth = roll_smooth(value, 8)) |>
  ungroup()

# Reference lines: 0 for disequilibrium / flux, 100 % for saturation.
href <- tibble(
  panel = factor(c("O2 saturation  (%)",
                   "surface O2 disequilibrium  dO2  (mmol m-3)",
                   "air-sea O2 flux  F_as  (mmol m-2 d-1)"),
                 levels = levels(diag$panel)),
  yint  = c(100, 0, 0)
)

p_diag <- ggplot(diag, aes(mtime, value)) +
  geom_hline(data = href, aes(yintercept = yint), colour = "grey50",
             linetype = "dashed") +
  geom_point(colour = "grey70", size = 0.5, alpha = 0.35) +
  geom_line(aes(y = value_smooth), colour = "#2166ac", linewidth = 0.8) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  labs(
    title = "Float 3902681: air-sea O2 flux diagnostics (Wanninkhof 2014, diffusive)",
    subtitle = "Per consecutive-profile pair. F_as > 0 = into ocean (winter undersaturation); F_as < 0 = outgassing (summer supersaturation). Blue = 8-pt rolling mean.",
    x = "Date", y = NULL
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 8))

ggsave(out_diag, p_diag, width = 10, height = 11, dpi = 200)
cat(sprintf("Wrote %s\n", out_diag))

# ---- printed physical-sanity summary ----------------------------------------
cat("\n--- Air-sea flux sanity (prev pairs) ---\n")
fin <- prev_tbl |> filter(is.finite(Fas_mmol_o2_m2_d))
sumr <- function(x) sprintf("min=%.1f  median=%.1f  mean=%.1f  max=%.1f",
                            min(x), median(x), mean(x), max(x))
cat(sprintf("k_m_d            : %s\n", sumr(fin$k_m_d)))
cat(sprintf("o2_eq_mmol_m3    : %s\n", sumr(fin$o2_eq_mmol_m3)))
cat(sprintf("delta_o2_mmol_m3 : %s\n", sumr(fin$delta_o2_mmol_m3)))
cat(sprintf("Fas_mmol_o2_m2_d : %s\n", sumr(fin$Fas_mmol_o2_m2_d)))
cat(sprintf("sign check (delta>0 => Fas<0): %s\n",
            ifelse(all(sign(fin$delta_o2_mmol_m3) == -sign(fin$Fas_mmol_o2_m2_d) |
                       fin$delta_o2_mmol_m3 == 0), "OK", "VIOLATED")))
