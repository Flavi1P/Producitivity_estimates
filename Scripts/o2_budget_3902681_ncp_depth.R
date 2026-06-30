# o2_budget_3902681_ncp_depth.R
# -----------------------------------------------------------------------------
# DEPTH-RESOLVED Net Community Production (NCP) from the dissolved-O2 tracer of
# BGC-Argo float 3902681, on a 10 m depth-bin grid (0-200 m).
#
# Depth-resolved companion to Scripts/o2_budget_3902681.R (the validated
# *integrated* 0-40 m / MLD budgets). Method:
#   1. Build per-profile O2 concentration on 10 m bins (0-200 m).
#   2. Grid to a regular WEEKLY x 10 m field O2(z,t) (mean of all profiles in the
#      week) and lightly smooth in time. Gridding first is essential: a fixed-
#      depth O2 difference between two individual 48 h profiles is dominated by
#      internal heave / mesoscale advection and optode noise (column-wide
#      stripes), not biology. Averaging ~14 profiles/week and differencing the
#      SMOOTH field isolates the seasonal biological tendency.
#   3. Local O2 tendency dO2/dt(z,t) = centred time derivative of the field.
#   4. Remove the FULL air-sea O2 flux (Wanninkhof 2014 diffusive + Liang 2013
#      bubble injection, ERA5-forced), weekly-averaged and spread over the mixed
#      layer (a surface boundary flux homogenised through the ML; Cornec &
#      Fassbender 2025 Eq. 5). Below the ML there is no direct gas-exchange term.
#
#   NCP(z) = dO2/dt(z) - F_as_total/MLD   (z within the mixed layer)
#   NCP(z) = dO2/dt(z)                     (z below the mixed layer)
#
# Sign: NCP > 0 = net autotrophy (biological O2 production). Volumetric units
# mmol O2 m-3 d-1 for the transect; the column integral is mmol O2 m-2 d-1.
# Kept in O2 (like the parent script); a carbon equivalent (PQ = 1.45) is
# printed and annotated.
#
# The column integral is over the PRODUCTIVE LAYER: 0 -> z_int = max(MLD, 40 m)
# (repo NCP convention pmax(MLD, zeu); this float carries no light/chla so zeu
# defaults to 40 m). Because z_int >= MLD, integrating the volumetric air-sea
# term over the layer removes exactly the full surface flux F_as_total, so the
# productive-layer integral is the air-sea-corrected NCP of the active layer.
#
# NOT resolved (documented, not corrected): interior vertical diffusion and
# advection between bins. The correction applied is the full air-sea flux only,
# as requested; sub-ML physical fluxes are assumed small relative to the
# biological signal over the productive season (standard float-NCP assumption,
# e.g. Plant et al. 2016, Bushinsky & Emerson 2015). In deep-winter mixing the
# layer is hundreds of m, so the winter integral is the noisiest part.
#
# Outputs:
#   Data/Processed/o2_budget_3902681_ncp_depth.csv      (weekly, per-bin field)
#   Data/Processed/o2_budget_3902681_ncp_prodlayer.csv  (weekly prod-layer integral)
#   Output/o2_budget_3902681_ncp_transect.png           (transect + integral)
#
# Run from repo root with R 4.5.2:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681_ncp_depth.R
# -----------------------------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(lubridate)
library(gsw)
library(zoo)
library(castr)       # mld() - same convention as the parent script
library(patchwork)

# ---- paths ------------------------------------------------------------------
nc_path   <- "Data/Raw/Floats/3902681_Sprof.nc"
era5_path <- "Data/Raw/ERA5/era5_wind_slp_3902681.nc"
out_csv_depth <- "Data/Processed/o2_budget_3902681_ncp_depth.csv"
out_csv_int   <- "Data/Processed/o2_budget_3902681_ncp_prodlayer.csv"
out_png       <- "Output/o2_budget_3902681_ncp_transect.png"

# ---- depth-bin grid ---------------------------------------------------------
# The field/integration grid reaches 500 m so the dynamic productive-layer
# integral z_int = max(MLD, 40) is fully captured in deep winter (weekly MLD
# peaks ~481 m; casts reach a median 904 m). The transect is displayed to 250 m
# where the biological structure sits.
Z_MAX     <- 500
Z_DISP    <- 250
BIN       <- 10
bin_edges <- seq(0, Z_MAX, by = BIN)               # 0,10,...,500
zc        <- head(bin_edges, -1) + BIN / 2         # bin centres 5,15,...,495
n_bins    <- length(zc)

# ---- constants (Liang 2013 bubble flux + optode precision) ------------------
SIG_O2  <- 1.0        # mmol O2 m-3, adjusted-optode precision
XG_O2   <- 0.20946    # O2 mole fraction in dry air
RHO_AIR <- 1.225      # reference air density, kg m-3
good_qc <- c(1L, 2L, 5L, 8L)
PQ      <- 1.45       # O2:C photosynthetic quotient (repo convention for O2->C)

vpress_atm <- function(S, Tc) {       # Weiss & Price 1980 sat. vapour press, atm
  TK <- Tc + 273.15
  exp(24.4543 - 67.4509 * (100 / TK) - 4.8489 * log(TK / 100) - 0.000544 * S)
}
parse_qc_col <- function(qc_entry, n_levels) {
  flags <- suppressWarnings(as.integer(strsplit(qc_entry, "")[[1]]))
  length(flags) <- n_levels
  flags
}

# =============================================================================
# Read float NetCDF
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

# =============================================================================
# Read ERA5 hourly wind + SLP; [longitude x latitude x time], "hours since 1900"
# =============================================================================
e_nc  <- nc_open(era5_path)
e_lon <- ncvar_get(e_nc, "longitude")
e_lat <- ncvar_get(e_nc, "latitude")
e_t   <- as.POSIXct(ncvar_get(e_nc, "time") * 3600, origin = "1900-01-01", tz = "UTC")
e_u10 <- ncvar_get(e_nc, "u10")
e_v10 <- ncvar_get(e_nc, "v10")
e_msl <- ncvar_get(e_nc, "msl")
nc_close(e_nc)
stopifnot(identical(dim(e_u10), c(length(e_lon), length(e_lat), length(e_t))))

era5_at <- function(t0, t1, lon0, lat0) {
  i  <- which.min(abs(e_lon - lon0))
  j  <- which.min(abs(e_lat - lat0))
  kt <- which(e_t >= t0 & e_t <= t1)
  if (length(kt) == 0) kt <- which.min(abs(as.numeric(e_t - (t0 + (t1 - t0) / 2))))
  u <- e_u10[i, j, kt]; v <- e_v10[i, j, kt]; m <- e_msl[i, j, kt]
  spd <- sqrt(u^2 + v^2)
  list(u2 = mean(u^2 + v^2, na.rm = TRUE), msl = mean(m, na.rm = TRUE),
       spd = spd[is.finite(spd)])
}

n_prof    <- length(juld)
time_prof <- as.POSIXct(juld * 86400, origin = "1950-01-01", tz = "UTC")

# =============================================================================
# Step 1 - per-profile O2 on the 10 m bin grid + MLD + surface state
# =============================================================================
bin_mean_o2 <- function(depth, o2_vol, z_top, z_bot, max_d) {
  if (z_bot > max_d) return(NA_real_)
  g  <- seq(z_top, z_bot, by = 1)
  yv <- approx(x = depth, y = o2_vol, xout = g, rule = 2)$y
  if (any(is.na(yv))) return(NA_real_)
  sum(diff(g) * (head(yv, -1) + tail(yv, -1)) / 2) / (z_bot - z_top)
}

o2_bins      <- matrix(NA_real_, nrow = n_prof, ncol = n_bins)
mld_vec      <- rep(NA_real_, n_prof)
o2_surf_vec  <- rep(NA_real_, n_prof)
o2_eq_vec    <- rep(NA_real_, n_prof)
sst_vec      <- rep(NA_real_, n_prof)
sss_vec      <- rep(NA_real_, n_prof)
rho_surf_vec <- rep(NA_real_, n_prof)

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
  o2_vol <- oo * rho / 1000                     # mmol O2 m-3
  depth  <- gsw_z_from_p(pr, lat[p]) * -1        # positive m down

  ord <- order(depth)
  depth <- depth[ord]; o2_vol <- o2_vol[ord]
  tt <- tt[ord]; ss <- ss[ord]; SA <- SA[ord]; CT <- CT[ord]; rho <- rho[ord]

  if (max(depth, na.rm = TRUE) < 40) next
  if (sum(depth <= 60) < 5) next
  max_d <- max(depth, na.rm = TRUE)

  for (b in seq_len(n_bins))
    o2_bins[p, b] <- bin_mean_o2(depth, o2_vol, bin_edges[b], bin_edges[b + 1], max_d)

  sigma0 <- gsw_sigma0(SA, CT)
  mld_vec[p] <- castr::mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03,
                           default.depth = 300)

  top <- depth <= 10
  if (sum(top) == 0) top <- seq_len(min(3, length(depth)))
  SA_s <- mean(SA[top]); CT_s <- mean(CT[top]); rho_s <- mean(rho[top])
  sst_vec[p]      <- mean(tt[top])
  sss_vec[p]      <- mean(ss[top])
  rho_surf_vec[p] <- rho_s
  o2_surf_vec[p]  <- bin_mean_o2(depth, o2_vol, 0, 10, max_d)
  o2_eq_vec[p]    <- gsw_O2sol(SA_s, CT_s, 0, lon[p], lat[p]) * rho_s / 1000
}

prof <- tibble(
  prof_index = seq_len(n_prof),
  time = time_prof, lat = lat, lon = lon,
  mld = mld_vec, o2_surf = o2_surf_vec, o2_eq = o2_eq_vec,
  sst = sst_vec, sss = sss_vec, rho_surf = rho_surf_vec
) |>
  mutate(has_surf = is.finite(o2_surf) & is.finite(mld)) |>
  filter(has_surf) |>
  arrange(time)

o2_bins_prof <- o2_bins[prof$prof_index, , drop = FALSE]

# trim first/last 5 profiles (irregular near-midday casts), as parent script
n_drop <- 5
stopifnot(nrow(prof) > 2 * n_drop)
keep_rows    <- (n_drop + 1):(nrow(prof) - n_drop)
prof         <- prof[keep_rows, ]
o2_bins_prof <- o2_bins_prof[keep_rows, , drop = FALSE]

message(sprintf("Profiles with valid surface+MLD: %d / %d (after trimming %d each end)",
                nrow(prof), n_prof, n_drop))
idx200 <- which(zc <= 200)
message(sprintf("Profiles fully covering 0-200 m / 0-500 m: %.0f%% / %.0f%%",
                100 * mean(apply(o2_bins_prof[, idx200, drop = FALSE], 1, function(r) all(is.finite(r)))),
                100 * mean(apply(o2_bins_prof, 1, function(r) all(is.finite(r))))))

# =============================================================================
# Step 2 - full air-sea O2 flux (W2014 diffusive + L13 bubble) per profile pair
# =============================================================================
airsea_total <- function(d, i, j) {
  T_s <- mean(c(d$sst[i], d$sst[j]))
  Sc  <- 1920.4 - 135.6 * T_s + 5.2122 * T_s^2 - 0.10939 * T_s^3 +
         0.00093777 * T_s^4
  w   <- era5_at(d$time[i], d$time[j],
                 mean(c(d$lon[i], d$lon[j])), mean(c(d$lat[i], d$lat[j])))
  k   <- 0.251 * w$u2 * (Sc / 660)^(-0.5) * 0.24            # m d-1 (W2014)

  o2_eq_p_i <- d$o2_eq[i] * (w$msl / 101325)
  o2_eq_p_j <- d$o2_eq[j] * (w$msl / 101325)
  dO2 <- mean(c(d$o2_surf[i], d$o2_surf[j])) - mean(c(o2_eq_p_i, o2_eq_p_j))
  Fas <- -k * dO2                                           # diffusive, +into ocean

  S_s    <- mean(c(d$sss[i], d$sss[j]))
  rho_w  <- mean(c(d$rho_surf[i], d$rho_surf[j]))
  o2eq1  <- mean(c(d$o2_eq[i], d$o2_eq[j]))
  o2surf <- mean(c(d$o2_surf[i], d$o2_surf[j]))
  ph2o   <- vpress_atm(S_s, T_s)
  slpc   <- ((w$msl / 101325) - ph2o) / (1 - ph2o)

  spd <- w$spd; if (length(spd) == 0) spd <- sqrt(w$u2)
  Cd <- ifelse(spd <= 11, 0.0012,
        ifelse(spd >= 20, 0.0018, 4.9e-4 + 6.5e-5 * spd))
  ustar_a <- spd * sqrt(Cd)
  ustar_w <- ustar_a / sqrt(rho_w / RHO_AIR)
  Kb_md   <- 5.5 * ustar_w^2.76 * (Sc / 660)^(-2 / 3) * 86400
  dP      <- 1.5244 * ustar_w^1.06
  Fc <- mean(XG_O2 * 5.56 * ustar_w^3.86 * 1000 * 86400)
  Fp <- mean(Kb_md * (o2eq1 * (1 + dP) * slpc - o2surf))

  fx_tot <- Fas + Fc + Fp
  tibble(mtime = d$time[i] + (d$time[j] - d$time[i]) / 2,
         Fas = Fas, Fas_bub = Fc + Fp, Fas_tot = fx_tot)
}

# per consecutive-pair flux (phase-independent surface term); weekly-averaged below
flux_pairs <- map_dfr(seq_len(nrow(prof) - 1), ~ airsea_total(prof, .x, .x + 1))

# =============================================================================
# Step 3 - weekly gridded O2(z,t) field, time derivative, air-sea correction
# =============================================================================
t0 <- min(prof$time)
week_of <- function(t) as.integer(floor(as.numeric(difftime(t, t0, units = "weeks"))))
prof <- prof |> mutate(week = week_of(time))

# long form of the per-profile bin matrix (column-major: bin-major blocks)
bins_long <- tibble(
  row = rep(seq_len(nrow(prof)), times = n_bins),
  z   = rep(zc, each = nrow(prof)),
  o2  = as.vector(o2_bins_prof)
) |>
  mutate(week = prof$week[row])

# weekly O2(z) field (mean of all profiles in the week) and weekly timestamp/MLD
field <- bins_long |>
  group_by(week, z) |>
  summarise(o2 = mean(o2, na.rm = TRUE), .groups = "drop")
wk_meta <- prof |>
  group_by(week) |>
  summarise(t = mean(time), mld = mean(mld, na.rm = TRUE), .groups = "drop")
flux_wk <- flux_pairs |>
  mutate(week = week_of(mtime)) |>
  group_by(week) |>
  summarise(Fas_tot = mean(Fas_tot, na.rm = TRUE), .groups = "drop")

# light 3-week centred smooth per depth, then centred time derivative dO2/dt
roll_t <- function(x, k = 3) {
  if (sum(is.finite(x)) < 2) return(x)
  zoo::rollapply(x, width = k, FUN = function(v) mean(v, na.rm = TRUE),
                 fill = NA, align = "center", partial = TRUE)
}
field <- field |>
  left_join(wk_meta, by = "week") |>
  left_join(flux_wk, by = "week") |>
  arrange(z, week) |>
  group_by(z) |>
  mutate(
    o2_s   = roll_t(o2, 3),
    dt_c   = as.numeric(difftime(lead(t), lag(t), units = "days")),
    dO2dt  = (lead(o2_s) - lag(o2_s)) / dt_c          # mmol m-3 d-1, centred
  ) |>
  ungroup() |>
  mutate(
    # full air-sea flux spread over the mixed layer (volumetric), ML bins only
    airsea_vol = ifelse(z <= mld, Fas_tot / mld, 0),
    ncp_vol    = dO2dt - airsea_vol                   # mmol O2 m-3 d-1, NCP>0=autotrophy
  )

ncp_depth <- field |>
  transmute(week, time = t, z_center_m = z, mld_m = mld,
            o2_mmol_m3 = o2_s,
            dO2dt_mmol_o2_m3_d = dO2dt,
            airsea_vol_mmol_o2_m3_d = airsea_vol,
            ncp_vol_mmol_o2_m3_d = ncp_vol) |>
  arrange(time, z_center_m)
write_csv(ncp_depth, out_csv_depth)
cat(sprintf("Wrote %s (%d rows)\n", out_csv_depth, nrow(ncp_depth)))

# =============================================================================
# Step 4 - productive-layer integral (areal) per week: 0 -> z_int = max(MLD,40)
# =============================================================================
# Integrate the volumetric NCP from the surface to the productive layer depth
# z_int = max(MLD, 40 m) (repo NCP convention pmax(MLD, zeu); no light data on
# this float, so zeu defaults to 40 m). The bottom bin is included fractionally
# (overlap thickness within [0, z_int]). Because z_int >= MLD, the full air-sea
# flux is removed from the column (integral of Fas_tot/MLD over the ML = Fas_tot).
# Require every overlapped bin to be finite (else the week is dropped, not
# extrapolated).
integ_prod <- function(ncp_vol, overlap) {
  ok <- overlap > 0
  if (!any(ok) || !all(is.finite(ncp_vol[ok]))) return(NA_real_)
  sum(ncp_vol[ok] * overlap[ok])
}

ncp_int <- field |>
  mutate(bin_top = z - BIN / 2,
         z_int   = pmax(mld, 40),
         overlap = pmax(pmin(z_int, bin_top + BIN) - bin_top, 0)) |>  # bin thickness within [0,z_int]
  group_by(week, t, mld, z_int, Fas_tot) |>
  summarise(ncp_prod = integ_prod(ncp_vol, overlap), .groups = "drop") |>
  transmute(week, time = t, mld_m = mld, z_int_m = z_int,
            Fas_total_mmol_o2_m2_d = Fas_tot,
            ncp_prod_mmol_o2_m2_d = ncp_prod) |>
  arrange(time) |>
  filter(is.finite(ncp_prod_mmol_o2_m2_d))
write_csv(ncp_int, out_csv_int)
cat(sprintf("Wrote %s (%d rows)\n", out_csv_int, nrow(ncp_int)))

# headline numbers (equal weight per week)
mean_ncp  <- mean(ncp_int$ncp_prod_mmol_o2_m2_d, na.rm = TRUE)
med_ncp   <- median(ncp_int$ncp_prod_mmol_o2_m2_d, na.rm = TRUE)
mean_prod <- mean(ncp_int$ncp_prod_mmol_o2_m2_d[month(ncp_int$time) %in% 4:9], na.rm = TRUE)

# =============================================================================
# Step 5 - the diagnostic figure: NCP transect (depth x time) + 0-200 m integral
# =============================================================================
# colour limit from the displayed (0-250 m) field; clip the field to Z_DISP
trans <- field |> filter(is.finite(ncp_vol), z <= Z_DISP)
lim   <- as.numeric(quantile(abs(trans$ncp_vol), 0.97, na.rm = TRUE))
trans <- trans |> mutate(ncp_clip = pmax(pmin(ncp_vol, lim), -lim))
# productive-layer depth line = max(MLD, 40 m) per week (the integration depth)
zline <- wk_meta |> transmute(t, z_int = pmax(mld, 40))

p_trans <- ggplot(trans, aes(t, z, fill = ncp_clip)) +
  geom_tile() +
  geom_line(data = zline, aes(t, z_int), inherit.aes = FALSE,
            colour = "black", linewidth = 0.6) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_datetime(expand = c(0, 0), date_labels = "%b %Y") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim, lim),
                       name = expression(NCP~"(mmol O"[2]~m^-3~d^-1*")")) +
  coord_cartesian(ylim = c(Z_DISP, 0)) +
  labs(
    title = "Float 3902681 - depth-resolved NCP from O2 (air-sea corrected: W2014 + Liang 2013 bubble)",
    subtitle = "Weekly O2(z,t) field, centred dO2/dt, minus full air-sea flux over the mixed layer. Red = autotrophy. Black = productive-layer depth max(MLD,40 m).",
    x = NULL, y = "Depth (m)"
  ) +
  theme_bw() + theme(legend.position = "right")

roll_smooth <- function(x, k = 5) {
  if (sum(is.finite(x)) < 2) return(x)
  zoo::rollapply(x, width = k, FUN = function(v) mean(v, na.rm = TRUE),
                 fill = NA, align = "center", partial = TRUE)
}
ncp_int <- ncp_int |> arrange(time) |> mutate(ncp_s = roll_smooth(ncp_prod_mmol_o2_m2_d, 5))

p_int <- ggplot(ncp_int, aes(time)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(aes(y = ncp_prod_mmol_o2_m2_d), colour = "grey70", size = 0.8) +
  geom_line(aes(y = ncp_s), colour = "#08519c", linewidth = 0.9) +
  geom_hline(yintercept = mean_ncp, colour = "#08519c", linetype = "dashed") +
  scale_x_datetime(date_labels = "%b %Y") +
  labs(
    subtitle = sprintf(
      "Productive-layer NCP, 0->max(MLD,40 m) (weekly). Grey = weekly value; blue = 5-pt mean; dashed = record mean %.0f mmol O2 m-2 d-1 (=%.0f mmol C, PQ %.2f).",
      mean_ncp, mean_ncp / PQ, PQ),
    x = "Date", y = expression("NCP"[prod.layer]~"(mmol O"[2]~m^-2~d^-1*")")
  ) +
  theme_bw()

fig <- p_trans / p_int + plot_layout(heights = c(1.4, 1))
ggsave(out_png, fig, width = 12, height = 9, dpi = 200)
cat(sprintf("Wrote %s\n", out_png))

# =============================================================================
# Summary
# =============================================================================
cat("\n--- Depth-resolved NCP summary (float 3902681, productive layer, O2 tracer) ---\n")
cat(sprintf("Productive layer: 0 -> z_int = max(MLD, 40 m); weekly z_int median=%.0f m, max=%.0f m\n",
            median(ncp_int$z_int_m), max(ncp_int$z_int_m)))
cat(sprintf("Weeks integrated: %d\n", nrow(ncp_int)))
cat(sprintf("Air-sea flux: full (W2014 diffusive + Liang 2013 bubble), spread over MLD\n"))
cat(sprintf("Prod-layer NCP record mean   = %7.1f mmol O2 m-2 d-1 (%6.1f mmol C, PQ %.2f)\n",
            mean_ncp, mean_ncp / PQ, PQ))
cat(sprintf("Prod-layer NCP record median = %7.1f mmol O2 m-2 d-1\n", med_ncp))
cat(sprintf("Prod-layer NCP Apr-Sep mean  = %7.1f mmol O2 m-2 d-1 (%6.1f mmol C)\n",
            mean_prod, mean_prod / PQ))
cat(sprintf("Transect colour limit (97th pct |NCP|, 0-%.0f m): %.2f mmol O2 m-3 d-1\n", Z_DISP, lim))
