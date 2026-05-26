# Spatial functional clustering of monthly NPP, adapted from
# Ballari et al. (2018, Int. J. Climatology).  See plan at
# C:/Users/petit/.claude/plans/read-the-intl-journal-eager-unicorn.md
#
# Pipeline:
#   1.  Load Output/cmems_cbpm_monthly_0p25deg.nc
#   2.  Per-cell 12-month climatology over 2015-2023
#   3.  Coverage gate over Feb-Oct, pad winter to zero
#   4.  Fourier basis fit (period 12, nbasis = 7)
#   5.  Optional spatial-trend removal
#   6.  L2-norm matrix between curves
#   7.  Empirical + parametric trace-variogram (geofd)
#   8.  Hierarchical clustering (complete + Ward), K = 2..10
#   9.  Silhouette to pick K
#  10.  Bootstrap ARI stability over years
#  11.  Plots + CSV in Output/
#
# Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/npp_bioregions.R

suppressPackageStartupMessages({
  library(ncdf4)
  library(tidyverse)
  library(fda)
  library(geofd)
  library(cluster)
  library(mclust)
  library(RColorBrewer)
})

set.seed(42)

NC_FILE   <- "Output/cmems_cbpm_monthly_0p25deg.nc"
OUT_DIR   <- "Output"
N_BASIS   <- 7                # Fourier basis count (3 harmonics + intercept)
COV_FRAC  <- 7 / 9            # >= 7 valid months out of Feb..Oct
K_RANGE   <- 2:8              # candidate cluster counts
N_BOOT    <- 100              # bootstrap iterations for stability
DETREND   <- TRUE             # remove spatial trend in coefficients
MONTH_MID <- seq(0.5, 11.5, by = 1)  # month midpoints on [0, 12]

# Optional: cell-level exclusion based on a prior cluster assignment.
# When non-empty, read the prior CSV at PRIOR_CSV, drop cells whose
# prior cluster is in EXCLUDE_CLUSTERS, suffix all outputs with TAG,
# and rerun the full pipeline (silhouette picks K freely on the reduced
# set).  Leave EXCLUDE_CLUSTERS empty for the standard run.
EXCLUDE_CLUSTERS <- c(3)
PRIOR_CSV        <- "Output/npp_bioregions.csv"
# When FORCE_K is a positive integer, override silhouette selection.
# Matches Ballari's habit of choosing a local silhouette max rather than
# the global one when the partition is more interpretable.
FORCE_K <- 4
TAG_EXCL <- if (length(EXCLUDE_CLUSTERS)) {
  paste0("_noR", paste(EXCLUDE_CLUSTERS, collapse = "_"))
} else {
  ""
}
TAG_K <- if (!is.null(FORCE_K)) sprintf("_K%d", FORCE_K) else ""
TAG <- paste0(TAG_EXCL, TAG_K)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load NetCDF ---------------------------------------------------------
message("Loading ", NC_FILE)
nc        <- nc_open(NC_FILE)
lons      <- ncvar_get(nc, "longitude")
lats      <- ncvar_get(nc, "latitude")
t_units   <- ncatt_get(nc, "time", "units")$value
t_raw     <- ncvar_get(nc, "time")
# parse "<unit> since YYYY-MM-DD ..."
unit_word <- sub(" since.*$", "", t_units)
origin_str <- sub("^[A-Za-z]+ since ", "", t_units)
origin <- as.Date(substr(origin_str, 1, 10))
times <- switch(unit_word,
                "days"    = origin + t_raw,
                "hours"   = origin + t_raw / 24,
                "seconds" = origin + t_raw / 86400,
                stop("Unsupported time unit: ", unit_word))
npp_int   <- ncvar_get(nc, "npp_int")   # (lon, lat, time)  in R order
nc_close(nc)

# Reorder so dims are (time, lat, lon) for ease
stopifnot(dim(npp_int) == c(length(lons), length(lats), length(times)))
npp_int <- aperm(npp_int, c(3, 2, 1))   # -> (time, lat, lon)
months  <- as.numeric(format(times, "%m"))
years   <- as.numeric(format(times, "%Y"))
message("Grid: ", length(times), " months x ", length(lats),
        " lat x ", length(lons), " lon")

# --- 2. Per-cell climatology ------------------------------------------------
clim <- array(NA_real_, dim = c(12, length(lats), length(lons)))
for (m in 1:12) {
  clim[m, , ] <- apply(npp_int[months == m, , , drop = FALSE],
                       c(2, 3), mean, na.rm = TRUE)
}
clim[!is.finite(clim)] <- NA_real_

# --- 3. Coverage gate + winter padding --------------------------------------
prod_months <- 2:10
n_valid_prod <- apply(!is.na(clim[prod_months, , ]), c(2, 3), sum)
keep <- n_valid_prod >= round(COV_FRAC * length(prod_months))
message("Cells passing coverage gate: ", sum(keep), " / ", length(keep))

# Drop cells assigned to EXCLUDE_CLUSTERS in the prior run (if requested).
if (length(EXCLUDE_CLUSTERS)) {
  stopifnot(file.exists(PRIOR_CSV))
  prior <- read_csv(PRIOR_CSV, show_col_types = FALSE)
  drop_xy <- prior |> filter(cluster %in% EXCLUDE_CLUSTERS) |>
    dplyr::select(lon, lat)
  lat_vec <- as.numeric(lats); lon_vec <- as.numeric(lons)
  for (r in seq_len(nrow(drop_xy))) {
    j <- which.min(abs(lat_vec - drop_xy$lat[r]))
    i <- which.min(abs(lon_vec - drop_xy$lon[r]))
    keep[j, i] <- FALSE
  }
  message("After excluding clusters ", paste(EXCLUDE_CLUSTERS, collapse = ","),
          ": ", sum(keep), " cells remain (",
          nrow(drop_xy), " dropped).")
}
stopifnot(sum(keep) >= 1500)

# winter (Jan, Nov, Dec) -> 0 for kept cells
winter_months <- c(1, 11, 12)
for (m in winter_months) {
  layer <- clim[m, , ]
  layer[keep & is.na(layer)] <- 0
  clim[m, , ] <- layer
}
# remaining NA in productive months for kept cells: linear interp across months
for (i in which(keep, arr.ind = TRUE) |> as.data.frame() |> nrow() |> seq_len()) {
  idx <- which(keep, arr.ind = TRUE)[i, ]
  cycle <- clim[, idx[1], idx[2]]
  if (any(is.na(cycle))) {
    cycle <- approx(MONTH_MID, cycle, xout = MONTH_MID, rule = 2)$y
    clim[, idx[1], idx[2]] <- cycle
  }
}

# Flatten kept cells into a (12 x ncells) matrix and align (lon, lat) tables
idx_mat <- which(keep, arr.ind = TRUE)              # cols: lat_i, lon_i
ncells  <- nrow(idx_mat)
cell_df <- tibble(
  cell_id = seq_len(ncells),
  lat_i   = as.integer(idx_mat[, 1]),
  lon_i   = as.integer(idx_mat[, 2]),
  lat     = as.numeric(lats)[as.integer(idx_mat[, 1])],
  lon     = as.numeric(lons)[as.integer(idx_mat[, 2])],
)
clim_mat <- matrix(NA_real_, nrow = 12, ncol = ncells)
for (k in seq_len(ncells)) {
  clim_mat[, k] <- clim[, cell_df$lat_i[k], cell_df$lon_i[k]]
}
stopifnot(all(is.finite(clim_mat)))
message("Climatology matrix: ", nrow(clim_mat), " months x ",
        ncol(clim_mat), " cells")

# --- 4. Fourier basis fit ---------------------------------------------------
fb <- create.fourier.basis(rangeval = c(0, 12), nbasis = N_BASIS,
                           period = 12)
fd_obj <- Data2fd(argvals = MONTH_MID, y = clim_mat, basisobj = fb)
M_inprod <- inprod(fb, fb)
message("Fourier basis: nbasis = ", N_BASIS, ", period = 12")

# --- 5. Optional spatial trend removal --------------------------------------
# Coefficient-wise OLS: each basis coefficient (per cell) is regressed
# on (1, lat, lon).  Equivalent to functional regression on scalar
# covariates, and avoids fda's fRegress version differences.
detrend_coefs <- function(coefs, lat_v, lon_v) {
  X <- cbind(1, lat_v, lon_v)               # (ncells x 3)
  XtX_inv <- solve(crossprod(X))
  fits <- X %*% XtX_inv %*% crossprod(X, t(coefs))   # (ncells x nbasis)
  coefs - t(fits)
}

if (DETREND) {
  resid_coef <- detrend_coefs(fd_obj$coefs, cell_df$lat, cell_df$lon)
  fd_for_cluster <- fd(resid_coef, fb)
  message("Detrending applied (coefficient-wise OLS on lat + lon).")
} else {
  fd_for_cluster <- fd_obj
}

# --- 6. L2 norm matrix ------------------------------------------------------
# Vectorised replacement for geofd::l2.norm (which is a pure-R double
# loop, ~minutes per call on 2545 cells).  Identity used:
#   ||c_i - c_j||^2_M = c_i' M c_i + c_j' M c_j - 2 c_i' M c_j.
fast_l2sq <- function(coefs, M) {
  MC    <- M %*% coefs                # nbasis x ncells
  cross <- crossprod(coefs, MC)       # ncells x ncells: c_i'M c_j
  a     <- diag(cross)
  outer(a, a, "+") - 2 * cross
}

message("Computing L2 norm matrix (", ncells, " x ", ncells, ") ...")
L2sq   <- fast_l2sq(fd_for_cluster$coefs, M_inprod)
L2sq[L2sq < 0] <- 0                   # floor numerical noise
diag(L2sq) <- 0
L2dist <- sqrt(L2sq)

# --- 7. Trace-variogram -----------------------------------------------------
# coords are (lon, lat) in degrees.  geofd computes Euclidean distance
# between coords; for a small region at 60 N this is acceptable
# (1 deg lon ~ 55 km, 1 deg lat ~ 111 km) -- variogram interpretation
# is in degrees-of-arc, not km.
coords <- as.matrix(cell_df[, c("lon", "lat")])
max_d  <- 0.6 * max(dist(coords))
message("Empirical trace-variogram, max.dist = ",
        sprintf("%.2f", max_d), " deg")
emp_tv <- trace.variog(coords = coords, L2norm = L2sq,
                       bin = TRUE, max.dist = max_d)
# geofd::trace.variog returns u (bin centres) and v (semivariance).
# Drop empty bins (NA) before passing to fit.tracevariog.
keep_bin <- is.finite(emp_tv$u) & is.finite(emp_tv$v)
emp_tv_clean <- emp_tv
emp_tv_clean$u <- emp_tv$u[keep_bin]
emp_tv_clean$v <- emp_tv$v[keep_bin]
emp_dist  <- emp_tv_clean$u
emp_gamma <- emp_tv_clean$v
sill0 <- mean(tail(emp_gamma, max(3, length(emp_gamma) %/% 4)),
              na.rm = TRUE)
fit_tv <- tryCatch(
  fit.tracevariog(emp.trace.vari = emp_tv_clean,
                  models = c("spherical", "exponential", "gaussian"),
                  sigma2.0 = sill0, phi.0 = max_d / 3),
  error = function(e) {
    message("fit.tracevariog failed: ", conditionMessage(e))
    NULL
  }
)
best_model <- if (!is.null(fit_tv)) fit_tv$best$cov.model else NA
best_range <- if (!is.null(fit_tv)) fit_tv$best$cov.pars[2] else NA
message("Best variogram model: ", best_model,
        " | range = ", sprintf("%.2f", best_range), " deg")

# --- 8. Hierarchical clustering ---------------------------------------------
dmat <- as.dist(L2dist)
hc_complete <- hclust(dmat, method = "complete")
hc_ward     <- hclust(dmat, method = "ward.D2")

# --- 9. Silhouette over K ---------------------------------------------------
sil_table <- tibble()
for (K in K_RANGE) {
  for (method in c("complete", "ward.D2")) {
    hc <- if (method == "complete") hc_complete else hc_ward
    cl <- cutree(hc, k = K)
    s  <- mean(silhouette(cl, dmat)[, "sil_width"])
    sil_table <- bind_rows(sil_table, tibble(K = K, method = method,
                                             mean_sil = s))
  }
}
if (!is.null(FORCE_K)) {
  # Pick the linkage that maximises silhouette at the forced K.
  best <- sil_table |> filter(K == FORCE_K) |>
    arrange(desc(mean_sil)) |> slice(1)
  message("FORCE_K = ", FORCE_K, ": using K = ", best$K,
          ", method = ", best$method,
          ", mean width = ", sprintf("%.3f", best$mean_sil))
} else {
  best <- sil_table |> arrange(desc(mean_sil)) |> slice(1)
  message("Best silhouette: K = ", best$K, ", method = ", best$method,
          ", mean width = ", sprintf("%.3f", best$mean_sil))
}

hc_best   <- if (best$method == "complete") hc_complete else hc_ward
cluster_v <- cutree(hc_best, k = best$K)
cell_df$cluster <- cluster_v

if (best$mean_sil < 0.30)
  warning("Mean silhouette < 0.30 -- bioregion partition is weak.")

# --- 10. Bootstrap stability ------------------------------------------------
message("Bootstrap stability (", N_BOOT, " iterations) ...")
ari_vec <- numeric(N_BOOT)
unique_years <- sort(unique(years))
for (b in seq_len(N_BOOT)) {
  yrs <- sample(unique_years, replace = TRUE)
  # rebuild climatology with resampled years
  npp_sub <- npp_int[years %in% yrs, , , drop = FALSE]
  m_sub   <- months[years %in% yrs]
  clim_b  <- array(NA_real_, dim = c(12, length(lats), length(lons)))
  for (m in 1:12) {
    clim_b[m, , ] <- apply(npp_sub[m_sub == m, , , drop = FALSE],
                           c(2, 3), mean, na.rm = TRUE)
  }
  # winter pad + interp for kept cells only
  clim_b[!is.finite(clim_b)] <- NA_real_
  for (m in winter_months) {
    layer <- clim_b[m, , ]
    layer[keep & is.na(layer)] <- 0
    clim_b[m, , ] <- layer
  }
  cmat_b <- matrix(NA_real_, nrow = 12, ncol = ncells)
  for (k in seq_len(ncells)) {
    cyc <- clim_b[, cell_df$lat_i[k], cell_df$lon_i[k]]
    if (any(is.na(cyc)))
      cyc <- approx(MONTH_MID, cyc, xout = MONTH_MID, rule = 2)$y
    cmat_b[, k] <- cyc
  }
  if (any(!is.finite(cmat_b))) { ari_vec[b] <- NA; next }
  fd_b <- Data2fd(argvals = MONTH_MID, y = cmat_b, basisobj = fb)
  if (DETREND) {
    fd_b <- fd(detrend_coefs(fd_b$coefs, cell_df$lat, cell_df$lon), fb)
  }
  L2_b <- fast_l2sq(fd_b$coefs, M_inprod)
  L2_b[L2_b < 0] <- 0
  d_b  <- as.dist(sqrt(L2_b))
  hc_b <- hclust(d_b, method = best$method)
  cl_b <- cutree(hc_b, k = best$K)
  ari_vec[b] <- adjustedRandIndex(cluster_v, cl_b)
  if (b %% 10 == 0) message("  boot ", b, "/", N_BOOT,
                            "  running ARI median = ",
                            sprintf("%.3f", median(ari_vec[1:b], na.rm = TRUE)))
}
ari_med <- median(ari_vec, na.rm = TRUE)
ari_q05 <- quantile(ari_vec, 0.05, na.rm = TRUE)
ari_q95 <- quantile(ari_vec, 0.95, na.rm = TRUE)
message(sprintf("Bootstrap ARI: median = %.3f, 5-95%% = [%.3f, %.3f]",
                ari_med, ari_q05, ari_q95))

writeLines(c(
  sprintf("K               : %d", best$K),
  sprintf("Linkage method  : %s", best$method),
  sprintf("Mean silhouette : %.3f", best$mean_sil),
  sprintf("Variogram model : %s (range = %.2f deg)", best_model, best_range),
  sprintf("Bootstrap ARI   : median = %.3f, 5%%-95%% = [%.3f, %.3f] (n=%d)",
          ari_med, ari_q05, ari_q95, N_BOOT)
), file.path(OUT_DIR, sprintf("npp_bioregions_stability%s.txt", TAG)))

# --- 11. Plots --------------------------------------------------------------
clu_palette <- brewer.pal(max(3, best$K), "Set1")[seq_len(best$K)]

# 11a. Cluster map
p_map <- ggplot(cell_df) +
  geom_tile(aes(x = lon, y = lat, fill = factor(cluster)),
            width = 0.25, height = 0.25) +
  borders("world", xlim = c(-45, -5), ylim = c(56, 67),
          fill = "grey90", colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = clu_palette, name = "Bioregion") +
  coord_quickmap(xlim = c(-41, -9), ylim = c(57, 66)) +
  labs(x = "Longitude", y = "Latitude",
       title = sprintf("NPP bioregions (K = %d, %s linkage)",
                       best$K, best$method),
       subtitle = sprintf("CbPM monthly 0.25 deg, 2015-2023 climatology"),
       caption = sprintf("Mean silhouette = %.3f | Bootstrap ARI = %.2f",
                         best$mean_sil, ari_med)) +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_clusters%s.png", TAG)),
       p_map, width = 9, height = 5, dpi = 140)

# 11b. Central curves per cluster (median + IQR)
# evaluate fd_obj (un-detrended) on a fine grid for plotting
month_grid <- seq(0, 12, length.out = 200)
curves_eval <- eval.fd(month_grid, fd_obj)  # 200 x ncells
curves_df <- as_tibble(curves_eval) |>
  setNames(as.character(cell_df$cell_id)) |>
  mutate(month = month_grid) |>
  pivot_longer(-month, names_to = "cell_id", values_to = "npp") |>
  mutate(cell_id = as.integer(cell_id)) |>
  left_join(cell_df |> dplyr::select(cell_id, cluster), by = "cell_id")

curve_summary <- curves_df |>
  group_by(cluster, month) |>
  summarise(med = median(npp),
            q25 = quantile(npp, 0.25),
            q75 = quantile(npp, 0.75),
            .groups = "drop")

month_lbl <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
p_curves <- ggplot(curve_summary) +
  geom_ribbon(aes(x = month, ymin = q25, ymax = q75,
                  fill = factor(cluster)), alpha = 0.25) +
  geom_line(aes(x = month, y = med, colour = factor(cluster)),
            linewidth = 1.0) +
  scale_colour_manual(values = clu_palette, guide = "none") +
  scale_fill_manual(values = clu_palette, name = "Bioregion") +
  scale_x_continuous(breaks = 0:11 + 0.5, labels = month_lbl,
                     limits = c(0, 12)) +
  labs(x = NULL, y = expression("NPP (mg C m"^-2*" d"^-1*")"),
       title = "Central seasonal cycle per bioregion",
       subtitle = "Median (line) and IQR (shaded) across cells in each cluster") +
  facet_wrap(~ cluster, ncol = best$K) +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_central_curves%s.png", TAG)),
       p_curves, width = 11, height = 4, dpi = 140)

# 11c. Silhouette vs K
p_sil <- ggplot(sil_table, aes(K, mean_sil, colour = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = best$K, linetype = "dashed", colour = "grey50") +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "Number of clusters K", y = "Mean silhouette width",
       title = sprintf("Cluster-count selection (best: K = %d, %s)",
                       best$K, best$method),
       colour = "Linkage") +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_silhouette%s.png", TAG)),
       p_sil, width = 7, height = 4, dpi = 140)

# 11d. Trace-variogram
tv_df <- tibble(distance = emp_dist, gamma = emp_gamma) |>
  filter(!is.na(distance), !is.na(gamma))
p_tv <- ggplot(tv_df, aes(distance, gamma)) +
  geom_point(alpha = 0.8, size = 2)
if (!is.null(fit_tv)) {
  d_grid <- seq(0, max(tv_df$distance), length.out = 200)
  sigma2 <- fit_tv$best$cov.pars[1]
  phi    <- fit_tv$best$cov.pars[2]
  nugget <- fit_tv$best$nugget
  cov.model <- fit_tv$best$cov.model
  cov_fn <- function(h) {
    switch(cov.model,
           "spherical"   = ifelse(h < phi,
                                  sigma2 * (1.5 * h/phi - 0.5 * (h/phi)^3),
                                  sigma2),
           "exponential" = sigma2 * (1 - exp(-h / phi)),
           "gaussian"    = sigma2 * (1 - exp(-(h / phi)^2)))
  }
  fit_df <- tibble(distance = d_grid, gamma = nugget + cov_fn(d_grid))
  p_tv <- p_tv + geom_line(data = fit_df, linewidth = 0.9, colour = "firebrick")
}
p_tv <- p_tv +
  labs(x = "Separation (deg)", y = expression(gamma(h)),
       title = "Empirical trace-variogram",
       subtitle = sprintf("Best model: %s, range = %.2f deg",
                          best_model, best_range)) +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_variogram%s.png", TAG)),
       p_tv, width = 7, height = 4, dpi = 140)

# --- 12. CSV ----------------------------------------------------------------
write_csv(cell_df |> dplyr::select(cell_id, lon, lat, cluster),
          file.path(OUT_DIR, sprintf("npp_bioregions%s.csv", TAG)))

message("Done. Wrote outputs to ", normalizePath(OUT_DIR))
