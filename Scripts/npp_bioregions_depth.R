# Depth-resolved multivariate spatial functional clustering.
# Annual-mean depth profiles of NPP, Chla, and Cphyto are fit with a
# B-spline basis on [0, 100] m, each variable's coefficients are
# z-scored across cells, then stacked into one super-vector per cell.
# Clustering follows the same Ballari sFDA recipe as npp_bioregions.R
# (trace-variogram + hclust + silhouette + bootstrap-ARI).
#
# Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/npp_bioregions_depth.R

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

NC_CBPM   <- "Output/cmems_cbpm_monthly_0p25deg.nc"
NC_BGC    <- "Output/cmems_chl_cphyto_monthly_0p25deg.nc"
OUT_DIR   <- "Output"
Z_GRID    <- seq(0, 100, by = 5)        # 21 levels
N_BASIS   <- 6                          # cubic B-spline coefficients
N_ORDER   <- 4
COV_FRAC  <- 7 / 9
K_RANGE   <- 2:8
N_BOOT    <- 100
DETREND   <- TRUE
FORCE_K   <- 5
TAG       <- if (!is.null(FORCE_K)) sprintf("_K%d", FORCE_K) else ""

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load both NetCDFs ---------------------------------------------------
message("Loading ", NC_CBPM)
nc1     <- nc_open(NC_CBPM)
lons    <- as.numeric(ncvar_get(nc1, "longitude"))
lats    <- as.numeric(ncvar_get(nc1, "latitude"))
z_cbpm  <- as.numeric(ncvar_get(nc1, "depth"))             # 0..199, 200 levels
t_units <- ncatt_get(nc1, "time", "units")$value
t_raw   <- ncvar_get(nc1, "time")
unit_word  <- sub(" since.*$", "", t_units)
origin_str <- sub("^[A-Za-z]+ since ", "", t_units)
origin     <- as.Date(substr(origin_str, 1, 10))
times <- switch(unit_word,
                "days"    = origin + t_raw,
                "hours"   = origin + t_raw / 24,
                "seconds" = origin + t_raw / 86400,
                stop("Unsupported time unit: ", unit_word))
npp_4d   <- ncvar_get(nc1, "npp")        # (lon, lat, depth, time)
npp_int  <- ncvar_get(nc1, "npp_int")    # (lon, lat, time)
nc_close(nc1)

message("Loading ", NC_BGC)
nc2     <- nc_open(NC_BGC)
z_bgc   <- as.numeric(ncvar_get(nc2, "depth"))             # ~22 native
chl_4d    <- ncvar_get(nc2, "chl")
cphyto_4d <- ncvar_get(nc2, "cphyto")
nc_close(nc2)

# Reorder to (time, depth, lat, lon)
npp_4d    <- aperm(npp_4d,    c(4, 3, 2, 1))
chl_4d    <- aperm(chl_4d,    c(4, 3, 2, 1))
cphyto_4d <- aperm(cphyto_4d, c(4, 3, 2, 1))
npp_int   <- aperm(npp_int,   c(3, 2, 1))

months <- as.numeric(format(times, "%m"))
years  <- as.numeric(format(times, "%Y"))
unique_years <- sort(unique(years))

message("Grid: ", length(times), " months x ", length(lats),
        " lat x ", length(lons), " lon  | depths cbpm=",
        length(z_cbpm), ", bgc=", length(z_bgc))

# --- 2. Coverage gate from npp_int climatology ------------------------------
clim_int <- array(NA_real_, dim = c(12, length(lats), length(lons)))
for (m in 1:12) {
  clim_int[m, , ] <- apply(npp_int[months == m, , , drop = FALSE],
                           c(2, 3), mean, na.rm = TRUE)
}
clim_int[!is.finite(clim_int)] <- NA_real_
prod_months   <- 2:10
n_valid_prod  <- apply(!is.na(clim_int[prod_months, , ]), c(2, 3), sum)
keep          <- n_valid_prod >= round(COV_FRAC * length(prod_months))
message("Cells passing npp_int coverage gate: ", sum(keep), " / ", length(keep))

# --- 3. Annual mean profiles + interpolate to Z_GRID ------------------------
mean_over_time <- function(arr4d, t_mask) {
  # arr4d: (time, depth, lat, lon)
  apply(arr4d[t_mask, , , , drop = FALSE], c(2, 3, 4), mean, na.rm = TRUE)
}

interp_profile <- function(prof, z_src, z_dst) {
  ok <- is.finite(prof)
  if (sum(ok) < 2) return(rep(NA_real_, length(z_dst)))
  approx(z_src[ok], prof[ok], xout = z_dst, rule = 2)$y
}

build_profiles <- function(t_mask) {
  # Returns list(npp, chl, cphyto): each a (length(Z_GRID), ncells) matrix
  npp_ann    <- mean_over_time(npp_4d, t_mask)        # (200, lat, lon)
  chl_ann    <- mean_over_time(chl_4d, t_mask)        # (22, lat, lon)
  cphyto_ann <- mean_over_time(cphyto_4d, t_mask)
  idx_mat <- which(keep, arr.ind = TRUE)
  ncells_local <- nrow(idx_mat)
  np <- matrix(NA_real_, nrow = length(Z_GRID), ncol = ncells_local)
  cl <- matrix(NA_real_, nrow = length(Z_GRID), ncol = ncells_local)
  cp <- matrix(NA_real_, nrow = length(Z_GRID), ncol = ncells_local)
  for (k in seq_len(ncells_local)) {
    li <- idx_mat[k, 1]; oi <- idx_mat[k, 2]
    np[, k] <- interp_profile(npp_ann[, li, oi],    z_cbpm, Z_GRID)
    cl[, k] <- interp_profile(chl_ann[, li, oi],    z_bgc,  Z_GRID)
    cp[, k] <- interp_profile(cphyto_ann[, li, oi], z_bgc,  Z_GRID)
  }
  list(npp = np, chl = cl, cphyto = cp, idx_mat = idx_mat)
}

message("Building annual-mean profiles on z = ", paste(range(Z_GRID), collapse = ".."),
        " m, step ", Z_GRID[2] - Z_GRID[1], " m")
prof_full <- build_profiles(rep(TRUE, length(times)))

# Drop cells with any NaN in any of the three profiles
finite_per_cell <- apply(is.finite(prof_full$npp), 2, all) &
                   apply(is.finite(prof_full$chl), 2, all) &
                   apply(is.finite(prof_full$cphyto), 2, all)
message("Cells fully finite across NPP/Chl/Cphyto profiles: ",
        sum(finite_per_cell), " / ", length(finite_per_cell))

# Re-index to the surviving cells
sel    <- which(finite_per_cell)
ncells <- length(sel)
prof_full$npp    <- prof_full$npp[, sel, drop = FALSE]
prof_full$chl    <- prof_full$chl[, sel, drop = FALSE]
prof_full$cphyto <- prof_full$cphyto[, sel, drop = FALSE]
prof_full$idx_mat <- prof_full$idx_mat[sel, , drop = FALSE]

cell_df <- tibble(
  cell_id = seq_len(ncells),
  lat_i = as.integer(prof_full$idx_mat[, 1]),
  lon_i = as.integer(prof_full$idx_mat[, 2]),
  lat   = lats[as.integer(prof_full$idx_mat[, 1])],
  lon   = lons[as.integer(prof_full$idx_mat[, 2])],
)

# --- 4. B-spline functional fit per variable --------------------------------
bb <- create.bspline.basis(rangeval = c(0, 100), nbasis = N_BASIS,
                           norder = N_ORDER)
fit_fd <- function(mat) Data2fd(argvals = Z_GRID, y = mat, basisobj = bb)
fd_npp    <- fit_fd(prof_full$npp)
fd_chl    <- fit_fd(prof_full$chl)
fd_cphyto <- fit_fd(prof_full$cphyto)
message("Spline basis: nbasis = ", N_BASIS, ", order = ", N_ORDER,
        ", range = [0, 100] m")

# --- 5. Z-score coefficients per (variable, basis-index) and stack ----------
zscore_rows <- function(coefs) {
  mu <- rowMeans(coefs)
  sdv <- apply(coefs, 1, sd)
  sdv[sdv == 0] <- 1
  (coefs - mu) / sdv
}
z_npp    <- zscore_rows(fd_npp$coefs)
z_chl    <- zscore_rows(fd_chl$coefs)
z_cphyto <- zscore_rows(fd_cphyto$coefs)
super_coef <- rbind(z_npp, z_chl, z_cphyto)              # (3*N_BASIS, ncells)
stopifnot(nrow(super_coef) == 3 * N_BASIS, ncol(super_coef) == ncells)
message("Super-coef matrix: ", nrow(super_coef), " x ", ncol(super_coef))

# Sanity: each row should have mean ~ 0, sd ~ 1.
row_mu <- rowMeans(super_coef)
row_sd <- apply(super_coef, 1, sd)
stopifnot(max(abs(row_mu)) < 1e-8, max(abs(row_sd - 1)) < 1e-8)

# --- 6. Spatial trend removal (coefficient-wise OLS on lat, lon) ------------
detrend_coefs <- function(coefs, lat_v, lon_v) {
  X <- cbind(1, lat_v, lon_v)
  XtX_inv <- solve(crossprod(X))
  fits <- X %*% XtX_inv %*% crossprod(X, t(coefs))
  coefs - t(fits)
}
if (DETREND) {
  super_coef <- detrend_coefs(super_coef, cell_df$lat, cell_df$lon)
  message("Detrending applied (coefficient-wise OLS on lat + lon).")
}

# --- 7. L2 distance (Euclidean on z-scored coefficients) --------------------
fast_l2sq <- function(coefs) {
  cross <- crossprod(coefs)                  # ncells x ncells
  a <- diag(cross)
  outer(a, a, "+") - 2 * cross
}
message("Computing L2 distance matrix (", ncells, " x ", ncells, ") ...")
L2sq <- fast_l2sq(super_coef)
L2sq[L2sq < 0] <- 0
diag(L2sq) <- 0
L2dist <- sqrt(L2sq)

# --- 8. Trace-variogram -----------------------------------------------------
coords <- as.matrix(cell_df[, c("lon", "lat")])
max_d  <- 0.6 * max(dist(coords))
message("Empirical trace-variogram, max.dist = ", sprintf("%.2f", max_d), " deg")
emp_tv <- trace.variog(coords = coords, L2norm = L2sq,
                       bin = TRUE, max.dist = max_d)
keep_bin <- is.finite(emp_tv$u) & is.finite(emp_tv$v)
emp_tv_clean <- list(u = emp_tv$u[keep_bin], v = emp_tv$v[keep_bin])
sill0 <- mean(tail(emp_tv_clean$v, max(3, length(emp_tv_clean$v) %/% 4)),
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

# --- 9. Hierarchical clustering + silhouette --------------------------------
dmat <- as.dist(L2dist)
hc_complete <- hclust(dmat, method = "complete")
hc_ward     <- hclust(dmat, method = "ward.D2")
sil_table <- tibble()
for (K in K_RANGE) {
  for (method in c("complete", "ward.D2")) {
    hc <- if (method == "complete") hc_complete else hc_ward
    cl <- cutree(hc, k = K)
    s  <- mean(silhouette(cl, dmat)[, "sil_width"])
    sil_table <- bind_rows(sil_table,
                           tibble(K = K, method = method, mean_sil = s))
  }
}
if (!is.null(FORCE_K)) {
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
  warning("Mean silhouette < 0.30 -- depth-bioregion partition is weak.")

# --- 10. Bootstrap stability over years -------------------------------------
message("Bootstrap stability (", N_BOOT, " iterations) ...")
ari_vec <- numeric(N_BOOT)
for (b in seq_len(N_BOOT)) {
  yrs    <- sample(unique_years, replace = TRUE)
  t_mask <- years %in% yrs
  prof_b <- build_profiles(t_mask)
  np <- prof_b$npp[,    sel, drop = FALSE]
  cl <- prof_b$chl[,    sel, drop = FALSE]
  cp <- prof_b$cphyto[, sel, drop = FALSE]
  if (any(!is.finite(np)) || any(!is.finite(cl)) || any(!is.finite(cp))) {
    ari_vec[b] <- NA; next
  }
  c_np <- zscore_rows(fit_fd(np)$coefs)
  c_cl <- zscore_rows(fit_fd(cl)$coefs)
  c_cp <- zscore_rows(fit_fd(cp)$coefs)
  sc_b <- rbind(c_np, c_cl, c_cp)
  if (DETREND) sc_b <- detrend_coefs(sc_b, cell_df$lat, cell_df$lon)
  L2_b <- fast_l2sq(sc_b)
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

# --- 10b. Cross-check vs. column-integrated bioregions ----------------------
ari_vs_int <- NA_real_
prior_csv <- "Output/npp_bioregions.csv"
if (file.exists(prior_csv)) {
  prior <- read_csv(prior_csv, show_col_types = FALSE)
  joined <- cell_df |>
    inner_join(prior |> dplyr::select(lon, lat, cluster_int = cluster),
               by = c("lon", "lat"))
  if (nrow(joined) > 50) {
    ari_vs_int <- adjustedRandIndex(joined$cluster, joined$cluster_int)
    message(sprintf("ARI vs column-integrated K=3 partition: %.3f (n=%d cells)",
                    ari_vs_int, nrow(joined)))
  }
}

writeLines(c(
  sprintf("K                : %d", best$K),
  sprintf("Linkage method   : %s", best$method),
  sprintf("Mean silhouette  : %.3f", best$mean_sil),
  sprintf("Variogram model  : %s (range = %.2f deg)", best_model, best_range),
  sprintf("Bootstrap ARI    : median = %.3f, 5%%-95%% = [%.3f, %.3f] (n=%d)",
          ari_med, ari_q05, ari_q95, N_BOOT),
  sprintf("ARI vs column-NPP: %.3f", ari_vs_int)
), file.path(OUT_DIR, sprintf("npp_bioregions_depth_stability%s.txt", TAG)))

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
       title = sprintf("Depth-resolved bioregions (K = %d, %s linkage)",
                       best$K, best$method),
       subtitle = "Annual-mean NPP / Chla / Cphyto profiles, 0-100 m",
       caption = sprintf("Mean silhouette = %.3f | Bootstrap ARI = %.2f",
                         best$mean_sil, ari_med)) +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_depth_clusters%s.png", TAG)),
       p_map, width = 9, height = 5, dpi = 140)

# 11b. Central profiles per cluster per variable (median + IQR)
prof_long <- bind_rows(
  as_tibble(prof_full$npp,    .name_repair = ~ as.character(seq_along(.x))) |>
    mutate(depth = Z_GRID, variable = "NPP (mg C m^-3 d^-1)"),
  as_tibble(prof_full$chl,    .name_repair = ~ as.character(seq_along(.x))) |>
    mutate(depth = Z_GRID, variable = "Chla (mg m^-3)"),
  as_tibble(prof_full$cphyto, .name_repair = ~ as.character(seq_along(.x))) |>
    mutate(depth = Z_GRID, variable = "Cphyto (mg C m^-3)")
) |>
  pivot_longer(-c(depth, variable), names_to = "cell_id", values_to = "val") |>
  mutate(cell_id = as.integer(cell_id)) |>
  left_join(cell_df |> dplyr::select(cell_id, cluster), by = "cell_id")

prof_summary <- prof_long |>
  group_by(variable, cluster, depth) |>
  summarise(med = median(val, na.rm = TRUE),
            q25 = quantile(val, 0.25, na.rm = TRUE),
            q75 = quantile(val, 0.75, na.rm = TRUE),
            .groups = "drop")

p_prof <- ggplot(prof_summary,
                 aes(x = med, y = depth, colour = factor(cluster),
                     fill = factor(cluster))) +
  geom_ribbon(aes(xmin = q25, xmax = q75), alpha = 0.20, colour = NA) +
  geom_path(linewidth = 0.9) +
  scale_y_reverse() +
  scale_colour_manual(values = clu_palette, name = "Bioregion") +
  scale_fill_manual(values = clu_palette, guide = "none") +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Depth (m)",
       title = "Central depth profile per bioregion",
       subtitle = "Median (line) and IQR (shaded) across cells in each cluster") +
  theme_bw()
ggsave(file.path(OUT_DIR,
                 sprintf("npp_bioregions_depth_central_profiles%s.png", TAG)),
       p_prof, width = 11, height = 5, dpi = 140)

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
ggsave(file.path(OUT_DIR,
                 sprintf("npp_bioregions_depth_silhouette%s.png", TAG)),
       p_sil, width = 7, height = 4, dpi = 140)

# 11d. Trace-variogram
tv_df <- tibble(distance = emp_tv_clean$u, gamma = emp_tv_clean$v) |>
  filter(!is.na(distance), !is.na(gamma))
p_tv <- ggplot(tv_df, aes(distance, gamma)) +
  geom_point(alpha = 0.8, size = 2)
if (!is.null(fit_tv)) {
  d_grid <- seq(0, max(tv_df$distance), length.out = 200)
  sigma2    <- fit_tv$best$cov.pars[1]
  phi       <- fit_tv$best$cov.pars[2]
  nugget    <- fit_tv$best$nugget
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
       title = "Empirical trace-variogram (z-scored super-coef space)",
       subtitle = sprintf("Best model: %s, range = %.2f deg",
                          best_model, best_range)) +
  theme_bw()
ggsave(file.path(OUT_DIR,
                 sprintf("npp_bioregions_depth_variogram%s.png", TAG)),
       p_tv, width = 7, height = 4, dpi = 140)

# --- 12. CSV ----------------------------------------------------------------
write_csv(cell_df |> dplyr::select(cell_id, lon, lat, cluster),
          file.path(OUT_DIR, sprintf("npp_bioregions_depth%s.csv", TAG)))

message("Done. Wrote outputs to ", normalizePath(OUT_DIR))
