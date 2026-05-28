# Bioregion clustering on column-integrated NPP using TWO functional
# inputs per cell:
#   (1) climatological 12-month curve (mean over years)
#   (2) per-month standard deviation across years (interannual variability)
#
# The two curves are fit on the same Fourier basis (period 12), each
# basis coefficient is z-scored across cells, the two coefficient
# blocks are stacked, and the rest of the Ballari sFDA pipeline is
# identical to Scripts/npp_bioregions.R.
#
# Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/npp_bioregions_var.R

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
K_RANGE   <- 2:8
N_BOOT    <- 100
DETREND   <- TRUE
FORCE_K   <- NULL
MONTH_MID <- seq(0.5, 11.5, by = 1)
TAG       <- if (!is.null(FORCE_K)) sprintf("_K%d", FORCE_K) else ""

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load ---------------------------------------------------------------
message("Loading ", NC_FILE)
nc      <- nc_open(NC_FILE)
lons    <- as.numeric(ncvar_get(nc, "longitude"))
lats    <- as.numeric(ncvar_get(nc, "latitude"))
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
npp_int <- ncvar_get(nc, "npp_int")     # (lon, lat, time)
nc_close(nc)
npp_int <- aperm(npp_int, c(3, 2, 1))   # -> (time, lat, lon)
months <- as.numeric(format(times, "%m"))
years  <- as.numeric(format(times, "%Y"))
unique_years <- sort(unique(years))
message("Grid: ", length(times), " months x ", length(lats),
        " lat x ", length(lons), " lon")

# --- 2. Build climatology AND per-month sd-across-years --------------------
build_two_curves <- function(t_mask) {
  # Returns list(clim, sd_curve): each (12 x lat x lon), winter -> 0/0.
  clim   <- array(NA_real_, dim = c(12, length(lats), length(lons)))
  sd_arr <- array(NA_real_, dim = c(12, length(lats), length(lons)))
  for (m in 1:12) {
    layer <- npp_int[t_mask & months == m, , , drop = FALSE]
    if (dim(layer)[1] == 0) next
    # per-year means within this month (collapse to 1 value per year per cell)
    yrs_in <- years[t_mask & months == m]
    yr_uniq <- unique(yrs_in)
    if (length(yr_uniq) < 2) {
      clim[m, , ] <- apply(layer, c(2, 3), mean, na.rm = TRUE)
      sd_arr[m, , ] <- 0
      next
    }
    per_yr <- array(NA_real_, dim = c(length(yr_uniq),
                                      length(lats), length(lons)))
    for (yi in seq_along(yr_uniq)) {
      per_yr[yi, , ] <- apply(layer[yrs_in == yr_uniq[yi], , , drop = FALSE],
                              c(2, 3), mean, na.rm = TRUE)
    }
    clim[m, , ]   <- apply(per_yr, c(2, 3), mean, na.rm = TRUE)
    sd_arr[m, , ] <- apply(per_yr, c(2, 3), sd,   na.rm = TRUE)
  }
  clim[!is.finite(clim)]     <- NA_real_
  sd_arr[!is.finite(sd_arr)] <- NA_real_
  list(clim = clim, sd_curve = sd_arr)
}

curves <- build_two_curves(rep(TRUE, length(times)))
clim     <- curves$clim
sd_curve <- curves$sd_curve

# --- 3. Coverage gate + winter padding -------------------------------------
prod_months   <- 2:10
n_valid_prod  <- apply(!is.na(clim[prod_months, , ]), c(2, 3), sum)
keep          <- n_valid_prod >= round(COV_FRAC * length(prod_months))
message("Cells passing coverage gate: ", sum(keep), " / ", length(keep))

# Winter (Jan, Nov, Dec) -> 0 for kept cells, in both curves
winter_months <- c(1, 11, 12)
for (m in winter_months) {
  layer <- clim[m, , ];     layer[keep & is.na(layer)]    <- 0
  clim[m, , ] <- layer
  layer <- sd_curve[m, , ]; layer[keep & is.na(layer)]    <- 0
  sd_curve[m, , ] <- layer
}
# Interp any remaining mid-year gaps for kept cells
for (i in seq_len(sum(keep))) {
  idx <- which(keep, arr.ind = TRUE)[i, ]
  cyc <- clim[, idx[1], idx[2]]
  if (any(is.na(cyc)))
    clim[, idx[1], idx[2]] <- approx(MONTH_MID, cyc, xout = MONTH_MID,
                                     rule = 2)$y
  cyc <- sd_curve[, idx[1], idx[2]]
  if (any(is.na(cyc)))
    sd_curve[, idx[1], idx[2]] <- approx(MONTH_MID, cyc, xout = MONTH_MID,
                                         rule = 2)$y
}

idx_mat <- which(keep, arr.ind = TRUE)
ncells  <- nrow(idx_mat)
cell_df <- tibble(
  cell_id = seq_len(ncells),
  lat_i = as.integer(idx_mat[, 1]),
  lon_i = as.integer(idx_mat[, 2]),
  lat   = lats[as.integer(idx_mat[, 1])],
  lon   = lons[as.integer(idx_mat[, 2])],
)
clim_mat <- matrix(NA_real_, nrow = 12, ncol = ncells)
sd_mat   <- matrix(NA_real_, nrow = 12, ncol = ncells)
for (k in seq_len(ncells)) {
  clim_mat[, k] <- clim[,     cell_df$lat_i[k], cell_df$lon_i[k]]
  sd_mat[, k]   <- sd_curve[, cell_df$lat_i[k], cell_df$lon_i[k]]
}
stopifnot(all(is.finite(clim_mat)), all(is.finite(sd_mat)))
message("Curve matrices: 12 months x ", ncells, " cells (clim + sd)")

# --- 4. Fourier basis fit per curve ----------------------------------------
fb <- create.fourier.basis(rangeval = c(0, 12), nbasis = N_BASIS, period = 12)
fd_clim <- Data2fd(argvals = MONTH_MID, y = clim_mat, basisobj = fb)
fd_sd   <- Data2fd(argvals = MONTH_MID, y = sd_mat,   basisobj = fb)

# --- 5. Z-score per (curve, basis-index) and stack -------------------------
zscore_rows <- function(coefs) {
  mu <- rowMeans(coefs)
  sdv <- apply(coefs, 1, sd)
  sdv[sdv == 0] <- 1
  (coefs - mu) / sdv
}
z_clim <- zscore_rows(fd_clim$coefs)
z_sd   <- zscore_rows(fd_sd$coefs)
super_coef <- rbind(z_clim, z_sd)    # (2 * N_BASIS, ncells)
stopifnot(nrow(super_coef) == 2 * N_BASIS)
message("Super-coef matrix: ", nrow(super_coef), " x ", ncol(super_coef))

# --- 6. Spatial trend removal ----------------------------------------------
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

# --- 7. L2 distance --------------------------------------------------------
fast_l2sq <- function(coefs) {
  cross <- crossprod(coefs)
  a <- diag(cross)
  outer(a, a, "+") - 2 * cross
}
message("Computing L2 distance matrix (", ncells, " x ", ncells, ") ...")
L2sq <- fast_l2sq(super_coef)
L2sq[L2sq < 0] <- 0
diag(L2sq) <- 0
L2dist <- sqrt(L2sq)

# --- 8. Trace-variogram ----------------------------------------------------
coords <- as.matrix(cell_df[, c("lon", "lat")])
max_d  <- 0.6 * max(dist(coords))
message("Empirical trace-variogram, max.dist = ", sprintf("%.2f", max_d), " deg")
emp_tv <- trace.variog(coords = coords, L2norm = L2sq,
                       bin = TRUE, max.dist = max_d)
keep_bin <- is.finite(emp_tv$u) & is.finite(emp_tv$v)
emp_tv_clean <- list(u = emp_tv$u[keep_bin], v = emp_tv$v[keep_bin])

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
} else {
  best <- sil_table |> arrange(desc(mean_sil)) |> slice(1)
}
message("Best silhouette: K = ", best$K, ", method = ", best$method,
        ", mean width = ", sprintf("%.3f", best$mean_sil))
hc_best   <- if (best$method == "complete") hc_complete else hc_ward
cluster_v <- cutree(hc_best, k = best$K)
cell_df$cluster <- cluster_v

# --- 10. Bootstrap stability over years ------------------------------------
message("Bootstrap stability (", N_BOOT, " iterations) ...")
ari_vec <- numeric(N_BOOT)
for (b in seq_len(N_BOOT)) {
  yrs    <- sample(unique_years, replace = TRUE)
  t_mask <- years %in% yrs
  cv_b <- build_two_curves(t_mask)
  clim_b <- cv_b$clim; sd_b <- cv_b$sd_curve
  clim_b[!is.finite(clim_b)] <- NA_real_
  sd_b[!is.finite(sd_b)] <- NA_real_
  for (m in winter_months) {
    layer <- clim_b[m, , ]; layer[keep & is.na(layer)] <- 0
    clim_b[m, , ] <- layer
    layer <- sd_b[m, , ];   layer[keep & is.na(layer)] <- 0
    sd_b[m, , ] <- layer
  }
  cmat_b <- matrix(NA_real_, nrow = 12, ncol = ncells)
  smat_b <- matrix(NA_real_, nrow = 12, ncol = ncells)
  for (k in seq_len(ncells)) {
    cyc <- clim_b[, cell_df$lat_i[k], cell_df$lon_i[k]]
    if (any(is.na(cyc))) cyc <- approx(MONTH_MID, cyc, xout = MONTH_MID, rule = 2)$y
    cmat_b[, k] <- cyc
    cyc <- sd_b[, cell_df$lat_i[k], cell_df$lon_i[k]]
    if (any(is.na(cyc))) cyc <- approx(MONTH_MID, cyc, xout = MONTH_MID, rule = 2)$y
    smat_b[, k] <- cyc
  }
  if (any(!is.finite(cmat_b)) || any(!is.finite(smat_b))) {
    ari_vec[b] <- NA; next
  }
  zc <- zscore_rows(Data2fd(MONTH_MID, cmat_b, fb)$coefs)
  zs <- zscore_rows(Data2fd(MONTH_MID, smat_b, fb)$coefs)
  sc_b <- rbind(zc, zs)
  if (DETREND) sc_b <- detrend_coefs(sc_b, cell_df$lat, cell_df$lon)
  L2_b <- fast_l2sq(sc_b); L2_b[L2_b < 0] <- 0
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

# --- 10b. Cross-check vs. previous partitions ------------------------------
ari_vs_clim <- NA_real_
prior_csv <- "Output/npp_bioregions.csv"
if (file.exists(prior_csv)) {
  prior <- read_csv(prior_csv, show_col_types = FALSE)
  joined <- cell_df |>
    inner_join(prior |> dplyr::select(lon, lat, cluster_int = cluster),
               by = c("lon", "lat"))
  if (nrow(joined) > 50) {
    ari_vs_clim <- adjustedRandIndex(joined$cluster, joined$cluster_int)
  }
}
message(sprintf("ARI vs climatology-only K=3 partition: %.3f", ari_vs_clim))

writeLines(c(
  sprintf("K                : %d", best$K),
  sprintf("Linkage method   : %s", best$method),
  sprintf("Mean silhouette  : %.3f", best$mean_sil),
  sprintf("Bootstrap ARI    : median = %.3f, 5%%-95%% = [%.3f, %.3f] (n=%d)",
          ari_med, ari_q05, ari_q95, N_BOOT),
  sprintf("ARI vs clim-only : %.3f", ari_vs_clim)
), file.path(OUT_DIR, sprintf("npp_bioregions_var_stability%s.txt", TAG)))

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
       title = sprintf("Bioregions on climatology + interannual variability (K = %d, %s)",
                       best$K, best$method),
       subtitle = "Two-curve sFDA: 12-month mean cycle + per-month sd-across-years",
       caption = sprintf("Mean silhouette = %.3f | Bootstrap ARI = %.2f",
                         best$mean_sil, ari_med)) +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_var_clusters%s.png", TAG)),
       p_map, width = 9, height = 5, dpi = 140)

# 11b. Central curves: two panels per cluster (mean cycle + sd cycle)
mk_long <- function(mat, label) {
  as_tibble(mat, .name_repair = ~ as.character(seq_along(.x))) |>
    mutate(month = MONTH_MID, curve = label) |>
    pivot_longer(-c(month, curve), names_to = "cell_id", values_to = "val") |>
    mutate(cell_id = as.integer(cell_id))
}
curves_long <- bind_rows(
  mk_long(clim_mat, "Climatology (mean cycle)"),
  mk_long(sd_mat,   "Interannual sd (variability)")
) |>
  left_join(cell_df |> dplyr::select(cell_id, cluster), by = "cell_id")

curve_summary <- curves_long |>
  group_by(curve, cluster, month) |>
  summarise(med = median(val), q25 = quantile(val, 0.25),
            q75 = quantile(val, 0.75), .groups = "drop")

month_lbl <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
p_curves <- ggplot(curve_summary,
                   aes(x = month, y = med,
                       colour = factor(cluster), fill = factor(cluster))) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = clu_palette, name = "Bioregion") +
  scale_fill_manual(values = clu_palette, guide = "none") +
  scale_x_continuous(breaks = 0:11 + 0.5, labels = month_lbl,
                     limits = c(0, 12)) +
  facet_grid(curve ~ cluster, scales = "free_y") +
  labs(x = NULL, y = expression("NPP (mg C m"^-2*" d"^-1*")"),
       title = "Central seasonal cycle and interannual sd per bioregion",
       subtitle = "Median (line) and IQR (shaded) across cells in each cluster") +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_var_central_curves%s.png", TAG)),
       p_curves, width = 12, height = 6, dpi = 140)

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
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_var_silhouette%s.png", TAG)),
       p_sil, width = 7, height = 4, dpi = 140)

# 11d. Variogram
tv_df <- tibble(distance = emp_tv_clean$u, gamma = emp_tv_clean$v) |>
  filter(!is.na(distance), !is.na(gamma))
p_tv <- ggplot(tv_df, aes(distance, gamma)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = "Separation (deg)", y = expression(gamma(h)),
       title = "Empirical trace-variogram (climatology + variability stack)") +
  theme_bw()
ggsave(file.path(OUT_DIR, sprintf("npp_bioregions_var_variogram%s.png", TAG)),
       p_tv, width = 7, height = 4, dpi = 140)

# --- 12. CSV ----------------------------------------------------------------
write_csv(cell_df |> dplyr::select(cell_id, lon, lat, cluster),
          file.path(OUT_DIR, sprintf("npp_bioregions_var%s.csv", TAG)))

message("Done. Wrote outputs to ", normalizePath(OUT_DIR))
