library(tidyverse)
library(zoo)

# ── Main NCP estimation function ──────────────────────────────────────────────

compute_ncp <- function(
    dat,
    time_step    = "15 days",
    lon_min      = -40, lon_max = -10,
    lat_min      =  58, lat_max =  65,
    date_start   = "2023-01-01",
    date_end     = "2024-01-01",
    zeu_default  = 40,
    mld_spar     = 0.3,
    label        = NULL        # label for plotting, defaults to time_step
) {
  
  if (is.null(label)) label <- time_step
  
  # ── Filter spatial and temporal extent ──────────────────────────────────────
  dat <- dat |>
    filter(
      lon >= lon_min & lon <= lon_max,
      lat >= lat_min & lat <= lat_max,
      date > as.Date(date_start),
      date < as.Date(date_end)
    )
  
  if (nrow(dat) == 0) stop("No data remaining after filtering for label: ", label)
  message(label, ": ", nrow(dat), " rows after spatial/temporal filter")
  
  # ── Euphotic depth ───────────────────────────────────────────────────────────
  dat <- dat |> mutate(zeu = zeu_default)
  
  # ── Per-profile MLD and Zeu ──────────────────────────────────────────────────
  dat_prof <- dat |>
    select(float_wmo, prof_number, date, MLD, zeu) |>
    distinct() |>
    arrange(date) |>
    filter(!is.na(MLD), !is.na(zeu))
  
  # ── Smooth MLD with loess ────────────────────────────────────────────────────
  dat_smooth <- dat_prof |> filter(!is.na(date), !is.na(MLD))
  
  fit_mld  <- loess(MLD ~ as.numeric(date), data = dat_smooth,
                    span = mld_spar, family = "symmetric")
  fit_zeu  <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$zeu, spar = 0.6)
  
  pred_mld <- pmax(predict(fit_mld, newdata = data.frame(date = as.numeric(dat_prof$date))), 0)
  pred_zeu <- predict(fit_zeu, as.numeric(dat_prof$date))$y
  
  dat_prof <- dat_prof |>
    mutate(mld_smooth = pred_mld, zeu_smooth = pred_zeu)
  
  dat <- left_join(dat,
                   select(dat_prof, prof_number, float_wmo, mld_smooth, zeu_smooth),
                   by = c("float_wmo", "prof_number"))
  
  # ── Integration depth on time grid ──────────────────────────────────────────
  time_grid <- tibble(
    date_grid = seq(min(dat_prof$date), max(dat_prof$date), by = time_step)
  )
  
  integration_depth_data <- dat |>
    group_by(float_wmo, prof_number, date) |>
    summarise(mld_smooth = unique(mld_smooth),
              zeu_smooth = unique(zeu_smooth),
              .groups = "drop")
  
  prof_dat_smoothed <- integration_depth_data |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid) |>
    summarise(mld = mean(mld_smooth, na.rm = TRUE),
              zeu = mean(zeu_smooth, na.rm = TRUE),
              .groups = "drop") |>
    mutate(
      prev_MLD              = lag(mld),
      prev_zeu              = lag(zeu),
      NCP_integration_depth = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
      next_integration_depth = lead(NCP_integration_depth)
    ) |>
    full_join(time_grid, by = "date_grid") |>
    arrange(date_grid) |>
    mutate(
      NCP_integration_depth  = zoo::na.approx(NCP_integration_depth,
                                              x = as.numeric(date_grid),
                                              na.rm = FALSE, rule = 2),
      next_integration_depth = zoo::na.approx(next_integration_depth,
                                              x = as.numeric(date_grid),
                                              na.rm = FALSE, rule = 2)
    ) |>
    filter(date_grid %in% time_grid$date_grid) |>
    mutate(
      dt = as.numeric(date_grid - lag(date_grid)),
      NCP_integration_depth_smooth  = rollapply(NCP_integration_depth,
                                                width = 3, FUN = mean,
                                                align = "center", fill = NA),
      next_integration_depth_smooth = rollapply(next_integration_depth,
                                                width = 3, FUN = mean,
                                                align = "center", fill = NA)
    ) |>
    na.omit()
  
  # ── Nitrate interpolation on time grid ──────────────────────────────────────
  dat_smoothed <- dat |>
    mutate(date_grid = as.Date(cut(date, breaks = time_step))) |>
    group_by(date_grid, depth) |>
    summarise(nitrate = mean(canyon_nitrate, na.rm = TRUE), .groups = "drop") |>
    full_join(time_grid, by = "date_grid") |>
    arrange(date_grid) |>
    group_by(depth) |>
    mutate(nitrate = zoo::na.approx(nitrate, x = as.numeric(date_grid),
                                    na.rm = FALSE, rule = 2)) |>
    ungroup() |>
    filter(date_grid %in% time_grid$date_grid) |>
    left_join(prof_dat_smoothed, by = "date_grid") |>
    na.omit()
  
  # ── Nitrate integration ──────────────────────────────────────────────────────
  integrate_nitrate <- function(nitrate, depth, from, to) {
    idx <- depth >= from & depth <= to
    if (sum(idx) < 2) return(NA_real_)
    approxfun(depth[idx], nitrate[idx])(seq(from, to, by = 1)) |>
      mean(na.rm = TRUE) * (to - from)
  }
  
  dat_final <- dat_smoothed |>
    group_by(date_grid, mld, zeu) |>
    filter(!is.na(next_integration_depth)) |>
    summarise(
      int_N_mmol_m2       = integrate_nitrate(nitrate, depth, 0,
                                              unique(NCP_integration_depth)),
      next_int_N_mmol_m2  = integrate_nitrate(nitrate, depth, 0,
                                              unique(next_integration_depth)),
      mld_concentration   = integrate_nitrate(nitrate, depth,
                                              unique(NCP_integration_depth) - 20,
                                              unique(NCP_integration_depth) - 10) / 10,
      sub_mld_concentration = integrate_nitrate(nitrate, depth,
                                                unique(NCP_integration_depth) + 20,
                                                unique(NCP_integration_depth) + 30) / 10,
      .groups = "drop"
    )
  
  # ── NCP computation ──────────────────────────────────────────────────────────
  dat_final_smoothed <- dat_final |>
    mutate(
      int_N_smooth      = rollapply(int_N_mmol_m2,      width = 3, FUN = mean,
                                    align = "center", fill = NA),
      next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean,
                                    align = "center", fill = NA)
    )
  
  ncp_results <- dat_final_smoothed |>
    arrange(date_grid) |>
    mutate(
      dt             = as.numeric(date_grid - lag(date_grid)),
      diff_mld       = mld - lag(mld),
      we             = pmax(0, diff_mld),
      delta          = sub_mld_concentration - mld_concentration,
      diff           = delta * we,
      nitrate_consumption = lag(next_int_N_smooth) - int_N_smooth,
      c_consumption  = nitrate_consumption * 6.625 / dt,
      NCP            = c_consumption + diff
    ) |>
    mutate(
      NCP       = case_when(NCP < -60 ~ NCP / 2, TRUE ~ NCP),
      datenum   = as.numeric(date_grid)
    ) |>
    filter(!is.na(NCP))
  
  fit_loess <- loess(NCP ~ datenum, data = ncp_results, span = 0.4, family = "symmetric")
  ncp_results$loess_ncp <- predict(fit_loess)
  
  # tag with label for multi-run plotting
  ncp_results$time_step_label <- label
  
  return(ncp_results)
}


# ── Run multiple timesteps and overlay ────────────────────────────────────────

plot_ncp_comparison <- function(dat, time_steps = c("10 days", "15 days", "30 days"), ...) {
  
  results <- map(time_steps, function(ts) {
    tryCatch(
      compute_ncp(dat, time_step = ts, label = ts, ...),
      error = function(e) { message("Failed for ", ts, ": ", e$message); NULL }
    )
  }) |>
    compact() |>
    bind_rows()
  
  ggplot(results) +
    geom_line(aes(x = date_grid, y = NCP, color = time_step_label), alpha = 0.4) +
    #geom_line(aes(x = date_grid, y = loess_ncp, color = time_step_label), linewidth = 1.2) +
    scale_color_brewer(palette = "Set1", name = "Time step") +
    labs(y = "mmol C m⁻² d⁻¹", x = "Date",
         title = "NCP estimates — sensitivity to time step") +
    theme_bw() +
    scale_x_date(date_labels = "%b %Y", breaks = "2 months")
}

# single run
ncp_15d <- compute_ncp(dat, time_step = "15 days")

# sensitivity analysis — three timesteps overlaid
p <- plot_ncp_comparison(dat, time_steps = c("5 days", "10 days", "15 days", "20 days"))
print(p)
ggsave("output/ncp_sensitivity_timestep.png", p, width = 10, height = 5)

# also easy to vary other parameters
ncp_strict_mld <- compute_ncp(dat, time_step = "15 days", mld_spar = 0.5, label = "strict MLD")
ncp_loose_mld  <- compute_ncp(dat, time_step = "15 days", mld_spar = 0.2, label = "loose MLD")

p2 <- bind_rows(ncp_strict_mld, ncp_loose_mld) |>
  ggplot() +
  geom_line(aes(x = date_grid, y = loess_ncp, color = time_step_label), linewidth = 1.2) +
  scale_color_brewer(palette = "Set1", name = "MLD smoothing") +
  labs(y = "mmol C m⁻² d⁻¹") +
  theme_bw()

print(p2)
