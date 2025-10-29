# ============================================================
#   Compute NCP from nitrate drawdown using 3-point regression
# ============================================================

# --- Libraries ---
library(tidyverse)
library(lubridate)
library(castr)
library(zoo)
library(slider)

# --- Load data ---
dat <- read_csv("output/float_nitrate_data_corrected.csv")

# --- 1. Smooth MLD & Zeu time series ------------------------
dat_prof <- dat |>
  select(float_wmo, prof_number, date, MLD, zeu) |>
  distinct() |>
  arrange(date)

dat_smooth <- dat_prof %>%
  filter(!is.na(date), !is.na(MLD))

fit_mld <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$MLD, spar = 0.6)
fit_zeu <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$zeu, spar = 0.6)

dat_prof <- dat_prof %>%
  mutate(
    mld_smooth = predict(fit_mld, as.numeric(date))$y,
    zeu_smooth = predict(fit_zeu, as.numeric(date))$y
  )

dat <- left_join(dat, select(dat_prof, float_wmo, prof_number, mld_smooth, zeu_smooth),
                 by = c("float_wmo", "prof_number"))

# --- 2. Integration depth computation ------------------------
integration_depth_data <- dat |>
  group_by(float_wmo, prof_number, date) |>
  summarise(
    mld_smooth = unique(mld_smooth, na.rm = TRUE),
    zeu_smooth = unique(zeu_smooth, na.rm = TRUE)
  ) |>
  ungroup()

# Regrid on a 12-day interval (approx 10-day)
time_grid <- tibble(date_10day = seq(min(integration_depth_data$date),
                                     max(integration_depth_data$date),
                                     by = "12 days"))

prof_dat_smoothed <- integration_depth_data %>%
  mutate(date_10day = as.Date(cut(date, breaks = "12 days"))) %>%
  group_by(date_10day) %>%
  summarise(
    mld = mean(mld_smooth, na.rm = TRUE),
    zeu = mean(zeu_smooth, na.rm = TRUE)
  ) %>%
  mutate(
    prev_MLD = lag(mld),
    prev_zeu = lag(zeu),
    NCP_integration_depth = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
    next_integration_depth = lead(NCP_integration_depth)
  ) %>%
  full_join(time_grid, by = "date_10day") %>%
  arrange(date_10day) %>%
  mutate(
    NCP_integration_depth = na.approx(NCP_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
    next_integration_depth = na.approx(next_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)
  ) %>%
  filter(date_10day %in% time_grid$date_10day) %>%
  mutate(
    dt = as.numeric(date_10day - lag(date_10day)),
    NCP_integration_depth_smooth = rollapply(NCP_integration_depth, width = 3, FUN = mean, align = "center", fill = NA),
    next_integration_depth_smooth = rollapply(next_integration_depth, width = 3, FUN = mean, align = "center", fill = NA)
  ) %>%
  na.omit()

# --- 3. Interpolate nitrate & NPP ---------------------------------------
dat_smoothed <- dat %>%
  mutate(date_10day = as.Date(cut(date, breaks = "12 days"))) %>%
  group_by(date_10day, depth) %>%
  summarise(
    nitrate = mean(nitrate_corrected, na.rm = TRUE),
    npp = mean(cbpm_npp, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  full_join(time_grid, by = "date_10day") %>%
  arrange(date_10day) %>%
  group_by(depth) %>%
  mutate(
    nitrate = na.approx(nitrate, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
    npp = na.approx(npp, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)
  ) %>%
  ungroup() %>%
  filter(date_10day %in% time_grid$date_10day)

dat_smoothed <- dat_smoothed %>%
  left_join(prof_dat_smoothed, by = "date_10day") %>%
  na.omit()

# --- 4. Integrate nitrate over depth ------------------------------------
dat_final <- dat_smoothed %>%
  group_by(date_10day, mld, zeu) %>%
  filter(!is.na(next_integration_depth)) %>%
  summarise(
    int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(NCP_integration_depth)),
    next_int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(next_integration_depth)),
    npp = integrate(npp, depth, from = 0, to = 200),
    mld_concentration = integrate(nitrate, depth,
                                  from = unique(NCP_integration_depth) - 20,
                                  to = unique(NCP_integration_depth) - 10) / 10,
    sub_mld_concentration = integrate(nitrate, depth,
                                      from = unique(NCP_integration_depth) + 20,
                                      to = unique(NCP_integration_depth) + 30) / 10
  ) %>%
  ungroup()

dat_final_smoothed <- dat_final %>%
  mutate(
    int_N_smooth = rollapply(int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA),
    next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA)
  )

# --- 5. Compute NCP using 3-point regression -----------------------------
ncp_results <- dat_final_smoothed %>%
  arrange(date_10day) %>%
  mutate(
    dt = as.numeric(date_10day - lag(date_10day)),
    nitrate_slope = slide_dbl(
      .x = cur_data(),
      .before = 1, .after = 1,
      .f = ~{
        if (nrow(.x) < 3 || any(is.na(.x$int_N_smooth))) return(NA_real_)
        mod <- lm(int_N_smooth ~ as.numeric(date_10day), data = .x)
        coef(mod)[2]  # slope (mmol N m⁻² per day)
      }
    ),
    # Convert nitrate slope to C units (mg C m⁻² d⁻¹)
    c_consumption = nitrate_slope * 6.625,
    NCP = -c_consumption,  # negative slope = production
    # Physical correction for entrainment
    diff_mld = mld - lag(mld),
    we = pmax(0, diff_mld),
    delta = sub_mld_concentration - mld_concentration,
    diff = delta * we,
    NCP_total = NCP + diff
  ) %>%
  mutate(
    ncp_smooth = rollapply(NCP_total, width = 6, FUN = mean, align = "center", fill = NA),
    ncp_sd = rollapply(NCP_total, width = 6, FUN = sd, align = "center", fill = NA),
    npp_smooth = rollapply(npp, width = 6, FUN = median, align = "center", fill = NA)
  ) %>%
  ungroup()

# --- 6. Plot results -----------------------------------------------------
ggplot(ncp_results) +
  geom_line(aes(x = date_10day, y = npp_smooth/12, color = "NPP")) +
  geom_ribbon(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd),
              alpha = 0.3, fill = "grey70") +
  geom_line(aes(x = date_10day, y = ncp_smooth, color = "NCP")) +
  labs(y = "mg C m⁻² d⁻¹", color = NULL) +
  theme_bw()
