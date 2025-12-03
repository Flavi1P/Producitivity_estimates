library(tidyverse)
library(lubridate)
library(castr)
library(zoo)

dat <- read_csv("Data/Processed/float_nitrate_data_corrected.csv")


# filter the dataset to stay with ICB -------------------------------------
#add iceland contour on map

ggplot(select(dat, lon, lat) |> distinct())+
  geom_point(aes(x = lon, y = lat))+
  geom_rect(aes(xmin = -45, xmax = -10, ymin = 58, ymax = 65), fill = NA, color = "red")+
  theme_minimal()+#add land contour
  borders("world", xlim = c(-60, 0), ylim = c(50, 70), fill = "lightgrey", color = "black")+
  coord_quickmap(xlim = c(-55, - 10), ylim = c(50, 65))

ggsave("Output/float_locations_icb.png", width = 8, height = 6)
dat <- filter(dat, lon >= -40 & lon <= -10 & lat >= 58 & lat <= 65)
range(dat$date)
#filter to keep data between 2019 and now
dat <- filter(dat, date >= as.Date("2019-01-01"))
# MLD and Zeu time series -------------------------------------------------

#Artificially set zeu as 40 if NA
dat <- dat |> 
   mutate(zeu = case_when(is.na(zeu) ~ 40,
                          TRUE ~ zeu))
#Create a dataset with MLD and Zeu for each profile 
dat_prof <- dat |> select(float_wmo, prof_number, date, MLD, zeu) |> 
  distinct() |> 
  arrange(date) |> 
  filter(!is.na(MLD) & !is.na(zeu))

#Smooth the MLD and Zeu time series using a spline
# Remove NAs before smoothing
dat_smooth <- dat_prof %>%
  filter(!is.na(date), !is.na(MLD), !is.na(zeu))

# Fit spline
fit_mld <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$MLD, spar = 0.6)
fit_zeu <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$zeu, spar = 0.6)
# Predict smoothed values for all dates
pred_mld <- predict(fit_mld, as.numeric(dat_prof$date))$y
pred_zeu <- predict(fit_zeu, as.numeric(dat_prof$date))$y


# Add back to main data frame
dat_prof <- dat_prof %>%
  mutate(mld_smooth = pred_mld,
         zeu_smooth = pred_zeu)


ggplot(dat_prof)+
  geom_line(aes(x = date, y = -mld_smooth, color = "MLD"))+
  geom_point(aes(x = date, y = -mld_smooth), shape = 1)+
  geom_line(aes(x = date, y = -zeu_smooth, color = "Zeu"))+
  labs(y = "Depth (m)", color = "Parameter")+
  theme_minimal()+
  xlim(as.Date("2024-01-01"), as.Date("2025-12-31"))
#ggsave("Output/mld_zeu_time_series_icb3.png", width = 10, height = 6)
  
dat <- left_join(dat, select(dat_prof, prof_number, float_wmo, mld_smooth, zeu_smooth), by = c("float_wmo", "prof_number"))


# Integration depth computation -------------------------------------------


#Define Integration depth (i.e. deepest value of MLD or Zeu between t and t-1)
integration_depth_data <- dat |>
  group_by(float_wmo, prof_number, date) |>
  summarise(
    mld_smooth = unique(mld_smooth, na.rm = TRUE),
    zeu_smooth = unique(zeu_smooth, na.rm = TRUE)) |> 
  ungroup()

#Regrid the data on a 10day interval timeseries
time_grid <- tibble(date_10day = seq(min(integration_depth_data$date),
                                     max(integration_depth_data$date),
                                     by = "15 days"))

prof_dat_smoothed <- integration_depth_data %>%
  mutate(date_10day = as.Date(cut(date, breaks = "15 days"))) %>%
  group_by(date_10day) %>%
  summarise(mld = mean(mld_smooth, na.rm = TRUE),
            zeu = mean(zeu_smooth, na.rm = TRUE)) |> 
  mutate(
    prev_MLD = dplyr::lag(mld),
    prev_zeu = dplyr::lag(zeu),
    NCP_integration_depth = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
    next_integration_depth = lead(NCP_integration_depth)
  ) |> 
  full_join(time_grid, by = "date_10day") |>
  arrange(date_10day) %>%
  mutate(
    NCP_integration_depth = zoo::na.approx(NCP_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
    next_integration_depth = zoo::na.approx(next_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)
  ) |> 
  mutate(dt = as.numeric(date_10day - dplyr::lag(date_10day)),
         NCP_integration_depth_smooth = rollapply(NCP_integration_depth, width = 3, FUN = mean, align = "center", fill = NA),
         next_integration_depth_smooth = rollapply(next_integration_depth, width = 3, FUN = mean, align = "center", fill = NA)
  )

ggplot(filter(prof_dat_smoothed, !is.na(mld)))+
  geom_line(aes(x = date_10day, y = -NCP_integration_depth_smooth))+
  geom_line(aes(x = date_10day, y = -next_integration_depth_smooth), linetype = "dashed")

#ggsave("Output/integration_depth_time_series_icb2.png", width = 10, height = 6)
prof_dat_smoothed <- na.omit(prof_dat_smoothed)


# Interpolation of data  --------------------------------------------------

#Create a dataset with the median value of nitrate integrals every 10 days for smoothing
dat_smoothed <- dat %>%
  mutate(date_10day = as.Date(cut(date, breaks = "15 days"))) %>%
  group_by(date_10day, depth) %>%
  summarise(
    nitrate = mean(nitrate_corrected, na.rm = TRUE),#canbe set either as nitrate_corrected or canyon_nitrate
    npp = mean(cbpm_npp, na.rm = TRUE),
    chla = mean(chla, na.rm = TRUE)) |> 
  ungroup() |> 
  full_join(time_grid, by = "date_10day") %>%     # ensure regular 10-day spacing
  arrange(date_10day) %>%
  group_by(depth) |> 
  mutate(
    nitrate = na.approx(nitrate, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
    npp = na.approx(npp, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)) %>%
  ungroup() |> 
  filter(date_10day %in% time_grid$date_10day)




ggplot(dat_smoothed)+
  geom_point(aes(x = nitrate, y = - depth, color = date_10day))
#ggsave("Output/nitrate_profiles_icb.png", width = 8, height = 6)
dat_smoothed <- dat_smoothed |> 
  left_join(prof_dat_smoothed, by = "date_10day") |> 
  na.omit()

ggplot(dat_smoothed)+
  geom_tile(aes(x = date_10day, y = - depth, fill = nitrate))+
  geom_line(aes(x = date_10day, y = - NCP_integration_depth), color = "white", size = 1)+
  scale_fill_viridis_c()+
  labs(y = "Depth (m)", x = "Date", fill = "Nitrate (µmol/kg)")+
  theme_minimal()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")


ggsave("Output/synthetic nitrate transect.png", width = 10, height = 6)
# Integration of nitrate over integration depth ---------------------------

dat_final <- dat_smoothed |>
  group_by(date_10day, mld, zeu) |>
  filter(!is.na(next_integration_depth)) |> 
  summarise(
    int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(NCP_integration_depth)),
    next_int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(next_integration_depth)),
    npp = integrate(npp, depth, from = 0, to = 200),
    mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) - 20, to = unique(NCP_integration_depth) - 10) / 10,
    sub_mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) + 20, to = unique(NCP_integration_depth) + 30) / 10,
    chla_integrated = integrate(chla, depth, from = 0, to = unique(NCP_integration_depth)),
    chla_mld = integrate(chla, depth, from = 0, to = unique(mld)),
    chla_surf = integrate(chla, depth, from = 0, to = 10)
  ) |> 
  ungroup()

ggplot(dat_final)+
  geom_line(aes(x = date_10day, y = int_N_mmol_m2))+
  geom_point(aes(x = date_10day, y = int_N_mmol_m2))+
  geom_line(aes(x = date_10day, y = next_int_N_mmol_m2), linetype = "dashed")+
  labs(y = "Integrated Nitrate (mmol N m-2)")+
  theme_minimal()+
  scale_x_date(date_labels = "%b %Y", breaks = "1 months")

#ggsave("Output/integrated_nitrate_time_series_icb2.png", width = 10, height = 6)

ggplot(dat_final)+
  geom_line(aes(x = date_10day, y = chla_integrated, color = "Integration depth"))+
  geom_point(aes(x = date_10day, y = chla_mld, color = "MLD"))+
  geom_line(aes(x = date_10day, y = chla_surf, color = "Surface (10m)"))+
  labs(y = "Integrated Chla (mg m-2)")+
  theme_minimal()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

#ggsave("Output/integrated_chla_time_series_icb2.png", width = 10, height = 6)


#Smoothing the integration
dat_final_smoothed <- dat_final %>%
  mutate(
    int_N_smooth = rollapply(int_N_mmol_m2, width = 3, FUN = mean, align = "left", fill = NA),
    next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean, align = "left", fill = NA)
  )

ggplot(dat_final_smoothed)+
  geom_line(aes(x = date_10day, y = int_N_mmol_m2), linetype = "dashed")+
  geom_line(aes(x = date_10day, y = int_N_smooth))+
  labs(y = "Integrated Nitrate (mmol N m-2)")+
  theme_minimal()+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  ggtitle("Smoothed Integrated Nitrate Time Series ICB")

#ggsave("Output/integrated_nitrate_smooth_time_series_icb2.png", width = 10, height = 6)
# Computing NCP -----------------------------------------------------------


ncp_results <- dat_final_smoothed |> 
  arrange(date_10day) |> 
  mutate(dt = as.numeric(date_10day - dplyr::lag(date_10day)),
         diff_mld = mld - lag(mld),
         we = pmax(0, diff_mld),
         delta = sub_mld_concentration - mld_concentration,
         diff = delta * we,
         nitrate_consumption =  dplyr::lag(next_int_N_mmol_m2) - int_N_mmol_m2,
         nitrate_consumption_smooth =  dplyr::lag(next_int_N_smooth) - int_N_smooth,
         c_consumption = nitrate_consumption_smooth * 6.625 / dt,
         NCP = c_consumption + diff) |> 
  ungroup()

ncp_results <- ncp_results |> 
  mutate(ncp_smooth = rollapply(NCP, width = 3, FUN = mean, align = "center", fill = NA),
         ncp_sd = rollapply(NCP, width = 3, FUN = sd, align = "center", fill = NA),
         npp_smooth = rollapply(npp, width = 3, FUN = median, align = "center", fill = NA),
         npp_sd = rollapply(npp, width = 3, FUN = sd, align = "center", fill = NA))

ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = NCP, color = "NCP"))+
  geom_point(aes(x = date_10day, y = NCP))+
  geom_path(aes(x = date_10day, y = npp/12, color = "NPP"))+
  theme_bw()+
  labs(y = "mmol C m-2 d-1")+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  xlab("Date")

ggsave("Output/ncp_time_series_icb_raw.png", width = 10, height = 6)
  

ggplot(filter(ncp_results, npp != 0))+
  geom_line(aes(x = date_10day, y = npp_smooth/12, color = "NPP"))+
  geom_errorbar(aes(x = date_10day, y = npp_smooth/12, ymin = npp_smooth/12 - npp_sd/12, ymax = npp_smooth/12 + npp_sd/12), alpha = 0.3)+
  geom_errorbar(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd), alpha = 0.3)+
  geom_line(aes(x = date_10day, y = ncp_smooth, color = "NCP"))+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  xlim(as.Date("2024-01-01"), as.Date("2025-12-01"))+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  xlab("Date")

ggsave("Output/ncp_time_series_icb_focus_2024_2025.png", width = 10, height = 6)

ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = npp_smooth/12, color = "NPP"))+
  geom_errorbar(aes(x = date_10day, y = npp_smooth/12, ymin = npp_smooth/12 - npp_sd/12, ymax = npp_smooth/12 + npp_sd/12), alpha = 0.3)+
  geom_errorbar(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd), alpha = 0.3)+
  geom_line(aes(x = date_10day, y = ncp_smooth, color = "NCP"))+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  xlab("Date") 

ggsave("Output/ncp_time_series_icb_full_length_ncp.png", width = 10, height = 6)

synthetic_ncp <- ncp_results |>
  filter(date_10day > as.Date("2024-01-01")) |>
  mutate(month_of_year = month(date_10day),
         npp = npp/12) |> 
  group_by(month_of_year) |> 
  summarise(mean_ncp = mean(NCP, na.rm = TRUE),
            sd_ncp = sd(NCP, na.rm = TRUE),
            mean_npp = mean(npp, na.rm = TRUE),
            sd_npp = sd(npp, na.rm = TRUE))

ggplot(synthetic_ncp)+
  geom_line(aes(x = month_of_year, y = mean_ncp, color = "NCP"))+
  geom_ribbon(aes(x = month_of_year, ymin = mean_ncp - sd_ncp, ymax = mean_ncp + sd_ncp), alpha = 0.2)+
  geom_line(aes(x = month_of_year, y = mean_npp, color = "NPP"))+
  geom_ribbon(aes(x = month_of_year, ymin = mean_npp - sd_npp, ymax = mean_npp + sd_npp), alpha = 0.2)+
  labs(x = "Month", y = "mmol C m-2 d-1", color = "Variable")+
  scale_x_continuous(breaks = 1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  theme_bw()+
  xlim(1,11.5)
ggsave("Output/mean_synthetic_ncp_npp_icb_clean.png", width = 10, height = 6)

ncp_results_bp <- ncp_results |> 
  filter(date_10day > as.Date("2024-01-01")) |>
  mutate(npp = npp/12,
         month_of_year = month(date_10day))
ggplot(ncp_results_bp)+
  geom_boxplot(aes(x = as.factor(month_of_year), y = npp, color = "NPP"))+
  geom_boxplot(aes(x = as.factor(month_of_year), y = NCP, color = "NCP"))+
  theme_bw()+
  scale_color_brewer(palette = "Set1")+
  xlab("Month")+
  scale_x_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  ylab("mmol C m-2 d-1")+
  ylim(-60, 150)


ggsave("Output/boxplot_ts_ncp_npp.png", width = 10, height = 6)

ncp_results <- ncp_results |> 
  mutate(status = case_when(ncp_smooth > 0 ~ "autotrophic",
                            ncp_smooth < 0 ~ "heterotrophic",
                            TRUE ~ "balanced"))

ggplot(filter(ncp_results, npp != 0))+
  geom_point(aes(x = date_10day, y = ncp_smooth, color = status))+
  geom_line(aes(x = date_10day, y = npp_smooth/12), linetype = "dashed")+
  labs(y = "mmol C m-2 d-1")+
  scale_color_manual(values = c("autotrophic" = "#a6d854", "heterotrophic" = "#fb8072"))+
  theme_minimal()+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  xlab("Date")

#ggsave("Output/ncp_status_time_series_icb3.png", width = 10, height = 6)


ggplot(synthetic_ncp|> arrange(month_of_year))+
  geom_path(aes(x = mean_npp, y = mean_ncp, color = month_of_year))+
  geom_point(aes(x = mean_npp, y = mean_ncp, color = month_of_year))+
  scale_color_viridis_c(option = "C", direction = -1)+
  theme_minimal()+
  xlab("Mean NPP (mmol C m-2 d-1)")+
  ylab("Mean NCP (mmol C m-2 d-1)")

ncp_results <- ncp_results |>
  mutate(trophic_index = ncp_smooth/npp_smooth)

ggplot(filter(ncp_results, npp != 0))+
  geom_line(aes(x = date_10day, y = trophic_index))+
  geom_point(aes(x = date_10day, y = trophic_index, color = status))+
  geom_line(aes(x = date_10day, y = chla_surf/100, color = "Surface chla"), data = filter(dat_final, npp !=0))+
  geom_line(aes(x = date_10day, y = chla_integrated/100, color = "Integrated chla"), data = filter(dat_final, npp !=0))+
  scale_color_manual(values = c("autotrophic" = "blue", "heterotrophic" = "red", "balanced" = "grey", "Integrated chla" = "green", "Surface chla" = "darkgreen"))+
  theme_minimal()+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  xlab("Date")+
  ylab("Trophic Index (NCP/NPP)")

#ggsave("Output/trophic_index_time_series_icb_with_chla2.png", width = 10, height = 6)

ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = -mld))+
  geom_point(aes(x = date_10day, y = -mld, color = diff))+
  scale_color_viridis_c()

ggplot(ncp_results)+
  geom_point(aes(x = date_10day, y = - mld, color = NCP))+
  scale_color_distiller(palette = "RdBu")+
  theme_dark()

# draft -------------------------------------------------------------------

# # Ensure data are ordered by profile and date
# dat <- dat |> arrange(prof_number, date, depth)
# 
# # Step 1: Compute integrated nitrate down to each possible depth limit
# integrate_to_depth <- function(df, depth_limit) {
#   sub <- df |> filter(depth <= depth_limit)
#   if (nrow(sub) < 2) return(NA)  # Need at least 2 points for integration
#   trapz(sub$depth, sub$nitrate)
# }
# 
# # Step 2: Compute integrated nitrate to 200 m for sanity check
# dat_integrated_200 <- dat |>
#   group_by(prof_number, date, MLD) |>
#   summarise(int_N_200_mmol_m2 = integrate_to_depth(cur_data(), 200), .groups = "drop")
# 
# # Step 3: Compute NCP integration depth dynamically
# dat_integrated <- dat_integrated_200 |>
#   arrange(date) |>
#   mutate(
#     prev_MLD = lag(MLD),
#     NCP_integration_depth = pmax(MLD, prev_MLD, 40, na.rm = TRUE),
#     next_integration_depth = lead(NCP_integration_depth)
#   )
# 
# dat_joined <- left_join(dat, dat_integrated)
# 
# dat_final <- dat_joined |>
#   group_by(prof_number, date, MLD, NCP_integration_depth, next_integration_depth) |>
#   summarise(
#     int_N_mmol_m2 = integrate_to_depth(pick(nitrate, depth), unique(NCP_integration_depth)),
#     next_int_N_mmol_m2 = integrate_to_depth(pick(nitrate, depth), unique(next_integration_depth)),
#     .groups = "drop"
#   )
# 
# dat_final <- dat_final |> 
#   mutate(nitrate_consumption = int_N_mmol_m2 - lag(next_int_N_mmol_m2),
#          c_consumption = nitrate_consumption * 6.625 / 10,
#          NCP = - c_consumption * 12)  # Convert to mmol C m-2
# 
# ggplot(filter(dat_final, NCP > -4000))+
#   geom_line(aes(x = date, y = NCP))
