library(tidyverse)
source("Scripts/cbpm_r/cbpm_argo.r")
dat <- read_csv("Data/Processed/ALR6_fluorometer_data_cleaned.csv") |> select(-day_night)

#Regulrised every profile_num to have one obs per meter between 0 and 200m
dat_reg <- dat %>%
  mutate(depth = round(depth, 0)) |> 
  group_by(profile_num) %>%
  do({
    df <- .
    tibble(depth = seq(0, 200, by = 1)) %>%
      left_join(df, by = "depth") %>%
      group_by(depth) %>%
      summarise_all(median, na.omit = TRUE)
  }) %>%
  ungroup()


df_plot <- dat_reg |> filter(profile_num == 516)

ggplot(df_plot)+
  geom_path(aes(x = fluo, y = -depth))

#apply cbpm

dat_reg <- dat_reg |> mutate(year = lubridate::year(datetime),
                             month = lubridate::month(datetime),
                             day = lubridate::day(datetime))


dat_all_pp <- dat_reg |> 
  filter(!is.na(profile_num)) |> 
  group_by(profile_num) |> 
  mutate(
    carbon = bbp_baseline / (470/440) * 12128 + 0.59,
    izeu_cbpm = cbpm_argo(fluo_unquenched, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), median(na.omit(day)), mean(na.omit(lat)))$mzeu,
    cbpm_npp = cbpm_argo(fluo_unquenched, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), median(na.omit(day)), mean(na.omit(lat)))$pp,
    cbpm_mu = cbpm_argo(fluo_unquenched, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), median(na.omit(day)), mean(na.omit(lat)))$mu_z) |> 
  ungroup()


df_plot <- dat_all_pp |> filter(profile_num == 516)

ggplot(dat_all_pp)+
  geom_point(aes(x = profile_num, y = -depth, color = cbpm_npp))+
  scale_color_viridis_c()
  
dat_integrated <- dat_all_pp |> 
  group_by(profile_num) |> 
  summarise(
    npp_0_100m = sum(cbpm_npp[depth <= 100], na.rm = TRUE),
    npp_0_200m = sum(cbpm_npp[depth <= 200], na.rm = TRUE),
    mu_0_100m = mean(cbpm_mu[depth <= 100], na.rm = TRUE),
    mu_0_200m = mean(cbpm_mu[depth <= 200], na.rm = TRUE),
    izeu_0_100m = mean(izeu_cbpm[depth <= 100], na.rm = TRUE),
    izeu_0_200m = mean(izeu_cbpm[depth <= 200], na.rm = TRUE),
    lon = median(lon, na.rm = TRUE),
    lat = median(lat, na.rm = TRUE)
  ) |> 
  ungroup()

ggplot(dat_integrated)+
  geom_point(aes(x = lon, y = lat, color = npp_0_200m))+
  scale_color_viridis_c()+
  coord_quickmap()
