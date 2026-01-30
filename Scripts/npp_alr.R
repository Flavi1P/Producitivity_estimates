library(tidyverse)
library(matrixStats)
source("Scripts/cbpm_r/cbpm_argo.r")
dat <- read_csv("Data/Processed/ALR6_fluorometer_data_cleaned.csv") |> select(-day_night)
dat2 <- read_csv("Data/Processed/final_csv_alr6.csv")

#Regulrised every profile_num to have one obs per meter between 0 and 200m
dat_reg <- dat2 %>%
  mutate(depth = round(depth, 0)) |> 
  group_by(profile_id) %>%
  do({
    df <- .
    tibble(depth = seq(0, 200, by = 1)) %>%
      left_join(df, by = "depth") %>%
      group_by(depth) %>%
      summarise_all(median, na.omit = TRUE)
  }) %>%
  ungroup()


df_plot <- dat_reg |> filter(profile_id == unique(dat2$profile_id)[25])

ggplot(df_plot)+
  geom_path(aes(x = fluo, y = -depth))

#apply cbpm

dat_reg <- dat_reg |> mutate(year = lubridate::year(date),
                             month = lubridate::month(date),
                             day = lubridate::day(date))

sat_data <- dat %>% select(datetime, satellite_par) %>% 
  mutate(dateonly = as.Date(datetime)) %>% 
  select(-datetime) %>% 
  distinct()

shelf_par <- data.frame("dateonly" = as.Date("2024-06-12"), "satellite_par" = 22)
sat_data = bind_rows(sat_data, shelf_par)
dat_reg <- dat_reg %>% mutate(dateonly = as.Date(date)) %>% 
  left_join(sat_data)

dat_all_pp <- dat_reg |> 
  filter(!is.na(profile_id)) |> 
  group_by(profile_id) |> 
  mutate(
    carbon = bbp / (470/440) * 12128 + 0.59,
    izeu_cbpm = cbpm_argo(FLUO_NPQ, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), median(na.omit(day)), mean(na.omit(lat)))$mzeu,
    cbpm_npp = cbpm_argo(FLUO_NPQ, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), median(na.omit(day)), mean(na.omit(lat)))$pp,
    cbpm_mu = cbpm_argo(FLUO_NPQ, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), median(na.omit(day)), mean(na.omit(lat)))$mu_z) |> 
  ungroup()


df_plot <- dat_all_pp |> filter(profile_id == unique(dat2$profile_id)[25])

ggplot(dat_all_pp)+
  geom_point(aes(x = date, y = -depth, color = cbpm_npp))+
  scale_color_viridis_c()
  
dat_integrated <- dat_all_pp |> 
  group_by(profile_id) |> 
  summarise(
    npp_0_100m = sum(cbpm_npp[depth <= 100], na.rm = TRUE),
    npp_0_200m = sum(cbpm_npp[depth <= 200], na.rm = TRUE),
    mu_0_100m = mean(cbpm_mu[depth <= 100], na.rm = TRUE),
    mu_0_200m = mean(cbpm_mu[depth <= 200], na.rm = TRUE),
    izeu_0_100m = mean(izeu_cbpm[depth <= 100], na.rm = TRUE),
    izeu_0_200m = mean(izeu_cbpm[depth <= 200], na.rm = TRUE),
    lon = median(lon, na.rm = TRUE),
    lat = median(lat, na.rm = TRUE),
    date = mean(date)
  ) |> 
  ungroup()

library(patchwork)

p1 <- ggplot(arrange(dat_integrated, date))+
  geom_path(aes(x = lon, y = lat, color = npp_0_200m/12), linewidth = 3)+
  scale_color_viridis_c(name = "Integrated NPP (mmol Cm-1")+
  coord_quickmap()+
  theme_minimal()

p2 <- ggplot(arrange(dat_integrated, date))+
  geom_line(aes(x = date, y = npp_0_200m/12), linewidth = 1)+
  scale_color_viridis_c()+
  theme_minimal()+
  ylab("Integrated NPP (mmol c m-2 d-1)")

p3 <- ggplot(dat_all_pp)+
  geom_point(aes(x = date, y = -depth, color = cbpm_npp), size = 3)+
  scale_color_viridis_c()+
  theme_minimal()+
  ylab("Depth (m)")+
  xlab("Date")+
  ylim(-50, 0)

p1 + (p2 / p3)

ggsave("output/cbpm_output_zoomed.png", width = 18, height = 10, dpi = 300)

ggplot(dat_all_pp %>% filter(depth_round < 25))+
  geom_boxplot(aes(x = as.factor(cluster), y = cbpm_npp, fill = as.factor(cluster)))+
  scale_color_brewer(palette = "Set1")+
  theme_minimal(base_size = 16)

ggsave("output/boxplot_npp.png", width = 18, height = 10, dpi = 300)
