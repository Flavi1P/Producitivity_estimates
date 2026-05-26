library(tidyverse)

dat <- read_csv("~/ncp_from_canyon/output/NorthAtlantic_20s/ncp/ncp_uncertainty.csv")

ggplot(dat)+
  geom_point(aes(x = date_grid, y = NCP_mean))+
  geom_errorbar(aes(x = date_grid, ymin = NCP_mean - NCP_sd, ymax = NCP_mean + NCP_sd))+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95), alpha = 0.1)+
  geom_line(aes(x = date_grid, y = tot))+
  theme_bw()

dat$tot <- cumsum(dat$NCP_mean)

library(arrow)

dat <- read_parquet("Data/Processed/glider_npp/all_gliders_npp.parquet")
doombar

dat_int <- 
  dat %>% 
  select(date, glider, depth, profile_index, npp, npp_vgpm) %>%
  group_by(date, profile_index, glider) %>% 
  summarise(npp_cbpm = castr::integrate(npp, depth, from =0, to = 200),
            npp_vgpm = mean(npp_vgpm, na.rm = TRUE),
            vgpm_check = sd(npp_vgpm, na.rm = TRUE)) %>% 
  ungroup() %>%
  select(-profile_index) %>% 
  group_by(date, glider) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup()

ggplot(dat_int)+
  geom_line(aes(x = date, y = npp_vgpm, color = glider))


ggplot(dat_int)+
  geom_line(aes(x = date, y = npp_cbpm, color = glider))

dat_avg <- dat_int %>% 
  select(-glider) %>% 
  group_by(date) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup()

ggplot(dat_avg)+
  geom_line(aes(x = date, y = npp_cbpm/12, color= "cbpm"), linewidth = 1)+
  geom_line(aes(x = date, y = npp_vgpm/12, color = "vgpm"), linewidth = 1)+
  geom_line(aes(x = date, y = npp_cbpm/12, group = glider, color = "cbpm"), alpha = 0.3, data = dat_int)+
  geom_line(aes(x = date, y = npp_vgpm/12, group = glider, color = "vgpm"), alpha = 0.3, data = dat_int)+
  ylab("NPP (mmol C-1 m-2 d-1)")+
  scale_color_discrete(name = "Model")+
  theme_bw()+
  ylim(0, 220)
ggsave("Output/glider_npp_plot.png", width = 16, height = 9, dpi = 300)

write_csv(dat_int, "Data/Processed/glider_npp/biocarbon_integrated_npp_per_glider.csv")
write_csv(dat_int, "Data/Processed/glider_npp/biocarbon_integrated_npp_averaged_glider.csv")

northatl <- read_csv("~/ncp_from_canyon/output/NorthAtlantic_1box/ncp/IcelandBasin/ncp_uncertainty.csv")

northatl <- northatl %>% mutate(region = "northatl")

ggplot()+
  geom_line(aes(x = date_grid, y = NCP_mean, color = "Iceland Basin"), dat = filter(northatl, time_step_label == "15 days"))+
  scale_color_brewer(palette = "Set1", name = "Basin")+
  ylab("NCP (mmol C m-2 d-1)")+
  theme_bw(base_size = 18)+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

ggsave("Output/NCP_NorthAtlantic_15days.png", width = 16, height = 9, dpi = 300)

ggplot()+
  geom_line(aes(x = date_grid, y = NCP_mean, color = time_step_label), dat = northatl)+
  scale_color_brewer(palette = "Set2", name = "Time bin")+
  ylab("NCP (mmol C m-2 d-1)")+
  theme_bw(base_size = 18)+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

ggsave("Output/NCP_NorthAtlantic_timewindow.png", width = 16, height = 9, dpi = 300)

ggplot(northatl)+
  geom_line(aes(x = date_grid, y = NCP_mean, color = time_step_label))+
  geom_errorbar(aes(x = date_grid, ymin = NCP_mean - NCP_sd, ymax = NCP_mean +NCP_sd, group = time_step_label), color = "Grey")+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = time_step_label), alpha = 0.2)+
  scale_color_brewer(palette = "Set2", name = "Time bin")+
  scale_fill_brewer(palette = "Set2", name = "Q05-Q95")+
  ylab("NCP (mmol C m-2 d-1)")+
  theme_bw(base_size = 18)+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  facet_wrap(.~time_step_label, ncol = 1, scales = "free_y")

ggsave("Output/NCP_NorthAtlantic_montecarlo.png", width = 16, height = 9, dpi = 300)

ancp_na <- northatl %>% 
  mutate(month = lubridate::month(date_grid),
         year = lubridate::year(date_grid),
         timestep = as.numeric(str_extract(time_step_label, "[0-9]{2}")),
         tot_ncp = NCP_mean * timestep,
         tot_q05 = NCP_q05 * timestep,
         tot_q95 = NCP_q95 * timestep,
         tot_sd = NCP_sd * timestep) %>%
  group_by(year, time_step_label) %>% 
  summarise(ANCP = sum(tot_ncp)/1000,
            ANCP_q05 = sum(tot_q05)/1000,
            ANCP_q95 = sum(tot_q95)/1000,
            ANCP_sd = sqrt(sum(tot_sd**2))/1000,
            conf_int = ANCP_sd * 1.96) %>% 
  ungroup()
ggplot(ancp_na)+
  geom_col(aes(x = year, y = ANCP, fill = time_step_label), position = "dodge")+
  #geom_errorbar(aes(x = year, ymin = ANCP_q05/1000, ymax = ANCP_q95/1000, group = time_step_label), position = "dodge")+
  geom_errorbar(aes(x = year, ymin = ANCP-ANCP_sd, ymax = ANCP + ANCP_sd, group = time_step_label), position = "dodge")+
  scale_fill_brewer(palette = "Set2", name = "Time bin")+
  ylab("ANCP (mol m-2 yr-1)")+
  theme_bw(base_size = 18)
ggsave("Output/ANCP_sd.png", width = 16, height = 9, dpi = 300)

ggplot(ancp_na)+
  geom_col(aes(x = year, y = ANCP, fill = time_step_label), position = "dodge")+
  #geom_errorbar(aes(x = year, ymin = ANCP_q05/1000, ymax = ANCP_q95/1000, group = time_step_label), position = "dodge")+
  geom_errorbar(aes(x = year, ymin = ANCP-conf_int, ymax = ANCP + conf_int, group = time_step_label), position = "dodge")+
  scale_fill_brewer(palette = "Set2", name = "Time bin")+
  ylab("ANCP (mol m-2 yr-1)")+
  theme_bw(base_size = 18)

ggsave("Output/ANCP_conf_int.png", width = 16, height = 9, dpi = 300)


icb <- read_csv("~/ncp_from_canyon/output/Single_float_test/ncp/IcelandBasin/ncp_uncertainty.csv")
irm <- read_csv("~/ncp_from_canyon/output/NorthAtlantic_mld_methods_comparison/ncp/IrmingerSea/ncp_uncertainty.csv")
lbr <- read_csv("~/ncp_from_canyon/output/NorthAtlantic_seas_comparison/ncp/LabradorSea/ncp_uncertainty.csv")
nws <- read_csv("~/ncp_from_canyon/output/NorthAtlantic_seas_comparison/ncp/NorwegianSea/ncp_uncertainty.csv")
float <- read_csv("~/ncp_from_canyon/output/Single_float_test/ncp_float/1900722/ncp_float.csv")

icb_rate <- read_csv("~/ncp_from_canyon/output/Single_float_test/ncp/IcelandBasin/ncp_results.csv")

icb <- icb %>% mutate(region = "icb")
irm  <- irm %>% mutate(region = "irm")
lbr  <- lbr %>% mutate(region = "lbr")
nws <- nws %>% mutate(region = "nws")


grouped_regions <- bind_rows(icb, irm, lbr, nws)


# Formatting and computing climatologies, ancp, timeseries.... ------------

ancp_basins <- grouped_regions %>% 
  group_by(region) %>% 
  mutate(month = lubridate::month(date_grid),
         year = lubridate::year(date_grid),
         timestep = as.numeric(str_extract(time_step_label, "[0-9]{2}")),
         tot_ncp = NCP_mean * timestep,
         tot_q05 = NCP_q05 * timestep,
         tot_q95 = NCP_q95 * timestep,
         tot_sd = NCP_sd * timestep) %>%
  ungroup() %>% 
  group_by(region, year, time_step_label) %>% 
  summarise(ANCP = sum(tot_ncp)/1000,
            ANCP_q05 = sum(tot_q05)/1000,
            ANCP_q95 = sum(tot_q95)/1000,
            ANCP_sd = sqrt(sum(tot_sd**2))/1000,
            conf_int = ANCP_sd * 1.96) %>% 
  ungroup()

clima_float <-   float %>% 
  mutate(NCP = case_when(NCP > 0 ~ NCP,
                         NCP <= 0 ~0)) %>% 
  mutate(month = lubridate::month(date),
         year = lubridate::year(date)) %>% 
  group_by(month, year) %>% 
  summarise(NCP = mean(NCP),
            NCP_sd = sd(NCP)) %>% 
  ungroup()

ANCP_float <- clima_float %>% 
  mutate(tot_ncp = NCP * 30.5,
         tot_sd = NCP_sd * 30.5) %>% 
  group_by(year) %>% 
  summarise(ancp = sum(tot_ncp)/1000,
            asd = sum(tot_sd)/1000)

climatology <- 
  grouped_regions %>% 
  mutate(month = lubridate::month(date_grid),
         year = lubridate::year(date_grid)) %>% 
  group_by(region, month, time_step_label) %>% 
  summarise(NCP = mean(NCP_mean),
            sd = mean(NCP_sd),
            q05 = mean(NCP_q05),
            q95 = mean(NCP_q95)) %>% 
  ungroup()
  

ANCP_clima <- climatology %>% 
  filter(time_step_label == "30 days") %>% 
  group_by(region) %>% 
  mutate(tot_ncp = NCP * 30.5,
         tot_q05 = q05 * 30.5,
         tot_q95 = q95 * 30.5,
         sd = sd * 30.5) %>% 
  summarise(ancp = sum(tot_ncp)/1000,
            aq05 = sum(tot_q05)/1000,
            aq95 = sum(tot_q95)/1000,
            sd = sqrt(sum(sd**2))/1000) %>% 
  ungroup()

ANCP_clima_timestep_positive <- climatology %>% 
  group_by(region, time_step_label) %>% 
  mutate(NCP = case_when(NCP > 0 ~ NCP,
                         NCP <= 0 ~0),
         q95 = case_when(q95 > 0 ~ q95,
                         q95 <= 0 ~ 0),
         q05 = case_when(q05 > 0 ~ q05,
                         q05 <= 0 ~ 0),
         tot_ncp = NCP * 30.5,
         tot_q05 = q05 * 30.5,
         tot_q95 = q95 * 30.5,
         sd = sd * 30.5) %>% 
  summarise(ancp = sum(tot_ncp)/1000,
            aq05 = sum(tot_q05)/1000,
            aq95 = sum(tot_q95)/1000,
            sd = sqrt(sum(sd**2))/1000) %>% 
  ungroup()

ANCP_clima_timestep <- climatology %>% 
  group_by(region, time_step_label) %>% 
  mutate(tot_ncp = NCP * 30.5,
         tot_q05 = q05 * 30.5,
         tot_q95 = q95 * 30.5,
         sd = sd * 30.5) %>% 
  summarise(ancp = sum(tot_ncp)/1000,
            aq05 = sum(tot_q05)/1000,
            aq95 = sum(tot_q95)/1000,
            sd = sqrt(sum(sd**2))/1000) %>% 
  ungroup()

rate_clima <- icb_rate %>% 
  mutate(rate = NCP/mld,
         month = lubridate::month(date_grid),
         year = lubridate::year(date_grid)) %>% 
  group_by(month) %>% 
  summarise(rate_mean = mean(rate),
            rate_sd = sd(rate)) %>% 
  ungroup()

# plotting just two regions -----------------------------------------------


ggplot(filter(grouped_regions, region %in% c("icb", "irm") & time_step_label == "30 days"))+
  geom_line(aes(x = date_grid, y = NCP_mean, color = region))+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = region), alpha = 0.2)+
  scale_fill_brewer(palette= "Set1", name = "Monte Carlo q05 - q95", guide = "none")+
  scale_color_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  ggtitle("30 days time window comparison")

ggsave("Output/ANCP_two_basin_30days.png", width = 16, height = 9, dpi = 300)

ggplot(filter(grouped_regions, region %in% c("icb", "irm")))+
  geom_line(aes(x = date_grid, y = NCP_mean, color = time_step_label))+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = time_step_label), alpha = 0.2)+
  scale_fill_brewer(palette= "Set2", name = "Monte Carlo q05 - q95", guide = "none")+
  scale_color_brewer(palette = "Set2", name = "Time bin")+
  theme_bw(base_size = 18)+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  facet_wrap(.~region, ncol  = 1, scales = "free_y")
ggsave("Output/ANCP_two_basin_timestep.png", width = 16, height = 9, dpi = 300)


ggplot(filter(ancp_basins, region %in% c("icb", "irm")))+
  geom_col(aes(x = year, y = ANCP, fill = time_step_label), position = "dodge")+
  geom_errorbar(aes(x = year, ymin = ANCP-conf_int, ymax = ANCP + conf_int, group = time_step_label), position = "dodge")+
  scale_fill_brewer(palette = "Set2", name = "Time bin")+
  ylab("ANCP (mol m-2 yr-1)")+
  theme_bw(base_size = 18)+
  facet_wrap(.~region, ncol = 1)
ggsave("Output/ANCP_barplot_two_basin_timestep.png", width = 16, height = 9, dpi = 300)


ggplot(filter(climatology, region %in% c("icb", "irm")))+
  geom_line(aes(x = month, y = NCP, color = time_step_label))+
  geom_ribbon(aes(x = month, ymin = q05, ymax = q95, fill = time_step_label), alpha = 0.2)+
  scale_color_brewer(palette = "Set2", name = "Time step")+
  scale_fill_brewer(palette = "Set2", guide = "none")+
  theme_bw(base_size = 18)+
  facet_wrap(.~ region, ncol = 1)
ggsave("Output/climatology_two_basin_timestep.png", width = 16, height = 9, dpi = 300)

ggplot(filter(climatology, region %in% c("icb", "irm") & time_step_label == "30 days"))+
  geom_line(aes(x = month, y = NCP, color = region))+
  geom_ribbon(aes(x = month, ymin = q05, ymax = q95, fill = region), alpha = 0.2)+
  scale_color_brewer(palette = "Set1", name = "Basin")+
  scale_fill_brewer(palette = "Set1", guide = "none")+
  ylab("NCP (mmol C m-2 d-1)")+
  theme_bw(base_size = 18)

ggsave("Output/climatology_two_basin_30days.png", width = 16, height = 9, dpi = 300)


ggplot(filter(ANCP_clima, region %in% c("icb", "irm")))+
  geom_col(aes(x = region, y = ancp, fill = region), alpha = 0.8)+
  geom_errorbar(aes(x = region, ymin = ancp -sd, ymax = ancp + sd))+
  scale_fill_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  ylab("ANCP (mol c m-2 yr-1)")+
  ggtitle("Iceland Basin = 0.73 mol C m-2 yr -1 \n Irminger Sea = 0.36 mol C m-2 yr-1")
ggsave("Output/climatology_ANCP_with_errobar.png", width = 16, height = 9, dpi = 300)

ggplot(filter(ANCP_clima, region %in% c("icb", "irm")))+
  geom_col(aes(x = region, y = ancp, fill = region), alpha = 0.8)+
  scale_fill_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  ylab("ANCP (mol c m-2 yr-1)")+
  ggtitle("Iceland Basin = 0.73 mol C m-2 yr -1 \n Irminger Sea = 0.36 mol C m-2 yr-1")
ggsave("Output/climatology_ANCP_without_errobar.png", width = 16, height = 9, dpi = 300)

ggplot(filter(ANCP_clima_timestep, region %in% c("icb", "irm")))+
  geom_col(aes(x = region, y = ancp, fill = time_step_label), alpha = 0.8, position = "dodge")+
  scale_fill_brewer(palette = "Set2", name = "Time bin")+
  theme_bw(base_size = 18)+
  ylab("ANCP (mol c m-2 yr-1)")+
  ggtitle("Iceland Basin = 3.51 mol C m-2 yr -1 \n Irminger Sea = 2.61 mol C m-2 yr-1")

ggplot(filter(ANCP_clima_timestep_positive, region %in% c("icb", "irm")))+
    geom_col(aes(x = region, y = ancp, fill = time_step_label), alpha = 0.8, position = "dodge")+
    geom_errorbar(aes(x = region, ymin = aq05, ymax = aq95, group = time_step_label), alpha = 0.8, position = "dodge")+
    scale_fill_brewer(palette = "Set2", name = "Time bin")+
    theme_bw(base_size = 18)+
    ylab("ANCP (mol c m-2 yr-1)")

ggplot(rate_clima)+
  geom_line(aes(x = month, y = rate_mean))+
  geom_ribbon(aes(x = month, ymin = rate_mean - rate_sd, ymax = rate_mean + rate_sd), alpha = 0.4)+
  ylab("mmol C m-3 d-1")+
  theme_bw()
# plotting float ----------------------------------------------------------

ggplot(filter(ANCP_float, year != 2006 & year != 2010))+
  geom_col(aes(x = year, y = ancp))

mean(filter(ANCP_float, year != 2006 & year != 2010) %>% pull(ancp))
sd(filter(ANCP_float, year != 2006 & year != 2010) %>% pull(ancp))

# Plotting four regions ---------------------------------------------------

ggplot(filter(grouped_regions, time_step_label == "30 days"))+
  geom_line(aes(x = date_grid, y = NCP_mean, color = region))+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = region), alpha = 0.2)+
  scale_fill_brewer(palette= "Set1", name = "Monte Carlo q05 - q95", guide = "none")+
  scale_color_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  ggtitle("30 days time window comparison")

ggsave("Output/ANCP_four_basins_30days.png", width = 16, height = 9, dpi = 300)

ggplot(grouped_regions)+
  geom_line(aes(x = date_grid, y = NCP_mean, color = time_step_label))+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = time_step_label), alpha = 0.2)+
  scale_fill_brewer(palette= "Set2", name = "Monte Carlo q05 - q95", guide = "none")+
  scale_color_brewer(palette = "Set2", name = "Time bin")+
  theme_bw(base_size = 18)+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  facet_wrap(.~region, ncol  = 1, scales = "free_y")
ggsave("Output/ANCP_four_basins_timestep.png", width = 16, height = 9, dpi = 300)


ggplot(ancp_basins)+
  geom_col(aes(x = year, y = ANCP, fill = time_step_label), position = "dodge")+
  geom_errorbar(aes(x = year, ymin = ANCP-conf_int, ymax = ANCP + conf_int, group = time_step_label), position = "dodge")+
  scale_fill_brewer(palette = "Set2", name = "Time bin")+
  ylab("ANCP (mol m-2 yr-1)")+
  theme_bw(base_size = 18)+
  facet_wrap(.~region, ncol = 1)
ggsave("Output/ANCP_barplot_four_basin_timestep.png", width = 16, height = 9, dpi = 300)


ggplot(climatology)+
  geom_line(aes(x = month, y = NCP, color = time_step_label))+
  geom_ribbon(aes(x = month, ymin = q05, ymax = q95, fill = time_step_label), alpha = 0.2)+
  scale_color_brewer(palette = "Set2", name = "Time step")+
  scale_fill_brewer(palette = "Set2", guide = "none")+
  theme_bw(base_size = 18)+
  facet_wrap(.~ region, ncol = 1)
ggsave("Output/climatology_four_basins_timestep.png", width = 16, height = 9, dpi = 300)

ggplot(filter(climatology, time_step_label == "30 days"))+
  geom_line(aes(x = month, y = NCP, color = region))+
  geom_ribbon(aes(x = month, ymin = q05, ymax = q95, fill = region), alpha = 0.2)+
  scale_color_brewer(palette = "Set1", name = "Basin")+
  scale_fill_brewer(palette = "Set1", guide = "none")+
  ylab("NCP (mmol C m-2 d-1)")+
  theme_bw(base_size = 18)

ggsave("Output/climatology_four_basins_30days.png", width = 16, height = 9, dpi = 300)


ggplot(ANCP_clima)+
  geom_col(aes(x = region, y = ancp, fill = region), alpha = 0.8)+
  geom_errorbar(aes(x = region, ymin = ancp -sd, ymax = ancp + sd))+
  scale_fill_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  ylab("ANCP (mol c m-2 yr-1)")+
  ggtitle("Iceland Basin: 0.73 mol C m-2 yr -1 \n Irminger Sea: 0.36 \n Labrador sea: 0.14 \n Norwegian Sea: -0.03")
ggsave("Output/climatology_ANCP_four_basins_with_errobar.png", width = 16, height = 9, dpi = 300)

ggplot(ANCP_clima)+
  geom_col(aes(x = region, y = ancp, fill = region), alpha = 0.8)+
  scale_fill_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  ylab("ANCP (mol c m-2 yr-1)")+
  ggtitle("Iceland Basin: 0.73 mol C m-2 yr -1 \n Irminger Sea: 0.36 \n Labrador sea: 0.14 \n Norwegian Sea: -0.03")
ggsave("Output/climatology_ANCP_four_basins_without_errobar.png", width = 16, height = 9, dpi = 300)

ggplot(ANCP_clima_timestep)+
  geom_col(aes(x = region, y = ancp, fill = time_step_label), position = "dodge", alpha = 0.8)+
  scale_fill_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  ylab("ANCP (mol c m-2 yr-1)")

literature <- data.frame(region = c("Quay et al 2012 (Subpolar NA)", "Jeansson et al 2015 (Iceland sea)", "Yoder et al 2025 (Irminger)", "Mork et al 2024 (Norwegian sea)"), ancp = c(2.8, 7.5, 9.8, 4))
ancp_tot <- bind_rows(ANCP_clima, literature)

ggplot(ancp_tot)+
  geom_col(aes(x = region, y = ancp, fill = region), alpha = 0.8)+
  scale_fill_brewer(palette = "Set1", name = "Basin")+
  theme_bw(base_size = 18)+
  ylab("ANCP (mol c m-2 yr-1)")+
  ggtitle("Iceland Basin: 0.73 mol C m-2 yr -1 \n Irminger Sea: 0.36 \n Labrador sea: 0.14 \n Norwegian Sea: -0.03")

ggsave("Output/ANCP_literature.png", width = 16, height = 9, dpi = 300)


# looking at DOC climatology ----------------------------------------------

DOC <- read_csv("C:/Users/petit/Downloads/111995.csv")

doc <- filter(DOC, lon < -12 & lon > -32 & lat > 55 & lat < 63) %>% 
  group_by(season) %>% 
  summarise_all(mean) %>% 
  mutate(rate = diff(c(surf_doc_avg[4], surf_doc_avg))) %>% 
  mutate(season = case_when(season == 1 ~ 1,
                   season == 2 ~ 4,
                   season == 3 ~ 7,
                   season == 4 ~ 10))

ggplot(doc)+
  geom_col(aes(x = season, y = rate/90))+
  geom_line(aes(x = month, y = rate_mean), dat = rate_clima)+
  ylim(-0.15, 0.15)

# Old plots ---------------------------------------------------------------
literature = data.frame(
  source = c("Alkire_14", "Alkire_14", "Alkire_12", "Alkire_12", "Kortzinger_08a", "Kortzinger_08b", "Quay_12", "Quay_12"),
  month = c(5, 6, 4, 5, 6, 6, 1, 5),
  rate = c(25, 43, 66, 115, 25, 60, 8, 62)
)
ggplot(literature)+
  geom_col(aes(x = month , y = rate, fill = source), position = "dodge")

ggplot()+
  geom_line(aes(x = date_grid, y = NCP_mean, color = "Iceland Basin"), dat = icb)+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = "Iceland Basin"), dat = icb, alpha = 0.2)+
  geom_line(aes(x = date_grid, y = NCP_mean, color = "Irminger Sea"), dat = irm)+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = "Irminger Sea"), dat = irm, alpha = 0.2)+
  geom_line(aes(x = date_grid, y = NCP_mean, color = "Labrador Sea"), dat = lbr)+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = "Labrador Sea"), dat = lbr, alpha = 0.2)+
  geom_line(aes(x = date_grid, y = NCP_mean, color = "Norvegian Sea"), dat = nws)+
  geom_ribbon(aes(x = date_grid, ymin = NCP_q05, ymax = NCP_q95, fill = "Norvegian Sea"), dat = nws, alpha = 0.2)+
  scale_fill_brewer(palette= "Set1", guide = "none")+
  scale_color_brewer(palette = "Set1", name = "Basin")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")

icb_mthly <- icb %>% mutate(month = lubridate::month(date_grid)) %>% 
  group_by(month) %>% 
  summarise_all(.funs = list(mean, sd), na.rm = TRUE) %>% 
  ungroup()

irm_mthly <- irm %>% mutate(month = lubridate::month(date_grid)) %>% 
  group_by(month) %>% 
  summarise_all(.funs = list(mean, sd), na.rm = TRUE) %>% 
  ungroup()

lbr_mthly <- lbr %>% mutate(month = lubridate::month(date_grid)) %>% 
  group_by(month) %>% 
  summarise_all(.funs = list(mean, sd), na.rm = TRUE) %>% 
  ungroup()

nws_mthly <- nws %>% mutate(month = lubridate::month(date_grid)) %>% 
  group_by(month) %>% 
  summarise_all(.funs = list(mean, sd), na.rm = TRUE) %>% 
  ungroup()


region_combined <- bind_rows(icb, irm, lbr, nws) %>% 
  mutate(month = lubridate::month(date_grid),
         year = lubridate::year(date_grid),
         daily = NCP_mean *30,
         daily_sd = NCP_sd * 30) %>%
  group_by(region, year) %>% 
  summarise(ANCP = sum(daily)/1000,
            ANCP_sd = sum(daily_sd)/2000)

lbr_mthly <- lbr %>% mutate(month = lubridate::month(date_grid)) %>% 
  group_by(month) %>% 
  summarise_all(.funs = list(mean, sd), na.rm = TRUE) %>% 
  ungroup()
icb_ancp <- icb %>% mutate(month = lubridate::month(date_grid),
                           year = lubridate::year(date_grid),
                           daily = NCP_mean *15,
                           daily_sd = NCP_sd * 15) %>%
  group_by(year) %>% 
  summarise(ANCP = sum(daily)/1000,
            ANCP_sd = sum(daily_sd)/1000)
  
ggplot(region_combined)+
  geom_point(aes(x = year, y = ANCP, color = region))+
  geom_line(aes(x = year, y = ANCP, color = region))

ggplot(region_combined)+
  geom_col(aes(x = year, y = ANCP, fill = region), position = "dodge")+
  geom_errorbar(aes(x = year, ymin = ANCP - ANCP_sd, ymax = ANCP + ANCP_sd, group = region), position = "dodge")+
  ylab("ANCP (mol m-2 yr-1)")

tot <- region_combined %>% 
  group_by(region) %>% 
  summarise(ANCP_mean = mean(ANCP),
          ANCP_sum = sum(ANCP),
          sd = sd(ANCP)) %>% 
  ungroup()

ggplot(tot)+
  geom_col(aes(x = region, y = ANCP_mean))+
  geom_errorbar(aes(x = region, ymin = ANCP_mean - sd, ymax = ANCP_mean +sd))


ggplot()+
  geom_line(aes(x = month, y = NCP_mean_fn1, color = "icb"), dat = icb_mthly)+
  geom_ribbon(aes(x = month, ymin = NCP_mean_fn1 - NCP_mean_fn2, ymax = NCP_mean_fn1 + NCP_mean_fn2, fill = "icb"), data = icb_mthly, alpha = 0.2)+
  geom_line(aes(x = month, y = NCP_mean_fn1, color = "irm"), dat = irm_mthly)+
  geom_ribbon(aes(x = month, ymin = NCP_mean_fn1 - NCP_mean_fn2, ymax = NCP_mean_fn1 + NCP_mean_fn2, fill = "irm"), data = irm_mthly, alpha = 0.2)+
  geom_line(aes(x = month, y = NCP_mean_fn1, color = "lbr"), dat = lbr_mthly)+
  geom_ribbon(aes(x = month, ymin = NCP_mean_fn1 - NCP_mean_fn2, ymax = NCP_mean_fn1 + NCP_mean_fn2, fill = "lbr"), data = lbr_mthly, alpha = 0.2)+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  theme_bw()+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")

library(ggplot2)
library(maps)
library(dplyr)

# Define polygons
iceland_basin <- data.frame(
  lon = c(-32, -22, -12, -23),
  lat = c(57.4, 63, 63, 55.7),
  region = "Iceland Basin"
) 

irminger_sea <- data.frame(
  lon = c(-36.5, -26, -32.7, -42.8),
  lat = c(57, 63, 65.2, 59.29),
  region = "Irminger Sea"
)

labrador_sea <- data.frame(
  lon = c(-50.5, -44.5, -52.8, -60.4),
  lat = c(53.8, 58.5, 62.5, 60.4),
  region = "Labrador Sea"
)

northatl_box <- data.frame(
  lon = c(-12, -24, -32.68, -42.45, -36.5, -23, -16.9, -12),
  lat = c(63, 63, 65.2, 59.1, 57, 55.7, 59.7, 63),
  region = "North Atlantic"
)

norwegian_sea <- data.frame(
  lon = c(3, 8.5, -8.2, -11, 3),
  lat = c(62, 70, 70, 60.4, 62),
  region = "Norwegian Sea"
)


# Combine
regions <- bind_rows(iceland_basin, irminger_sea, labrador_sea, northatl_box, norwegian_sea)

# World map data
world <- map_data("world")

# Plot
ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey60", linewidth = 0.2) +
  
  geom_polygon(data = filter(regions, region != "North Atlantic"),
               aes(x = lon, y = lat, fill = region, group = region),
               color = "black", alpha = 0.4) +
  
  coord_fixed(xlim = c(-61, 10), ylim = c(50, 71), expand = FALSE) +
  
  scale_fill_manual(values = c("Iceland Basin" = "#e31a1c",
                               "Irminger Sea" = "#1f78b4",
                               "Labrador Sea" = "#4daf4a",
                               "North Atlantic" = "#984ea3",
                               "Norwegian Sea" = "#ff7f00")) +
  
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Region") +
  
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    legend.position = "right"
  )

ggsave("Output/Multiple_box_map.png", width = 16, height = 9, dpi = 300)

ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey60", linewidth = 0.2) +
  
  geom_polygon(data = filter(regions, region %in% c("Iceland Basin", "Irminger Sea")),
               aes(x = lon, y = lat, fill = region, group = region),
               color = "black", alpha = 0.4) +
  
  coord_fixed(xlim = c(-60, -10), ylim = c(50, 67), expand = FALSE)+
  
  scale_fill_manual(values = c("Iceland Basin" = "#e31a1c",
                               "Irminger Sea" = "#1f78b4",
                               "Labrador Sea" = "#4daf4a",
                               "North Atlantic" = "#984ea3",
                               "Norwegian Sea" = "#ff7f00")) +
  
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Region") +
  
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    legend.position = "right"
  )

ggsave("Output/Icb_Irm_map.png", width = 16, height = 9, dpi = 300)

ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey60", linewidth = 0.2) +
  
  geom_polygon(data = filter(regions, region == "North Atlantic"),
               aes(x = lon, y = lat, fill = region, group = region),
               color = "black", alpha = 0.4) +
  
  coord_fixed(xlim = c(-60, -10), ylim = c(50, 67), expand = FALSE) +
  
  scale_fill_manual(values = c("Iceland Basin" = "#e31a1c",
                               "Irminger Sea" = "#1f78b4",
                               "Labrador Sea" = "#4daf4a",
                               "North Atlantic" = "#984ea3")) +
  
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Region") +
  
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    legend.position = "right"
  )

ggsave("Output/NorthAtlantic1box_map.png", width = 16, height = 9, dpi = 300)

alex <- read_csv("C:/Users/petit/Downloads/plume_events_spatio_temp_extent.csv")

alex <- alex %>% group_by(WMO) %>% 
  select(WMO, longitude_max, longitude_min, latitude_min, latitude_max) %>% 
  summarise(lon_min = min(longitude_min),
             lon_max = max(longitude_max),
             lat_min = min(latitude_min),
             lat_max = max(latitude_max)) %>% 
  ungroup()

ggplot()+
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey60", linewidth = 0.2) +
  geom_rect(aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max, fill = as.factor(WMO)), data = alex)+
  coord_fixed(xlim = c(-60, 0), ylim = c(40, 75), expand = FALSE) +
  theme_minimal()
