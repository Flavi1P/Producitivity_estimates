library(tidyverse)
library(zoo)
library(gsw)
library(castr)
library(slider)
library(patchwork)
source("Scripts/cbpm_r/cbpm_argo.r")
source("Scripts/cbpm_r/despike.r")

files <- list.files("Data/Intermediate/Full_doxy_db/", full.names = TRUE)

dat <- data.frame()
for(file in files){
  dat_temp <- read_csv(file)
  dat_temp$float_wmo = strsplit(basename(file), "_")[[1]][2]
  dat <- bind_rows(dat, dat_temp)
}

dat <- dat |> mutate(sa = gsw_SA_from_SP(sal, depth, lon ,lat), ct = gsw_CT_from_t(sa, temp, depth), sigma0 = gsw_sigma0(sa, ct))

range(dat$date)

dat = filter(dat, !is.na(lon))

# Prof dat ---------------------------------------------------------------

prof_dat <- dat |>
  group_by(prof_number, float_wmo, date) |>
  arrange(depth) |> 
  summarise(
    MLD = mld(sigma0, depth, ref.depths = 0:5, criteria = 0.03, default.depth = 300)
  ) |>
  ungroup() 


ggplot(prof_dat)+
  geom_point(aes(x = date, y = - MLD))

dat_all_pp <- left_join(dat, prof_dat)

ggplot(select(dat_all_pp, float_wmo, lon, lat, date) |> distinct())+
  geom_point(aes(x = lon, y = lat, color= date))+
  geom_rect(aes(xmin = -45, xmax = -10, ymin = 58, ymax = 65), fill = NA, color = "red")+
  coord_quickmap()
  

data_to_share <- dat_all_pp |> 
  select(float_wmo, prof_number, lon, lat, date, oxygen, depth, temp, sal, MLD) %>% 
  filter(lon >= -40 & lon <= -10 & lat >= 58 & lat <= 65)

write_csv(data_to_share, "Data/Processed/float_full_db_oxygen_data_without_pp.csv")
