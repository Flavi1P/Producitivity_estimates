library(tidyverse)
library(zoo)
library(gsw)
library(castr)
library(slider)
library(patchwork)
source("Scripts/cbpm_r/cbpm_argo.r")
source("Scripts/cbpm_r/despike.r")

files <- list.files("Data/Intermediate/Floats/", full.names = TRUE)

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
  geom_point(aes(x = date, y = - MLD, color = float_wmo))

dat_all_pp <- left_join(dat, prof_dat)

ggplot(select(dat_all_pp, float_wmo, lon, lat, date) |> distinct() |> filter(lon > -25 & lat  > 59))+
  geom_point(aes(x = lon, y = lat, color= date))+
  coord_quickmap()
  

data_to_share <- dat_all_pp |> 
  select(prof_number, lon, lat, date, oxygen, depth, temp, sal, MLD)

write_csv(data_to_share, "Data/Processed/float_oxygen_data_without_pp.csv")
