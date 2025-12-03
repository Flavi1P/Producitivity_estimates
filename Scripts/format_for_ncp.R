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

# NPP data ---------------------------------------------------------------

library(vroom)

npp <- vroom("Data/Processed/argo_pp_estimations_floats_11_2025.csv")

npp <- npp |> mutate(date = lubridate::date(JULD),
                     float_wmo = as.character(PLATFORM_NUMBER)) |> 
  janitor::clean_names() |> 
  select(float_wmo, date, day, month, year, satellite_par) |> 
  unique()


dat_all_pp <- left_join(dat, npp, by = c("float_wmo", "date"))


dat_all_pp <- dat_all_pp |>
  filter(!is.na(satellite_par))

  
dat_all_pp <- 
  dat_all_pp |> 
  group_by(prof_number, float_wmo) |> 
  mutate(
    bbp_baseline = despike(bbp, 10)$baseline,
    bbp_baseline = case_when(is.na(bbp_baseline) ~ bbp,
                             TRUE ~ bbp_baseline),
    bbp_dark = mean(bbp_baseline[which(depth >= 900)], na.rm = TRUE),
    bb_baseline = bbp_baseline - bbp_dark,
    carbon = bbp_baseline / (470/440) * 12128 + 0.59,
    izeu_cbpm = cbpm_argo(chla, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), unique(na.omit(day)), unique(na.omit(lat)))$mzeu,
    cbpm_npp = cbpm_argo(chla, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), unique(na.omit(day)), unique(na.omit(lat)))$pp,
    cbpm_mu = cbpm_argo(chla, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), unique(na.omit(day)), unique(na.omit(lat)))$mu_z)


# Prof dat ---------------------------------------------------------------


zeu_calc <- function(par, depth){
  # find the depth where par first exceeds the threshold
  if(all(is.na(par)) | all(par < 1)){
    return(NA)
  } else {
    par_surf = par[which(depth == 1)]
    threshold = par_surf * 0.01
    return(depth[which(par <= threshold)[1]])
  }
}


prod_depth <- function(npp, depth){
  if(all(is.na(npp))){
    return(NA)
  } else {
    return(depth[which(npp < 0.01)[1]])
  }
}
prof_dat <- dat_all_pp |>
  group_by(prof_number, float_wmo, date) |>
  arrange(depth) |> 
  summarise(
    MLD = mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03, default.depth = 300),
    chl_cline = clined(chla, depth),
    zeu_cbpm = mean(izeu_cbpm, na.rm = TRUE),
    #zeu = zeu_calc(iPAR, depth),
    zeu = unique(zeu_cbpm),
    z_prod = prod_depth(cbpm_npp, depth)
  ) |>
  ungroup() 


ggplot(prof_dat)+
  geom_point(aes(x = date, y = - MLD, color = float_wmo))

dat_all_pp <- left_join(dat_all_pp, prof_dat)

ggplot(select(dat_all_pp, float_wmo, lon, lat, date) |> distinct())+
  geom_point(aes(x = lon, y = lat, color= date))+
  coord_quickmap()
  

data_to_share <- dat_all_pp |> 
  select(lon, lat, date, nitrate, oxygen, depth, chla, bbp, temp, sal, cbpm_npp, zeu_cbpm, zeu, MLD)

write_csv(data_to_share, "Data/Processed/float_nitrate_data_with_pp.csv")
