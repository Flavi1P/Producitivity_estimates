library(ncdf4)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)

# open the Sprof NetCDF
#nc <- nc_open("Data/Raw/Floats/3901586_Sprof.nc")

files <- list.files("Data/Raw/Floats/", pattern = "_Sprof.nc", full.names = TRUE)
for(file in files){
  print(paste(file, "started..."))
  nc <- nc_open(file)
  # extract metadata
  lon   <- ncvar_get(nc, "LONGITUDE")
  lat   <- ncvar_get(nc, "LATITUDE")
  juld  <- ncvar_get(nc, "JULD")  # days since 1950-01-01 00:00:00
  date  <- as.Date(juld, origin = "1950-01-01")
  
  # extract depth and variables (dimensions usually [levels x profiles])
  depth <- ncvar_get(nc, "PRES")  
  #ipar  <- ncvar_get(nc, "DOWNWELLING_PAR")
  temp  <- ncvar_get(nc, "TEMP")
  sal   <- ncvar_get(nc, "PSAL")
  oxygen <- ncvar_get(nc, "DOXY_ADJUSTED")
  
  qc <- ncvar_get(nc, "DOXY_ADJUSTED_QC")
  #create filename
  wmo <- ncvar_get(nc, "PLATFORM_NUMBER")[1]
  #remove whitespace
  wmo <- gsub(" ", "", wmo)
  filename <- paste0("Data/Intermediate/Floats/Doxy_db/argo_", wmo, "_interp.csv")
  
  # close connection
  nc_close(nc)
  
  # get dimensions
  n_levels   <- dim(depth)[1]
  n_profiles <- dim(depth)[2]
  
  # repeat metadata across depth levels
  lon_rep  <- rep(lon, each = n_levels)
  lat_rep  <- rep(lat, each = n_levels)
  date_rep <- rep(date, each = n_levels)
  
  # convert QC characters to numeric flags (1–9)
  doxy_qc <- as.numeric(charToRaw(paste(qc, collapse = "")))
  oxygen <- as.vector(oxygen)
  # keep only good QC (1 or 2)
  oxygen[!doxy_qc %in% c(49, 50)] <- NA
  
  # flatten arrays into vectors
  df <- tibble(
    lon   = lon_rep,
    lat   = lat_rep,
    date  = date_rep,
    depth = as.vector(depth),
    temp  = as.vector(temp),
    sal   = as.vector(sal),
    oxygen = oxygen
  )
  
  df_nona <- filter(df, !is.na(temp))
  df_nona <- df_nona |> mutate(depth = round(depth))
  
  df_interp <- df_nona %>%
    group_by(lon, lat, date) %>%
    arrange(depth, .by_group = TRUE) %>%
    reframe(
      depth = seq(0,
                  2000, 
                  by = 1)
    ) %>%
    left_join(df_nona, by = c("lon", "lat", "date", "depth")) %>%
    arrange(lon, lat, date, depth) %>%
    group_by(lon, lat, date) %>%
    mutate(across(
      c(oxygen, temp, sal),
      ~ if (sum(!is.na(.x)) > 10) {
        na.approx(.x, x = depth, na.rm = FALSE, rule = 2)
      } else {
        rep(NA_real_, length(.x))  # if no interpolation possible
      }
    )) %>%
    ungroup() |> 
    group_by(lon, lat, date, depth) |> 
    summarise_all(mean, na.rm = TRUE) |> 
    ungroup()
  
  df_interp <- df_interp %>%
    arrange(date) %>% 
    mutate(prof_number = dense_rank(date))
  
  t <- filter(df_interp, prof_number == 75)
  
  # save the csv ------------------------------------------------------------
  
  write_csv(df_interp, filename)
  print(paste0(filename, " written"))
}

ggplot(t)+
  geom_point(aes(x = oxygen, y = -depth))
