library(ncdf4)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)

file = "Data/Processed/ALR6_CTD_with_ecopuk.nc"
nc <- nc_open(file)
# extract metadata
profile_id <- ncvar_get(nc, "event_id")
profile_dir <- ncvar_get(nc, "event_direction")
lon   <- ncvar_get(nc, "LON")
lat   <- ncvar_get(nc, "LAT")
juld  <- ncvar_get(nc, "nav_timestamp")  # days since 1950-01-01 00:00:00
date  <- as.POSIXct(juld, origin = "1970-01-01", tz = "UTC")

# extract depth and variables (dimensions usually [levels x profiles])
depth <- ncvar_get(nc, "DEPTH_SCI_INTERP")  
chla  <- ncvar_get(nc, "fluo_corrected")
bbp   <- ncvar_get(nc, "bbp_baseline")   # sometimes "BBP700_ADJUSTED"
#ipar  <- ncvar_get(nc, "DOWNWELLING_PAR")
temp  <- ncvar_get(nc, "Temperature")
sal   <- ncvar_get(nc, "Conductivity")
fluo <- ncvar_get(nc, "fluo")
oxygen <- ncvar_get(nc, "Oxygen")

filename <- paste0("Data/Processed/ALR6_csv_from_nc.csv")

# close connection
nc_close(nc)

# flatten arrays into vectors
df <- tibble(
  profile_id = as.vector(profile_id),
  profile_dir = as.vector(profile_dir),
  lon   = as.vector(lon),
  lat   = as.vector(lat),
  date  = as.vector(date),
  depth = as.vector(depth),
  chla  = as.vector(chla),
  bbp   = as.vector(bbp),
  fluo = as.vector(fluo),
  temp  = as.vector(temp),
  sal   = as.vector(sal),
  oxygen = as.vector(oxygen)
)

df_nona <- filter(df, !is.na(temp))
df_nona <- df_nona |> mutate(depth = round(depth))

write_csv(df, filename)

# Plotting raw data -------------------------------------------------------

# # create folders
# dir.create("output/plots/temp_sal", recursive = TRUE, showWarnings = FALSE)
# dir.create("output/plots/chla_bbp", recursive = TRUE, showWarnings = FALSE)
# 
# for (p in unique(df_interp$prof_number)) {
#   
#   df_prof <- df_interp %>% filter(prof_number == p)
#   profile_date <- unique(df_prof$date)
#   
#   ## ---- Plot 1: Temperature + Salinity ----
#   scale_factor_ts <- max(df_prof$temp, na.rm = TRUE) / 
#     max(df_prof$sal,  na.rm = TRUE)
#   
#   p1 <- ggplot(df_prof, aes(y = depth)) +
#     geom_path(aes(x = temp, colour = "Temperature")) +
#     geom_path(aes(x = sal * scale_factor_ts, colour = "Salinity")) +
#     scale_y_reverse(expand = c(0,0)) +
#     scale_x_continuous(
#       name = "Temperature (°C)",
#       sec.axis = sec_axis(~ . / scale_factor_ts, name = "Salinity (PSU)")
#     ) +
#     labs(
#       title = paste("Profile", p, "-", profile_date),
#       y = "Depth (m)",
#       colour = ""
#     ) +
#     theme_minimal()
#   
#   ## ---- Plot 2: Chla + BBP (dual x axis, top 100 m) ----
#   df_prof100 <- df_prof %>% filter(depth <= 300)
#   
#   # scaling factor between Chla and BBP ranges
#   scale_factor <- max(df_prof100$chla, na.rm = TRUE) / 
#     max(df_prof100$bbp,  na.rm = TRUE)
#   
#   p2 <- ggplot(df_prof100, aes(y = depth)) +
#     geom_path(aes(x = chla, colour = "Chla")) +
#     geom_path(aes(x = bbp * scale_factor, colour = "BBP")) +
#     scale_y_reverse(expand = c(0,0)) +
#     scale_x_continuous(
#       name = "Chla (mg/m³)",
#       sec.axis = sec_axis(~ . / scale_factor, name = "BBP (1/m)")
#     ) +
#     labs(
#       title = paste("Profile", p, "-", profile_date, "(0–300 m)"),
#       y = "Depth (m)",
#       colour = ""
#     ) +
#     theme_minimal()
#   
#   ggsave(
#     filename = file.path("output/plots/chla_bbp", paste0("profile_", p, "_chla_bbp.png")),
#     plot = p2,
#     width = 6, height = 8
#   )
# }

