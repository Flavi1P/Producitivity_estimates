library(tidyverse)
library(lubridate)
library(castr)
library(zoo)

dat <- read_csv2("Data/Processed/productivity_float_with_canyon.csv")

# filter the dataset to stay with ICB -------------------------------------
#add iceland contour on map

ggplot(select(dat, lon, lat) |> distinct())+
  geom_point(aes(x = lon, y = lat))+
  geom_rect(aes(xmin = -45, xmax = -10, ymin = 58, ymax = 65), fill = NA, color = "red")+
  theme_minimal()+#add land contour
  borders("world", xlim = c(-60, 0), ylim = c(50, 70), fill = "lightgrey", color = "black")+
  coord_quickmap(xlim = c(-55, - 10), ylim = c(50, 65))

ggsave("Output/float_locations_icb.png", width = 8, height = 6)
dat <- filter(dat, lon >= -40 & lon <= -10 & lat >= 58 & lat <= 65) %>% 
  mutate(float_wmo = "prodctivity_float")
range(dat$date)

# MLD and Zeu time series -------------------------------------------------

#Artificially set zeu as 40 if NA
dat <- dat |> 
   mutate(zeu = 40)
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

fit_mld <- loess(
  MLD ~ as.numeric(date),
  data = dat_smooth,
  span = 0.3,
  family = "symmetric"   # <-- makes it robust to deep outliers
)

pred_mld <- predict(fit_mld, newdata = data.frame(date = as.numeric(dat_prof$date)))
pred_mld <- pmax(pred_mld, 0)   # physical constraint

# Predict smoothed values for all dates
#pred_mld <- predict(fit_mld, as.numeric(dat_prof$date))$y
pred_zeu <- predict(fit_zeu, as.numeric(dat_prof$date))$y


# Add back to main data frame
dat_prof <- dat_prof %>%
  mutate(mld_smooth = pred_mld,
         zeu_smooth = pred_zeu)


ggplot(dat_prof)+
  geom_line(aes(x = date, y = -mld_smooth, color = "MLD"))+
  geom_point(aes(x = date, y = -MLD), shape = 1)+
  geom_line(aes(x = date, y = -zeu_smooth, color = "Zeu"))+
  labs(y = "Depth (m)", color = "Parameter")+
  theme_minimal()
#ggsave("Output/mld_zeu_time_series_icb2.png", width = 10, height = 6)
  
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
                                     by = "5 days"))

prof_dat_smoothed <- integration_depth_data %>%
  mutate(date_10day = as.Date(cut(date, breaks = "5 days"))) %>%
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
  ) %>%
  filter(date_10day %in% time_grid$date_10day) |> 
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
  mutate(date_10day = as.Date(cut(date, breaks = "5 days"))) %>%
  group_by(date_10day, depth) %>%
  summarise(
    nitrate = mean(canyon_nitrate, na.rm = TRUE)) |> 
  ungroup() |> 
  full_join(time_grid, by = "date_10day") %>%     # ensure regular 10-day spacing
  arrange(date_10day) %>%
  group_by(depth) |> 
  mutate(
    nitrate = na.approx(nitrate, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)) %>%
  ungroup() |> 
  filter(date_10day %in% time_grid$date_10day)

#ggplot(dat_smoothed)+
#  geom_point(aes(x = nitrate, y = - depth, color = date_10day))
#ggsave("Output/nitrate_profiles_icb.png", width = 8, height = 6)
dat_smoothed <- dat_smoothed |> 
  left_join(prof_dat_smoothed, by = "date_10day") |> 
  na.omit()

ggplot(dat_smoothed)+
  geom_raster(aes(x = date_10day, y = -depth, fill = nitrate))+
  geom_line(aes(x = date_10day, y = - mld), color = "white")+
  scale_fill_viridis_c()+
  ylim(c(-500, 0))
# Integration of nitrate over integration depth ---------------------------

dat_final <- dat_smoothed |>
  group_by(date_10day, mld, zeu) |>
  filter(!is.na(next_integration_depth)) |> 
  summarise(
    int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(NCP_integration_depth)),
    next_int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(next_integration_depth)),
    mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) - 20, to = unique(NCP_integration_depth) - 10) / 10,
    sub_mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) + 20, to = unique(NCP_integration_depth) + 30) / 10
  ) |> 
  ungroup()

ggplot(dat_final)+
  geom_line(aes(x = date_10day, y = int_N_mmol_m2))+
  geom_point(aes(x = date_10day, y = int_N_mmol_m2))+
  geom_line(aes(x = date_10day, y = next_int_N_mmol_m2), linetype = "dashed")+
  labs(y = "Integrated Nitrate (mmol N m-2)")+
  theme_minimal()+
  scale_x_date(date_labels = "%b %Y", breaks = "1 months")

ggsave("Output/integrated_nitrate_time_series_icb2.png", width = 10, height = 6)

#ggsave("Output/integrated_chla_time_series_icb2.png", width = 10, height = 6)


#Smoothing the integration
dat_final_smoothed <- dat_final %>%
  mutate(
    int_N_smooth = rollapply(int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA),
    next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA)
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
         nitrate_consumption =  dplyr::lag(next_int_N_smooth) - int_N_smooth,
         c_consumption = nitrate_consumption * 6.625 / dt,
         NCP = c_consumption + diff) |> 
  ungroup()

ncp_results <- ncp_results |> 
  mutate(ncp_smooth = rollapply(NCP, width = 6, FUN = mean, align = "center", fill = NA),
         ncp_sd = rollapply(NCP, width = 6, FUN = sd, align = "center", fill = NA),
         ncp_signif = rollapply(ncp_smooth, width = 6, FUN = mean, align = "center", fill = NA)) %>% 
  mutate(integration = case_when(zeu > mld ~ "zeu",
                                 TRUE ~ "mld"))

peak_ncp_date = ncp_results %>% filter(ncp_smooth == max(ncp_smooth, na.rm = TRUE))

ncp_results = filter(ncp_results, !is.na(NCP)) %>% 
  mutate(datenum = as.numeric(date_10day))

ncp_results <- ncp_results %>% 
  mutate(NCP= case_when(NCP < -60 ~ NCP /2,
                        TRUE ~ NCP))

fit <- loess(NCP ~datenum, data = ncp_results, span = 0.4, family = "symmetric")

ncp_results$loess_ncp = predict(fit)

ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = NCP, color = "Raw 5 days NCP"), alpha = 0.8)+
  geom_line(aes(x = date_10day, y = ncp_smooth, color = "30 days avg"), linewidth = 1.5)+
  geom_ribbon(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd), alpha = 0.2)+
  geom_line(aes( x = date_10day, y = loess_ncp, color = "Loess"), linewidth = 1.5)+
  scale_color_brewer(palette = "Set1")+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  xlab("Date")+
  ggtitle("NCP estimates from Nitrate drawdown\npredicted from CANYON-B applied on productivity float DOXY observation")



# plotting Nathan outputs -------------------------------------------------

o2_net_change <- read_csv("Data/Processed/O2_float_net_change.csv")%>% 
  mutate(date = as.POSIXct((mtime - 719529) * 86400, origin = "1970-01-01", tz = "UTC"),
         datenum = as.numeric(date),
         o2_net_smooth = rollapply(Int_dO2dt_40m_mmol_m2_d/2, width = 6, FUN = mean, align = "center", fill = NA)) %>% 
  na.omit()


fit <- loess(o2_net_smooth~datenum, data = o2_net_change, span = 0.2, family = "symmetric")

o2_net_change$o2_net_change = predict(fit)

o2_float_night_loss <- read_csv("Data/Processed/O2_float_night_loss.csv")%>% 
  mutate(date = as.POSIXct((mtime - 719529) * 86400, origin = "1970-01-01", tz = "UTC"),
         datenum = as.numeric(date),
         o2_night_smooth = rollapply(Int_O2_loss_40m_mmol_m2_d/2, width = 10, FUN = mean, align = "center", fill = NA)) %>% 
  na.omit()


fit <- loess(o2_night_smooth~datenum, data = o2_float_night_loss, span = 0.2, family = "symmetric")

o2_float_night_loss$o2_night_change = predict(fit)


npp <- read_csv("Data/Processed/icb_npp_ncp.csv") %>% 
  mutate(npp = npp/12,
         datenum = as.numeric(date_10day)) %>% 
  na.omit()

ggplot(ncp_results) +
  # NCP pair
  geom_line(aes(x = date_10day, y = NCP,
                color = "NCP light"), alpha = 0.8) +
  geom_line(aes(x = date_10day, y = loess_ncp,
                color = "NCP dark"), linewidth = 1) +
  
  # O2 NIGHT CHANGE pair
  geom_line(aes(x = date, y = o2_night_smooth,
                color = "O2 night light"),
            data = o2_float_night_loss, alpha = 0.8) +
  geom_line(aes(x = date, y = o2_night_change,
                color = "O2 night dark"),
            data = o2_float_night_loss, linewidth = 1) +
  
  # O2 NET CHANGE pair
  geom_line(aes(x = date, y = o2_net_smooth,
                color = "O2 net light"),
            data = o2_net_change, alpha = 0.8) +
  geom_line(aes(x = date, y = o2_net_change,
                color = "O2 net dark"),
            data = o2_net_change, linewidth = 1) +
  
  #NPP
  geom_line(aes(x = date_10day, y = npp, color = "NPP"), data = npp)+
  
  
  # ---- COLOR MANUAL WITH MATCHING PAIRS ----
scale_color_manual(
  name = "",
  values = c(
    # NCP
    "NCP light"     = "#6baed6",
    "NCP dark"      = "#08519c",
    "NCP dark2" = "#810f7c",
    
    # O2 NIGHT CHANGE
    "O2 night light" = "#fc9272",
    "O2 night dark"  = "#cb181d",
    
    # O2 NET CHANGE
    "O2 net light"   = "#31a354",
    "O2 net dark"    = "#006837",
    "NPP" = "#252525"
  ),
  labels = c(
    "NCP light"        = "NCP (30 days avg)",
    "NCP dark"         = "NCP (loess)",
    "NCP dark2" = "NCP (other floats)",
    "O2 night light"   = "O2 night change (30 days avg)",
    "O2 night dark"    = "O2 night change (loess)",
    "O2 net light"     = "O2 net change (18 days avg)",
    "O2 net dark"      = "O2 net change (loess)",
    "NPP" = "NPP (other floats)"
    )
  ) +
  labs(y = "mmol C m-2 d-1") +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
  xlab("Date") +
  ggtitle("NCP estimation compared to net and night O2 changes")

ggsave("Output/productivity_time_series_total.png", width = 10, height = 6)


# NAB08 data --------------------------------------------------------------

nab <- read_csv("Data/Processed/NAB08_PP.csv")%>% 
  mutate(date = as.POSIXct((mtime - 719529) * 86400, origin = "1987-01-01", tz = "UTC"))


ggplot(ncp_results) +
  # NCP pair
  geom_line(aes(x = date_10day, y = NCP,
                color = "NCP light"), alpha = 0.5, linetype = "dashed") +
  geom_line(aes(x = date_10day, y = loess_ncp,
                color = "NCP dark"), alpha = 0.5, linetype = "dashed") +
  
  geom_line(aes(x = date, y = NPP_mmol_m2, color = "NPP"), data = nab)+
  
  # O2 NIGHT CHANGE pair
  geom_line(aes(x = date, y = o2_smooth,
                color = "O2 night light"),
            data = o2_float_night_loss, alpha = 0.5, linetype = "dashed") +
  geom_line(aes(x = date, y = o2_night_change,
                color = "O2 night dark"),
            data = o2_float_night_loss, alpha = 0.5, linetype = "dashed") +
  
  # O2 NET CHANGE pair
  geom_line(aes(x = date, y = o2_smooth,
                color = "O2 net light"),
            data = o2_net_change, alpha = 0.5, linetype = "dashed") +
  geom_line(aes(x = date, y = o2_net_change,
                color = "O2 net dark"),
            data = o2_net_change, alpha = 0.5, linetype = "dashed") +
  geom_line(aes(x = date, y = GPP_mmol_m2, color = "NCP dark2"), data = nab)+
  
  #NPP
  geom_line(aes(x = date_10day, y = npp, color = "NPP"), data = npp, linetype = "dashed", alpha = 0.5)+
  
  
  # ---- COLOR MANUAL WITH MATCHING PAIRS ----
scale_color_manual(
  name = "",
  values = c(
    # NCP
    "NCP light"     = "#6baed6",
    "NCP dark"      = "#08519c",
    "NCP dark2" = "#810f7c",
    
    # O2 NIGHT CHANGE
    "O2 night light" = "#fc9272",
    "O2 night dark"  = "#cb181d",
    
    # O2 NET CHANGE
    "O2 net light"   = "#31a354",
    "O2 net dark"    = "#006837",
    "NPP" = "#252525"
  ),
  labels = c(
    "NCP light"        = "NCP (30 days avg)",
    "NCP dark"         = "NCP (loess)",
    "NCP dark2" = "GPP NAB08",
    "O2 night light"   = "O2 night change (30 days avg)",
    "O2 night dark"    = "O2 night change (loess)",
    "O2 net light"     = "O2 net change (18 days avg)",
    "O2 net dark"      = "O2 net change (loess)",
    "NPP" = "NPP NAB08"
  )
) +
  labs(y = "mmol C m-2 d-1") +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
  xlab("Date") +
  ggtitle("NCP estimation compared to net and night O2 changes")+
  xlim(as.Date("2025-03-01"), as.Date("2025-07-01"))

  ggsave("Output/productivity_time_series_compared_nab08_zoomed.png", width = 10, height = 6)
  
ggplot(nab)+
  geom_point(aes(x = NPP_mmol_m2, y = GPP_mmol_m2))+
  geom_abline(intercept = 0, slope = 1)

o2_tot <- full_join(o2_float_night_loss, o2_net_change)


#na pparox o2_net_change and o2_night_smooth
o2_tot <- o2_tot %>% 
  arrange(datenum) %>% 
  select(date, datenum, o2_net_smooth, o2_night_smooth) %>% 
  mutate(o2_net_smooth = na.approx(o2_net_smooth, date, na.rm = FALSE),
         o2_night_smooth = na.approx(o2_night_smooth, date, na.rm = FALSE)) %>% 
  na.omit() %>% 
  mutate(GPP = o2_net_smooth + o2_night_smooth)
         

ggplot()+
  geom_line(aes(x = date, y = GPP, color = "GPP"), data = o2_tot)+
  geom_line(aes(x = date_10day, y = npp, color = "NPP"), data = npp)+
  labs(y = "mmol C m-2 d-1") +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
  xlab("Date") 

library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare data: keep distinct point locations and extract month
dat2 <- dat %>%
  select(date, lon, lat) %>%
  distinct() %>%
  mutate(month = month(date))   # abbreviated month label

# Calculate bounding box around your points, with a small margin
lon_range <- range(dat2$lon, na.rm = TRUE)
lat_range <- range(dat2$lat, na.rm = TRUE)

lon_margin <- diff(lon_range) * 0.1
lat_margin <- diff(lat_range) * 0.1

ggplot() +
  geom_sf(data = world, fill = "gray95", color = "black") +
  geom_point(
    data = dat2,
    aes(x = lon, y = lat, colour = month),
    size = 2
  ) +
  scale_color_viridis_c()+
  coord_sf(
    xlim = c(lon_range[1] - lon_margin, lon_range[2] + lon_margin),
    ylim = c(lat_range[1] - lat_margin, lat_range[2] + lat_margin),
    expand = FALSE
  ) +
  theme_minimal() +
  labs(colour = "Month")

ggsave("Output/productivity_floaat.png", width = 10, height = 6)

# statistical analysis ----------------------------------------------------

ggplot(ncp_results)+ 
  geom_line(aes( x = date_10day, y = loess_ncp, color = "Loess"), linewidth = 1.5)+
  scale_color_brewer(palette = "Set1")+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  xlab("Date")+
  ggtitle("NCP estimates from Nitrate drawdown\npredicted from CANYON-B applied on productivity float DOXY observation")

library(forecast)

ts_data <- ts(ncp_results$loess_ncp, frequency = 73)
fit <- mstl(ts_data)
autoplot(fit)

ncp_results$seasonal <- fit[, 2]

model = lm(seasonal~date_10day, data = ncp_results)
summary(model)

ggplot(ncp_results)+ 
  geom_line(aes( x = date_10day, y = loess_ncp, color = "Loess"), linewidth = 1.5)+
  geom_abline(intercept = 88.13, slope = -0.0046)+
  scale_color_brewer(palette = "Set1")+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  xlab("Date")+
  ggtitle("NCP estimates from Nitrate drawdown\npredicted from CANYON-B applied on productivity float DOXY observation")


ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = NCP, color = "Raw 5 days NCP"), alpha = 0.8)+
  geom_line(aes(x = date_10day, y = ncp_smooth, color = "30 days avg"), linewidth = 1.5)+
  geom_ribbon(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd), alpha = 0.2)+
  geom_line(aes( x = date_10day, y = seasonal, color = "Loess"), linewidth = 1.5)+
  scale_color_brewer(palette = "Set1")+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "6 months")+
  xlab("Date")+
  ylim(c(-100, 100))+
  ggtitle("NCP estimates from Nitrate drawdown\npredicted from CANYON-B applied on productivity float DOXY observation")

synthetic_ncp <- ncp_results |> 
  mutate(month_of_year = yday(date_10day)) |> 
  group_by(month_of_year) |> 
  summarise(mean_ncp = mean(NCP, na.rm = TRUE),
            sd_ncp = sd(NCP, na.rm = TRUE))

ggplot(synthetic_ncp)+
  geom_line(aes(x = month_of_year, y = mean_ncp, color = "NCP"))+
  geom_ribbon(aes(x = month_of_year, ymin = mean_ncp - sd_ncp, ymax = mean_ncp + sd_ncp), alpha = 0.2)+
  labs(x = "Month", y = "mmol C m-2 d-1", color = "Variable")+
  scale_x_continuous(breaks = 1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  theme_bw()+
  xlim(1,12.1)


ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = NCP, color = "Raw 5 days NCP"), alpha = 0.8)+
  #geom_line(aes(x = date_10day, y = ncp_smooth, color = "30 days avg"), linewidth = 1.5)+
  #geom_ribbon(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd), alpha = 0.2)+
  geom_line(aes( x = date_10day, y = loess_ncp, color = "Loess"), linewidth = 1.5)+
  geom_curve(aes(x = peak_ncp_date$date_10day + 50, y = peak_ncp_date$ncp_smooth + 50, xend = peak_ncp_date$date_10day, yend = peak_ncp_date$ncp_smooth),
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black",
             size = 1.2,
             curvature = 0.3
  )+
  geom_text(aes(x = peak_ncp_date$date_10day + 50, y = peak_ncp_date$ncp_smooth + 70, label = "Max NCP = 47 mmol c m-2, 5th May"))+
  scale_color_brewer(palette = "Set1")+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  xlab("Date")+
  ggtitle("NCP estimates from Nitrate drawdown\npredicted from CANYON-B applied on productivity float DOXY observation")

ggsave("Output/ncp_time_series_prod_float.png", width = 10, height = 6)
ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = -mld))+
  geom_vline(aes(xintercept = peak_ncp_date$date_10day))

ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = NCP, color = "Raw 5 days NCP"), alpha = 0.8)+
  geom_line(aes(x = date_10day, y = ncp_smooth, color = "30 days avg"), linewidth = 1.5)+
  geom_ribbon(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd), alpha = 0.2)+
  geom_line(aes( x = date_10day, y = ncp_signif, color = "double avg"), linewidth = 1.5)+
  geom_curve(aes(x = peak_ncp_date$date_10day + 50, y = peak_ncp_date$ncp_smooth + 50, xend = peak_ncp_date$date_10day, yend = peak_ncp_date$ncp_smooth),
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "black",
             size = 1.2,
             curvature = 0.3
  )+
  geom_text(aes(x = peak_ncp_date$date_10day + 50, y = peak_ncp_date$ncp_smooth + 70, label = "Max NCP = 47 mmol c m-2, 5th May"))+
  scale_color_brewer(palette = "Set1")+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  xlab("Date")+
  ylim(c(-100, 100))+
  ggtitle("NCP estimates from Nitrate drawdown\npredicted from CANYON-B applied on productivity float DOXY observation\nzoomed on Y axis to visually remove outliers")


ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = ncp_smooth, color = integration, group = float_wmo), linewidth = 1.5)+
  geom_line(aes(x = date_10day, y = ncp_signif, color = "cleaned", group = float_wmo), linewidth = 1.5)+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  xlab("Date")

ggplot(filter(ncp_results, npp != 0))+
  geom_line(aes(x = date_10day, y = npp_smooth/12, color = "NPP"))+
  geom_errorbar(aes(x = date_10day, y = npp_smooth/12, ymin = npp_smooth/12 - npp_sd/12, ymax = npp_smooth/12 + npp_sd/12), alpha = 0.3)+
  geom_errorbar(aes(x = date_10day, ymin = ncp_smooth - ncp_sd, ymax = ncp_smooth + ncp_sd), alpha = 0.3)+
  geom_line(aes(x = date_10day, y = ncp_smooth, color = "NCP"))+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()+
  scale_x_date(date_labels = "%b %Y", breaks = "2 months")+
  xlab("Date")

#ggsave("Output/ncp_time_series_icb3.png", width = 10, height = 6)

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
