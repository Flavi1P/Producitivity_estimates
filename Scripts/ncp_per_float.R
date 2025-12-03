library(tidyverse)
library(lubridate)
library(castr)
library(zoo)

data <- read_csv("Data/Processed/float_nitrate_data_corrected.csv") |> 
  filter(lon >= -40 & lon <= -10 & lat >= 58 & lat <= 63) |> 
  filter(date >= as.Date("2019-01-01"))

smoothing = TRUE
smoothing_final = FALSE

regrid_time = TRUE

time_window = 15 #days

floats <- unique(data$float_wmo)

ncp_results_full <- tibble()

for(float in floats){
  # filter the dataset to stay with ICB -------------------------------------
  print(float)

  dat <- filter(data, float_wmo == float)
  #filter to keep data between 2019 and now
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
  
  if(smoothing == FALSE){
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
    
  }else{
    dat_prof <- dat_prof |> 
      mutate(mld_smooth = MLD,
             zeu_smooth = zeu)
  }
  
  
  dat <- left_join(dat, select(dat_prof, prof_number, float_wmo, mld_smooth, zeu_smooth), by = c("float_wmo", "prof_number"))
  
  
  # Integration depth computation -------------------------------------------
  
  
  #Define Integration depth (i.e. deepest value of MLD or Zeu between t and t-1)
  integration_depth_data <- dat |>
    group_by(float_wmo, prof_number, date) |>
    summarise(
      mld_smooth = unique(mld_smooth, na.rm = TRUE),
      zeu_smooth = unique(zeu_smooth, na.rm = TRUE)) |> 
    ungroup()
  
  time_space = paste(time_window, "days", sep = " ")
  
  #Regrid the data on a 10day interval timeseries
  time_grid <- tibble(date_10day = seq(min(integration_depth_data$date),
                                       max(integration_depth_data$date),
                                       by = time_space))
  if(regrid_time == TRUE){
    prof_dat_smoothed <- integration_depth_data %>%
      mutate(date_10day = as.Date(cut(date, breaks = time_space))) %>%
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
             NCP_integration_depth_smooth = case_when(smoothing == TRUE ~ rollapply(NCP_integration_depth, width = 3, FUN = mean, align = "center", fill = NA),
                                                      smoothing == FALSE ~ NCP_integration_depth),
             next_integration_depth_smooth = case_when(smoothing == TRUE ~ rollapply(next_integration_depth, width = 3, FUN = mean, align = "center", fill = NA),
                                                       smoothing == FALSE ~ next_integration_depth))
    
    #Create a dataset with the median value of nitrate integrals every 10 days for smoothing
    dat_smoothed <- dat %>%
      mutate(date_10day = as.Date(cut(date, breaks = time_space))) %>%
      group_by(date_10day, depth) %>%
      summarise(
        nitrate = mean(nitrate_corrected, na.rm = TRUE),
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
    
    dat_smoothed <- dat_smoothed |> 
      left_join(prof_dat_smoothed, by = "date_10day") |> 
      na.omit()
    
    
  }else{
    prof_dat_smoothed <- integration_depth_data %>%
      mutate(date_10day = date) %>%
      group_by(date_10day) %>%
      summarise(mld = mean(mld_smooth, na.rm = TRUE),
                zeu = mean(zeu_smooth, na.rm = TRUE)) |> 
      mutate(
        prev_MLD = dplyr::lag(mld),
        prev_zeu = dplyr::lag(zeu),
        NCP_integration_depth = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
        next_integration_depth = lead(NCP_integration_depth)
      ) |>
      arrange(date_10day) %>%
      mutate(
        NCP_integration_depth = zoo::na.approx(NCP_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
        next_integration_depth = zoo::na.approx(next_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)
      ) |> 
      mutate(dt = as.numeric(date_10day - dplyr::lag(date_10day)),
             NCP_integration_depth_smooth = case_when(smoothing == TRUE ~ rollapply(NCP_integration_depth, width = 3, FUN = mean, align = "center", fill = NA),
                                                      smoothing == FALSE ~ NCP_integration_depth),
             next_integration_depth_smooth = case_when(smoothing == TRUE ~ rollapply(next_integration_depth, width = 3, FUN = mean, align = "center", fill = NA),
                                                       smoothing == FALSE ~ next_integration_depth))
    
    #Create a dataset with the median value of nitrate integrals every 10 days for smoothing
    dat_smoothed <- dat %>%
      mutate(date_10day = date) %>%
      group_by(date_10day, depth) %>%
      summarise(
        nitrate = mean(nitrate_corrected, na.rm = TRUE),
        npp = mean(cbpm_npp, na.rm = TRUE),
        chla = mean(chla, na.rm = TRUE)) |> 
      ungroup() |>    # ensure regular 10-day spacing
      arrange(date_10day) %>%
      group_by(depth) |> 
      mutate(
        nitrate = na.approx(nitrate, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
        npp = na.approx(npp, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)) %>%
      ungroup()
    
    dat_smoothed <- dat_smoothed |> 
      left_join(prof_dat_smoothed, by = "date_10day") |> 
      na.omit()
    
  }

  prof_dat_smoothed <- na.omit(prof_dat_smoothed)
  
  
  # Interpolation of data  --------------------------------------------------
  

  
  
  # Integration of nitrate over integration depth ---------------------------
  
  dat_final <- dat_smoothed |>
    group_by(date_10day, mld, zeu) |>
    filter(!is.na(next_integration_depth)) |> 
    summarise(
      int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(NCP_integration_depth)),
      next_int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(next_integration_depth)),
      int_N_200 = integrate(nitrate, depth, from = 0, to = 200),
      npp = integrate(npp, depth, from = 0, to = 200),
      mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) - 20, to = unique(NCP_integration_depth) - 10) / 10,
      sub_mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) + 20, to = unique(NCP_integration_depth) + 30) / 10,
      chla_integrated = integrate(chla, depth, from = 0, to = unique(NCP_integration_depth)),
      chla_mld = integrate(chla, depth, from = 0, to = unique(mld)),
      chla_surf = integrate(chla, depth, from = 0, to = 10)
    ) |> 
    ungroup()

  #Smoothing the integration
  if(smoothing_final == TRUE){
    dat_final_smoothed <- dat_final %>%
      mutate(
        int_N_smooth = rollapply(int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA),
        next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA),
        int_N_200 = rollapply(int_N_200, width = 3, FUN = mean, align = "center", fill = NA)
      )
  }
  if(smoothing_final == FALSE){
    dat_final_smoothed <- dat_final |> 
      mutate(int_N_smooth = int_N_mmol_m2,
             next_int_N_smooth = next_int_N_mmol_m2,
             int_N_200 = int_N_200)
  }

  # Computing NCP -----------------------------------------------------------
  
  
  ncp_results <- dat_final_smoothed |> 
    arrange(date_10day) |> 
    mutate(dt = as.numeric(date_10day - dplyr::lag(date_10day)),
           diff_mld = mld - lag(mld),
           we = pmax(0, diff_mld),
           delta = sub_mld_concentration - mld_concentration,
           diff = delta * we,
           nitrate_consumption =  dplyr::lag(next_int_N_smooth) - int_N_smooth,
           nitrate_consumption200 = dplyr::lag(int_N_200) - int_N_200,
           c_consumption = nitrate_consumption * 6.625 / dt,
           c_consumption200 = nitrate_consumption200 * 6.625 / dt,
           NCP = c_consumption + diff,
           NCP200 = c_consumption200 + diff) |> 
    ungroup()
  
  ncp_results <- ncp_results |> 
    mutate(ncp_smooth = rollapply(NCP, width = 6, FUN = mean, align = "center", fill = NA),
           ncp_sd = rollapply(NCP, width = 6, FUN = sd, align = "center", fill = NA),
           npp_smooth = rollapply(npp, width = 6, FUN = median, align = "center", fill = NA),
           npp_sd = rollapply(npp, width = 6, FUN = sd, align = "center", fill = NA),
           ncp200_smooth = rollapply(NCP200, width = 6, FUN = mean, align = "center", fill = NA))
  
  ncp_results_full <- bind_rows(ncp_results_full,
                                    ncp_results |> mutate(float_wmo = float))
}

ggplot(filter(ncp_results_full, date_10day > as.Date("2024-01-01") & !float_wmo %in% c(3901581, 3901586)))+
  geom_line(aes(x =date_10day, y = NCP, color = as.factor(float_wmo), linetype = "NCP"))+
  #geom_point(aes(x = date_10day, y = NCP200, color = as.factor(float_wmo), shape = "NCP 200"))+
  geom_line(aes(x = date_10day, y = npp_smooth/12, color = as.factor(float_wmo), linetype = "NPP"))+
  theme_bw()+
  scale_color_brewer(palette = "Set1")+
  labs(x = "Date", y = "mmol C m-2 d-1", color = "Float WMO", linetype = "Variable")+
  theme(legend.position = "bottom")

ggsave("Output/ncp_npp_timeseries_per_float_filtered_smoothed_regrid.png", width = 10, height = 8, dpi = 300)

ncp_results_full <- ncp_results_full |> 
  mutate(npp = npp/12) |> 
  mutate(month_of_year = month(date_10day))

synthetic_ncp <- ncp_results_full |> 
  filter(float_wmo != 3901586 & float_wmo != 3901581) |> 
  mutate(month_of_year = month(date_10day)) |> 
  group_by(month_of_year) |> 
  summarise(mean_ncp = mean(NCP, na.rm = TRUE),
            sd_ncp = sd(NCP, na.rm = TRUE),
            mean_npp = mean(npp, na.rm = TRUE),
            sd_npp = sd(npp, na.rm = TRUE),
            mean_ncp200 = mean(NCP200, na.rm = TRUE),
            sd_ncp200 = sd(NCP200, na.rm = TRUE))

ggplot(synthetic_ncp)+
  geom_line(aes(x = month_of_year, y = mean_ncp, color = "NCP"))+
  geom_ribbon(aes(x = month_of_year, ymin = mean_ncp - sd_ncp, ymax = mean_ncp + sd_ncp), alpha = 0.2)+
  geom_line(aes(x = month_of_year, y = mean_npp, color = "NPP"))+
  geom_ribbon(aes(x = month_of_year, ymin = mean_npp - sd_npp, ymax = mean_npp + sd_npp), alpha = 0.2)+
  geom_line(aes(x = month_of_year, y = mean_ncp200, color = "NCP 200"))+
  geom_ribbon(aes(x = month_of_year, ymin = mean_ncp200 - sd_ncp200, ymax = mean_ncp200 + sd_ncp200), alpha = 0.2)+
  labs(x = "Month", y = "mmol C m-2 d-1", color = "Variable")+
  scale_x_continuous(breaks = 1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  theme_bw()+
  xlim(1,12.1)
ggsave("Output/synthetic_ncp_npp_seasonality_smoothed_regrid.png", width = 10, height = 8, dpi = 300)


ggplot(ncp_results_full)+
  geom_boxplot(aes(x = as.factor(month_of_year), y = NCP200, color = "NCP 200"))+
  geom_boxplot(aes(x = as.factor(month_of_year), y = NCP, color = "NCP"))+
  geom_boxplot(aes(x = as.factor(month_of_year), y = npp, color = "NPP"))+
  theme_bw()+
  scale_color_brewer(palette = "Set1")+
  xlab("Month")+
  scale_x_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  ylab("mmol C m-2 d-1")+
  ylim(c(-100, 200))

ggsave("Output/ncp_npp_boxplot_per_month_smoothed_regrid.png", width = 10, height = 8, dpi = 300)


