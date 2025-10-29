library(tidyverse)
library(patchwork)
library(zoo)
#Read the data
dat <- read_csv2("Data/Processed/float_icb_ts_with_canyon.csv")

#Ensure each profile is reaching at least 1000m depth

dat_filtered <- dat %>%
  group_by(float_wmo, prof_number) %>%
  filter(max(depth, na.rm = TRUE) >= 1000) %>%
  ungroup()

#Interpolate NA values of nitrate within each profile
dat_interpolated <- dat_filtered %>%  group_by(float_wmo, prof_number) %>%
  arrange(depth) %>%
  mutate(nitrate_interp = zoo::na.approx(nitrate, depth, na.rm = FALSE, rule = 2)) %>%
  ungroup()

#remove NAs after interpolation and remove profile with less than 1000 observation of nitrate_intepolated
dat_interpolated <- dat_interpolated %>%
  group_by(float_wmo, prof_number) %>%
  filter(sum(!is.na(nitrate_interp)) >= 1000) %>%
  ungroup()
#Calculate the nitrate correction offset for each profile as being the average value of the last ten meters of each profile
nitrate_offsets <- dat_interpolated %>%
  group_by(float_wmo, prof_number, date) %>%
  filter(depth >= 1100) %>%
  summarise(nitrate_offset = median(nitrate_interp, na.rm = TRUE),
            nitrate_ref = median(canyon_nitrate, na.rm = TRUE)) %>%
  ungroup()


#Plot the offset timeseries of each float
ggplot(nitrate_offsets) +
  geom_line(aes(x = date, y = nitrate_offset, color = as.factor(float_wmo))) +
  geom_line(aes(x = date, y = nitrate_ref, color = as.factor(float_wmo)), linetype = "dashed") +
  labs(x = "Date", y = "Nitrate Offset (µmol/kg)", color = "Float WMO") +
  theme_minimal()


#Apply substract the difference between nitrate_offset and nitrate_ref to the nitrate_interp values
dat_corrected <- dat_interpolated %>%
  left_join(nitrate_offsets, by = c("float_wmo", "prof_number", "date")) %>%
  mutate(nitrate_corrected = nitrate_interp - (nitrate_offset - nitrate_ref))

ggplot(dat_corrected)+
  geom_boxplot(aes(x = as.factor(float_wmo), y = nitrate_interp))+
  labs(x = "Float WMO", y = "Nitrate (µmol/kg)")+
  theme_minimal()+
  ggtitle("Before Correction")+

ggplot(dat_corrected)+
  geom_boxplot(aes(x = as.factor(float_wmo), y = nitrate_corrected))+
  labs(x = "Float WMO", y = "Nitrate (µmol/kg)")+
  theme_minimal()+
  ggtitle("After Correction")

ggsave("Output/nitrate_correction_boxplots.png", width = 29, height = 21, units = "cm", dpi = 300)
#Making sure each profile has one value of nitrate_corrected every meter from 0 to 1000m
dat_final <- dat_corrected %>%
  group_by(float_wmo, prof_number) %>%
  arrange(depth) %>%
  reframe(
    depth = seq(0, 1000, by = 1)
  ) %>%
  left_join(dat_corrected, by = c("float_wmo", "prof_number", "depth")) %>%
  arrange(float_wmo, prof_number, depth) %>%
  group_by(float_wmo, prof_number) %>%
  mutate(across(
    c(nitrate_corrected),
    ~ if (sum(!is.na(.x)) > 1) {
        na.approx(.x, x = depth, na.rm = FALSE, rule = 2)
      } else {
        rep(NA_real_, length(.x))  # if no interpolation possible   
      }
  )) %>%
  ungroup()

dat_final <- dat_final |> filter(!is.na(nitrate_corrected))

#keep only profiles with 1000 rows
dat_final <- dat_final %>%
  group_by(float_wmo, prof_number) %>%
  filter(n() == 1001) %>%
  ungroup()
#Save the final data
write_csv(dat_final, "Data/Processed/float_nitrate_data_corrected.csv")
