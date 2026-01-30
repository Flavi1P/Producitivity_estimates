library(tidyverse)
library(gsw)
dat <- read_csv("Data/Processed/ALR6_csv_from_nc.csv") %>% 
  filter(profile_dir %in% c("ascent", "descent"))

pos <- dat %>% select(lon, lat, chla) %>% distinct() %>% filter(!is.na(chla))

dat$date = as.POSIXct(dat$date, origin = "1970-01-01", tz = "UTC")

dat$SP <- gsw_SP_from_C(
  C = dat$sal,      # conductivity in mS/cm → must convert to S/m
  t = dat$temp,
  p = dat$depth   # dbar (0 if unknown)
)

# Compute IQR limits
Q1  <- quantile(dat$SP, 0.25, na.rm = TRUE)
Q3  <- quantile(dat$SP, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

lower <- Q1 - 1.5 * IQR
upper <- Q3 + 1.5 * IQR

# Filter dataset
df_clean <- dat[dat$SP >= lower & dat$SP <= upper, ]

df_clean_reduced <- df_clean %>% 
  mutate(depth_round = round(depth)) %>% 
  group_by(profile_id, depth_round) %>% 
  summarise_all(mean, na.rm = TRUE)
# ggplot(df_clean_reduced, aes(x = SP, y = temp, color = lat)) +
#   geom_point(alpha = 0.2, size = 2) +
#   theme_minimal(base_size = 14) +
#   labs(
#     title = "TS Diagram",
#     x = "Prac. Salinity",
#     y = "Temperature (°C)"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )+
#   scale_color_viridis_c(name = "Latitude")
# 
# ggsave("output/TS_diagram_lat.png", width = 10, height = 8, dpi = 300)
# 
# ggplot(df_clean_reduced, aes(x = SP, y = temp, color = depth_round)) +
#   geom_point(alpha = 0.2, size = 2) +
#   theme_minimal(base_size = 14) +
#   labs(
#     title = "TS Diagram",
#     x = "Prac. Salinity",
#     y = "Temperature (°C)"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )+
#   scale_color_viridis_c(name = "Depth")
# 
# ggsave("output/TS_diagram_depth.png", width = 10, height = 8, dpi = 300)
# 
# ggplot(df_clean_reduced %>% filter(depth_round < 100), aes(x = SP, y = temp, color = lat)) +
#   geom_point(alpha = 0.5, size = 2) +
#   theme_minimal(base_size = 14) +
#   labs(
#     title = "TS Diagram (first 100m)",
#     x = "Prac. Salinity",
#     y = "Temperature (°C)"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )+
#   scale_color_viridis_c(name = "Latitude")
# 
# ggsave("output/TS_diagram_lat_surf.png", width = 10, height = 8, dpi = 300)
# 
# 
ggplot(df_clean_reduced %>% filter(depth_round < 100), aes(x = SP, y = temp, color = depth_round)) +
  geom_point(alpha = 0.5, size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "TS Diagram (first 100m)",
    x = "Prac. Salinity",
    y = "Temperature (°C)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+
  scale_color_viridis_c(name = "Depth")

ggsave("output/TS_diagram_lat_surf_depth.png", width = 18, height = 10, dpi = 300)



# clustering --------------------------------------------------------------

# ---- 1. Select variables for clustering ----
df_ts <- df_clean_reduced %>%
  ungroup() %>% 
  filter(depth_round <= 100) %>%
  select(temp, SP, lat, depth)

# Remove NAs
df_ts <- na.omit(df_ts)

# ---- 2. Scale data (important for T–S clustering) ----
ts_scaled <- scale(df_ts)

# ---- 3. Choose number of clusters (water masses) ----
k <- 4   # adjust to your dataset

set.seed(123)
km <- kmeans(ts_scaled, centers = k)

df_ts$cluster <- factor(km$cluster)

# ---- 4. Plot clustered TS diagram ----
ggplot(df_ts, aes(x = SP, y = temp, color = cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 16) +
  labs(
    title = "TS Diagram with Clustered Water Masses",
    x = "Practical Salinity",
    y = "Temperature (°C)",
    color = "Water Mass"
  )
ggsave("output/TS_clustered.png", width = 18, height = 10, dpi = 300)



df_to_plot = left_join(df_clean_reduced, df_ts) %>% 
  filter(depth_round < 100)

ggplot(df_to_plot)+
  geom_point(aes(x = date, y = -depth_round, color = cluster))+
  scale_color_manual(
    values = c(
      "1" = "#e41a1c",
      "2" = "#4daf4a",
      "3" = "#377eb8",
      "4" = "#984ea3"
    )
  ) +
  theme_minimal(base_size = 16) +
  ylab("Depth")

ggsave("output/cluster_transect.png", width = 18, height = 10, dpi = 300)


ggplot(df_to_plot)+
  geom_point(aes(x = date, y = -depth_round, color = bbp))+
  scale_color_viridis_c()

ggsave("output/cluster_bbp.png", width = 10, height = 8, dpi = 300)



ggplot(df_to_plot)+
  geom_point(aes(x = date, y = -depth_round, color = fluo))+
  scale_color_viridis_c()



ggplot(df_to_plot)+
  geom_boxplot(aes(x = cluster, y = chla/bbp, fill = cluster))+
  ylim(0, 2000)


df_corr <- df_clean_reduced %>%
  arrange(profile_id, depth_round) %>%       # ensure proper order
  group_by(profile_id) %>%             # do this per profile
  mutate(
    # ---- 1. Compute MLD (0.2°C threshold from surface) ----
    T0 = first(temp),
    mld = {
      idx <- which(abs(temp - T0) > 0.2)[1]
      if (is.na(idx)) max(depth_round) else depth_round[idx]
    },
    
    # ---- 2. Find max fluorescence *below* the MLD (where NPQ is minimal) ----
    fluo_sub_mld_max = max(fluo[depth_round > mld], na.rm = TRUE),
    
    # ---- 3. Detect NPQ: fluo decreasing upward inside MLD ----
    # inside MLD, if fluo < deep-max → NPQ happening
    npq_flag = ifelse(depth_round <= mld & fluo < fluo_sub_mld_max, TRUE, FALSE),
    
    # ---- 4. Apply Xin 2012 NPQ correction ----
    FLUO_NPQ = ifelse(npq_flag, fluo_sub_mld_max, fluo)
  ) %>%
  ungroup()

ggplot(df_corr %>% filter(depth_round < 50) %>% arrange(date))+
  geom_path(aes(x = date, y = -depth_round, color = FLUO_NPQ), linewidth = 4)+
  geom_line(aes(x = date, y = -mld), color = "white")+
  scale_color_viridis_c(name = "FChla")+
  ylim(-50, 0)+
  theme_minimal(base_size = 16)+
  ylab("Depth (m)")

ggsave("output/transect_chla_zoomed.png", width = 18, height = 10, dpi = 300)

ggplot(df_corr %>% filter(depth_round < 50 & bbp < 0.006) %>% arrange(date))+
  geom_path(aes(x = date, y = -depth_round, color = bbp), linewidth = 4)+
  geom_line(aes(x = date, y = -mld), color = "white")+
  scale_color_viridis_c(name = "bbp")+
  ylim(-50, 0)+
  theme_minimal(base_size = 16)+
  ylab("Depth (m)")
ggsave("output/transect_bbp_zoomed.png", width = 18, height = 10, dpi = 300)

df_final <- df_corr %>% 
  filter(depth_round < 200) %>% 
  left_join(df_ts)

ggplot(df_final %>% filter(depth_round < 25))+
  geom_boxplot(aes(x = cluster, y = bbp/FLUO_NPQ, fill = cluster))+
  ylim(0, 0.0025)+
  ylab("bbp/FChla")+
  theme_minimal(base_size = 16)
ggsave("output/boxplot_refractive_index.png", width = 18, height = 10, dpi = 300)


ggplot(df_final %>% filter(depth_round < 25))+
  geom_boxplot(aes(x = date, y = FLUO_NPQ, fill = cluster))+
  ylab("Chla")+
  theme_minimal()

ggplot(df_final %>% filter(depth_round <10) %>% arrange(date))+
  geom_point(aes(x = lon, y =lat , color = cluster), size = 3)+
  scale_color_brewer(palette = "Set1")+
  coord_quickmap()+
  theme_minimal()

ggplot(df_final %>% filter(depth_round < 25))+
  geom_boxplot(aes(x = date, y = bbp, fill = cluster))+
  ylab("Chla")+
  theme_minimal()

write_csv(df_final, "Data/Processed/final_csv_alr6.csv")

