# Oxygen budget (diel net & night O2 change) overlaid on the NCP and NPP
# time series for the dawn/dusk productivity float 3902681.
# Replicates the overlay style of Scripts/ncp_on_canyon_without_F_data.R.
#
# Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/ncp_o2budget_float_3902681.R

library(tidyverse)
library(lubridate)
library(zoo)

# --- NCP (nitrate drawdown, 20-day window) for float 3902681 ----------------
ncp_results <- read_csv("Data/Processed/ncp_float_smoothed_30d.csv") %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  filter(!is.na(NCP))

# loess fit on the NCP series (same settings as the source script)
ncp_results$datenum <- as.numeric(ncp_results$date)
fit_ncp <- loess(NCP ~ datenum, data = ncp_results, span = 0.4,
                 family = "symmetric")
ncp_results$loess_ncp <- predict(fit_ncp)

# --- Diel oxygen budget -----------------------------------------------------
# mtime is MATLAB datenum; convert to date, then collapse to Date for a
# single Date x-axis. O2 fluxes divided by 2 (O2 -> C) as in the source script.
o2_net_change <- read_csv("Data/Processed/O2_float_net_change.csv") %>%
  mutate(date = as.Date(as.POSIXct((mtime - 719529) * 86400,
                                    origin = "1970-01-01", tz = "UTC")),
         datenum = as.numeric(date),
         o2_net_smooth = rollapply(Int_dO2dt_40m_mmol_m2_d / 2, width = 6,
                                   FUN = mean, align = "center", fill = NA)) %>%
  na.omit()
fit <- loess(o2_net_smooth ~ datenum, data = o2_net_change, span = 0.2,
             family = "symmetric")
o2_net_change$o2_net_change <- predict(fit)

o2_float_night_loss <- read_csv("Data/Processed/O2_float_night_loss.csv") %>%
  mutate(date = as.Date(as.POSIXct((mtime - 719529) * 86400,
                                    origin = "1970-01-01", tz = "UTC")),
         datenum = as.numeric(date),
         o2_night_smooth = rollapply(Int_O2_loss_40m_mmol_m2_d / 2, width = 10,
                                     FUN = mean, align = "center", fill = NA)) %>%
  na.omit()
fit <- loess(o2_night_smooth ~ datenum, data = o2_float_night_loss, span = 0.2,
             family = "symmetric")
o2_float_night_loss$o2_night_change <- predict(fit)

# --- NPP (mg C -> mmol C) ---------------------------------------------------
npp <- read_csv("Data/Processed/icb_npp_ncp.csv") %>%
  mutate(npp = npp / 12,
         date_10day = as.Date(date_10day)) %>%
  filter(!is.na(npp))

# --- Cut all series at Nov 2025 ---------------------------------------------
cutoff <- as.Date("2025-11-01")
ncp_results        <- filter(ncp_results, date <= cutoff)
o2_net_change      <- filter(o2_net_change, date <= cutoff)
o2_float_night_loss <- filter(o2_float_night_loss, date <= cutoff)
npp                <- filter(npp, date_10day <= cutoff)

# --- Plot: O2 budget on top of NCP + NPP (raw / rolling-mean series only) ----
ggplot(ncp_results) +
  # NCP
  geom_line(aes(x = date, y = NCP, color = "NCP light"), alpha = 0.9) +

  # O2 NIGHT CHANGE
  geom_line(aes(x = date, y = o2_night_smooth, color = "O2 night light"),
            data = o2_float_night_loss, alpha = 0.9) +

  # O2 NET CHANGE
  geom_line(aes(x = date, y = o2_net_smooth, color = "O2 net light"),
            data = o2_net_change, alpha = 0.9) +

  # NPP
  geom_line(aes(x = date_10day, y = npp, color = "NPP"), data = npp) +

  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +

  scale_color_manual(
    name = "",
    values = c(
      "NCP light"      = "#08519c",
      "O2 night light" = "#cb181d",
      "O2 net light"   = "#006837",
      "NPP"            = "#252525"
    ),
    labels = c(
      "NCP light"      = "NCP (30 days avg)",
      "O2 night light" = "O2 night change (30 days avg)",
      "O2 net light"   = "O2 net change (18 days avg)",
      "NPP"            = "NPP"
    )
  ) +
  labs(y = expression(mmol~C~m^-2~d^-1), x = "Date") +
  theme_bw() +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
  ylim(-100, 175)+
  ggtitle("Float 3902681: NCP and NPP with diel O2 budget (net & night change)")

ggsave("Output/ncp_o2budget_float_3902681_v2.png", width = 11, height = 6, dpi = 200)
cat("saved -> Output/ncp_o2budget_float_3902681.png\n")
