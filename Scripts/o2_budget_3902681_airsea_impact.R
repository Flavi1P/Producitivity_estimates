# Simple figures assessing the impact of the air-sea O2 flux correction on the
# 0-40 m budget of float 3902681. Reads the budget CSV written by
# Scripts/o2_budget_3902681.R and overlays, for (1) nighttime O2 loss
# (respiration) and (2) the net O2 budget (NCP), three series:
#   - uncorrected
#   - diffusive-corrected (Wanninkhof 2014)
#   - diffusive + bubble corrected (W2014 + Liang 2013 bubble injection)
#
# Convention reminder (see .claude/o2_budget_3902681_explained.md):
#   - CSV rate_mmol_o2_m2_d            = signed dO2/dt, uncorrected
#   - CSV rate_corr_mmol_o2_m2_d       = rate - F_as       (diffusive only)
#   - CSV rate_corr_total_mmol_o2_m2_d = rate - (F_as + bubble)  (diff + bubble)
#   - night LOSS = -rate (loss is plotted as a positive consumption magnitude)
#   - net change is used directly (signed, no flip)
#
# Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/o2_budget_3902681_airsea_impact.R

library(tidyverse)

in_csv  <- "Data/Processed/o2_budget_3902681.csv"
out_night <- "Output/o2_budget_3902681_airsea_impact_night.png"
out_net   <- "Output/o2_budget_3902681_airsea_impact_net.png"

# centered 8-point rolling mean, same helper as the main script
roll_smooth <- function(x, k = 8) {
  if (length(x) < k) return(x)
  zoo::rollapply(x, width = k, FUN = mean, na.rm = TRUE, fill = NA, align = "center")
}

budget <- read_csv(in_csv, show_col_types = FALSE) |>
  mutate(mtime = as.POSIXct(mtime, tz = "UTC")) |>
  arrange(mtime)

series_levels <- c("uncorrected", "diffusive (W2014)", "diffusive + bubble (W2014 + L13)")

# ---- Figure 1: nighttime O2 loss (respiration), 0-40 m ---------------------
# loss convention: plot -rate so consumption is positive.
night <- budget |>
  filter(type == "night") |>
  transmute(
    mtime,
    uncorrected           = -rate_mmol_o2_m2_d,
    `diffusive (W2014)`   = -rate_corr_mmol_o2_m2_d,
    `diffusive + bubble (W2014 + L13)` = -rate_corr_total_mmol_o2_m2_d
  ) |>
  pivot_longer(-mtime, names_to = "series", values_to = "loss") |>
  group_by(series) |>
  arrange(mtime, .by_group = TRUE) |>
  mutate(loss_smooth = roll_smooth(loss, 8)) |>
  ungroup() |>
  mutate(series = factor(series, levels = series_levels))

p_night <- ggplot(night, aes(mtime, loss_smooth, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = setNames(c("grey45", "#fb6a4a", "#a50f15"), series_levels)) +
  coord_cartesian(ylim = c(-200, 400)) +
  labs(
    title = "Float 3902681 - nighttime O2 loss (respiration), 0-40 m",
    subtitle = "Impact of the air-sea flux correction: diffusive (W2014) then + bubble injection (Liang 2013). 8-pt rolling mean.",
    x = NULL, y = expression("Night O"[2]*" loss  (mmol O"[2]*" m"^-2*" d"^-1*")"),
    colour = NULL
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(out_night, p_night, width = 11, height = 6, dpi = 200)

# ---- Figure 2: net O2 budget (NCP), 0-40 m ---------------------------------
# net change is signed dO2/dt, used directly (no flip).
net <- budget |>
  filter(type == "net") |>
  transmute(
    mtime,
    uncorrected           = rate_mmol_o2_m2_d,
    `diffusive (W2014)`   = rate_corr_mmol_o2_m2_d,
    `diffusive + bubble (W2014 + L13)` = rate_corr_total_mmol_o2_m2_d
  ) |>
  pivot_longer(-mtime, names_to = "series", values_to = "rate") |>
  group_by(series) |>
  arrange(mtime, .by_group = TRUE) |>
  mutate(rate_smooth = roll_smooth(rate, 8)) |>
  ungroup() |>
  mutate(series = factor(series, levels = series_levels))

p_net <- ggplot(net, aes(mtime, rate_smooth, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = setNames(c("grey45", "#74a9cf", "#08519c"), series_levels)) +
  coord_cartesian(ylim = c(-200, 400)) +
  labs(
    title = "Float 3902681 - net O2 budget (NCP), 0-40 m",
    subtitle = "Impact of the air-sea flux correction: diffusive (W2014) then + bubble injection (Liang 2013). 8-pt rolling mean.",
    x = NULL, y = expression("Net dO"[2]*"/dt  (mmol O"[2]*" m"^-2*" d"^-1*")"),
    colour = NULL
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(out_net, p_net, width = 11, height = 6, dpi = 200)

cat("Wrote:\n  ", out_night, "\n  ", out_net, "\n", sep = "")
