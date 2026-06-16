# Monthly NCP climatology (Irminger Sea, 30-day integration window)
#
# Reads Data/Processed/NCP/ncp_results_irminger.csv (output of compute_ncp()),
# keeps only the 30-day time window, and collapses the 2020-2025 series to a
# 12-month climatology grouped by calendar month, following the convention in
# Scripts/ncp_clean.R (group_by(month) -> mean_ncp / sd_ncp).
#
# NCP is net community production from nitrate drawdown, in mmol C m-2 d-1.
#
# Run from repo root:
#   & "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Scripts/ncp_climatology.R

suppressPackageStartupMessages({
  library(tidyverse)
})

in_file  <- "Data/Processed/NCP/ncp_results_irminger.csv"
out_file <- "Output/ncp_climatology_icb_30d.csv"
window   <- "30 days"

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

ncp <- read_csv(in_file, show_col_types = FALSE) |>
  filter(time_step_label == window) |>
  mutate(month = month(date_grid))

cat(sprintf("Window %s: %d rows, %s..%s\n",
            window, nrow(ncp),
            format(min(ncp$date_grid)), format(max(ncp$date_grid))))

clim <- ncp |>
  group_by(month) |>
  summarise(
    mean_ncp   = mean(NCP, na.rm = TRUE),
    sd_ncp     = sd(NCP, na.rm = TRUE),
    median_ncp = median(NCP, na.rm = TRUE),
    min_ncp    = min(NCP, na.rm = TRUE),
    max_ncp    = max(NCP, na.rm = TRUE),
    n_years    = sum(!is.na(NCP)),
    .groups = "drop"
  ) |>
  # fill missing calendar months with NA rows so the climatology is a full year
  right_join(tibble(month = 1:12), by = "month") |>
  arrange(month) |>
  mutate(
    month_name = factor(month_names[month], levels = month_names),
    se_ncp     = sd_ncp / sqrt(n_years),
    units      = "mmol C m-2 d-1"
  ) |>
  relocate(month, month_name)

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write_csv(clim, out_file)

cat(sprintf("Wrote %s\n", out_file))
print(clim |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = 12)
