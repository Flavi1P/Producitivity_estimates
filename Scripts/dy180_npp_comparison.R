# DY180 — Float CbPM/VGPM vs in-situ CTD NPP
#
# Pipeline:
#   1. Parse Data/DY180_NPP_profiles_compare_4_Flavien.xls           -> long-form CTD NPP
#   2. Parse Data/station_summary_dy180_all.csv                       -> CTD station lat/lon/time
#   3. Same-day matchup: for each (float, target CTD station)         -> float profile match
#   4. Build CbPM-ready profiles (despike, Cphyto, interp 0..199 m)
#   5. Surface PAR (MODIS daily L3m) lookup                            (requires 2024 PAR files)
#   6. Run CbPM and VGPM
#   7. Compare to CTD NPP, save plots and CSV
#
# Run from repo root:
#   & "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" Scripts/dy180_npp_comparison.R

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(ncdf4)
  library(geosphere)
  library(ggplot2)
})

# Paths -----------------------------------------------------------------------
EXCEL_PATH    <- "Data/DY180_NPP_profiles_compare_4_Flavien.xls"
STATION_CSV   <- "Data/station_summary_dy180_all.csv"
FLOAT_DIR     <- "Data/Raw/Floats/Biocarbon"
PAR_DIR       <- "Data/Raw/Satellite/PAR"
CTD_CALIB_DIR <- "Data/Processed/DY180_CTD_Calib"
OUT_INTERMED  <- "Data/Intermediate/Floats"
OUT_DIR       <- "Output"

# Python helper for pvlib clear-sky daily PAR (cbpm_ipar pathway)
PYTHON_BIN      <- "C:/Users/petit/anaconda3/envs/cmts_learn_olci/python.exe"
CLEARSKY_SCRIPT <- "Scripts/clearsky_daily_par.py"

dir.create(OUT_INTERMED, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Step 1 — Parse Excel --------------------------------------------------------
parse_npp_excel <- function(path) {
  raw <- as.data.frame(read_excel(path, sheet = "dpp", col_names = FALSE,
                                  .name_repair = "minimal"))
  n_blocks <- ncol(raw) %/% 4
  blocks <- vector("list", n_blocks)
  for (b in seq_len(n_blocks)) {
    cols <- ((b - 1) * 4 + 1):((b - 1) * 4 + 4)
    block_id   <- as.character(raw[1, cols[1]])              # e.g. "dpp.700"
    raw_id     <- as.character(raw[2, cols[1]])              # e.g. "CTD007", "CTD36"
    sample_id  <- suppressWarnings(as.numeric(raw[2, cols[2]]))
    rep_label  <- as.character(raw[2, cols[3]])
    # numeric data starts row 3
    # Per-block layout (confirmed by user):
    #   col 1 = depth_pe       (depth grid for PE-curve NPP, every 0.5 m)
    #   col 2 = npp_pe         (NPP from PE curve, depth-resolved)
    #   col 3 = depth_incub    (5 standard incubation depths, identical across blocks)
    #   col 4 = carbon_uptake  (carbon uptake from incubation at those 5 depths)
    depth_pe     <- suppressWarnings(as.numeric(raw[3:nrow(raw), cols[1]]))
    npp_pe       <- suppressWarnings(as.numeric(raw[3:nrow(raw), cols[2]]))
    depth_incub  <- suppressWarnings(as.numeric(raw[3:nrow(raw), cols[3]]))
    carbon_upt   <- suppressWarnings(as.numeric(raw[3:nrow(raw), cols[4]]))
    stn_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", raw_id)))
    common <- list(
      block_id  = block_id, ctd_raw = raw_id, stn_num = stn_num,
      ctd_id    = ifelse(is.na(stn_num), NA_character_, sprintf("CTD%03d", stn_num)),
      sample_id = sample_id, rep_label = rep_label
    )
    pe <- tibble(!!!common, src = "pe",
                 depth_m = abs(depth_pe), value = npp_pe) |>
          filter(!is.na(depth_m), !is.na(value))
    inc <- tibble(!!!common, src = "incub",
                  depth_m = abs(depth_incub), value = carbon_upt) |>
           filter(!is.na(depth_m), !is.na(value))
    blocks[[b]] <- bind_rows(pe, inc)
  }
  bind_rows(blocks)
}

excel_long <- parse_npp_excel(EXCEL_PATH)
cat("== Step 1: Excel parsed ==\n")
cat("  rows:", nrow(excel_long), "  blocks:", length(unique(excel_long$block_id)), "\n")
cat("  unique stations:", paste(unique(excel_long$ctd_id), collapse = ", "), "\n")

# Step 2 — Parse station_summary CSV -----------------------------------------
parse_deg_min <- function(s) {
  # "59 11.382 N" or "022 25.322 W"  ->  signed decimal degrees
  s <- str_trim(s)
  parts <- str_split_fixed(s, "\\s+", 3)
  deg <- suppressWarnings(as.numeric(parts[, 1]))
  minutes <- suppressWarnings(as.numeric(parts[, 2]))
  hemi <- toupper(parts[, 3])
  dec <- deg + minutes / 60
  ifelse(hemi %in% c("S", "W"), -dec, dec)
}

parse_station_csv <- function(path) {
  lines <- read_lines(path)
  # Header is line 1; data lines have a station number in col 1.
  # Each row is comma-separated; station rows have a non-blank first field.
  rows <- lapply(lines[-1], function(L) str_split(L, ",", simplify = TRUE)[1, ])
  ncols <- max(vapply(rows, length, integer(1)))
  mat <- do.call(rbind, lapply(rows, function(r) {
    length(r) <- ncols; r
  }))
  base_names <- c("stn", "datetime", "lat_str", "lon_str", "cordep", "maxd",
                  "minalt", "resid", "maxw", "maxp", "ndpths", "nsal", "noxy",
                  "nnut", "nco2")
  cn <- character(ncols)
  cn[seq_len(min(ncols, length(base_names)))] <-
    base_names[seq_len(min(ncols, length(base_names)))]
  if (ncols > length(base_names))
    cn[(length(base_names) + 1):ncols] <-
      paste0("extra_", seq_len(ncols - length(base_names)))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- cn
  df |>
    mutate(across(everything(), str_trim)) |>
    filter(stn != "" & !is.na(stn) & stn != "stn") |>
    mutate(stn_num  = suppressWarnings(as.integer(stn)),
           ctd_time = ymd_hm(datetime, tz = "UTC"),
           lat      = parse_deg_min(lat_str),
           lon      = parse_deg_min(lon_str),
           ctd_id   = sprintf("CTD%03d", stn_num)) |>
    filter(!is.na(stn_num))
}

stations_all <- parse_station_csv(STATION_CSV)
cat("== Step 2: station summary parsed ==\n")
cat("  rows:", nrow(stations_all), "  date range:", format(range(stations_all$ctd_time)), "\n")

excel_stations <- unique(na.omit(excel_long$ctd_id))
# One rep_label per ctd_id (first non-NA from any block of that station).
# CTD036 has two blocks but both share rep "S4", so collapsing is safe here;
# block-level rep_label is preserved in excel_long for facet titles further down.
station_rep <- excel_long |>
  filter(!is.na(ctd_id), !is.na(rep_label), rep_label != "NA") |>
  distinct(ctd_id, rep_label) |>
  group_by(ctd_id) |>
  summarise(rep_label = paste(unique(rep_label), collapse = "/"), .groups = "drop")

target_stations <- stations_all |>
  filter(ctd_id %in% excel_stations) |>
  select(ctd_id, stn_num, ctd_time, lat, lon) |>
  left_join(station_rep, by = "ctd_id") |>
  mutate(stn_label = sprintf("%s (%s)", ctd_id, rep_label))
cat("  matched targets:", nrow(target_stations), "/", length(excel_stations), "\n")
print(target_stations)

# Step 3 — Same-day float ↔ CTD matchup --------------------------------------
read_float_meta <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  juld <- ncvar_get(nc, "JULD")
  lat  <- ncvar_get(nc, "LATITUDE")
  lon  <- ncvar_get(nc, "LONGITUDE")
  wmo  <- as.integer(sub("_Sprof.nc$", "", basename(path)))
  tibble(float_wmo = wmo,
         prof_idx  = seq_along(juld),
         prof_time = as.POSIXct(as.numeric(juld) * 86400,
                                origin = "1950-01-01", tz = "UTC"),
         prof_lat  = as.numeric(lat),
         prof_lon  = as.numeric(lon))
}

float_files <- list.files(FLOAT_DIR, pattern = "_Sprof\\.nc$", full.names = TRUE)
float_meta  <- map_dfr(float_files, read_float_meta)
cat("== Step 3a: float metadata ==\n")
cat("  total profiles:", nrow(float_meta), "\n")
print(float_meta |> group_by(float_wmo) |>
        summarise(n = n(), tmin = min(prof_time), tmax = max(prof_time)))

# For each (float, target station), find profiles on the same calendar date.
# If multiple matches, pick the spatially closest (haversine).
matchup <- target_stations |>
  mutate(ctd_date = as.Date(ctd_time)) |>
  tidyr::expand_grid(float_wmo = unique(float_meta$float_wmo)) |>
  left_join(float_meta, by = "float_wmo", relationship = "many-to-many") |>
  filter(as.Date(prof_time) == ctd_date) |>
  rowwise() |>
  mutate(dist_km = distHaversine(c(lon, lat), c(prof_lon, prof_lat)) / 1000) |>
  ungroup() |>
  group_by(ctd_id, float_wmo) |>
  slice_min(dist_km, n = 1, with_ties = FALSE) |>
  ungroup()

cat("== Step 3b: same-day matchup ==\n")
cat("  matchups:", nrow(matchup), "\n")
print(matchup |> select(ctd_id, ctd_time, float_wmo, prof_idx, prof_time,
                        prof_lat, prof_lon, dist_km))

write_csv(matchup, file.path(OUT_INTERMED, "dy180_ctd_float_matchup.csv"))
cat("  -> saved", file.path(OUT_INTERMED, "dy180_ctd_float_matchup.csv"), "\n")

if (nrow(matchup) == 0) {
  cat("\n*** NO SAME-DAY MATCHES. Stopping here. Consider widening the time window. ***\n")
  quit(status = 0)
}

# Step 4 — Build CbPM-ready profiles -----------------------------------------
# For each matched (float, profile), extract PRES, CHLA, BBP700, TEMP; despike;
# derive Cphyto = 12128 * (bbp700 / (470/400)) + 0.59 (Graff et al. 2015,
# matching Scripts/3d_products.ipynb cell 42); interpolate to 0..199 m.

source("Scripts/cbpm_r/despike.r")
source("Scripts/cbpm_r/cbpm_argo.r")
source("Scripts/cbpm_r/daylength.r")
source("Scripts/cbpm_r/vgpm.r")

read_float_profile_vars <- function(float_path, prof_idx) {
  nc <- nc_open(float_path)
  on.exit(nc_close(nc))
  pick <- function(name_adj, name_raw) {
    v_adj <- tryCatch(ncvar_get(nc, name_adj), error = function(e) NULL)
    v_raw <- tryCatch(ncvar_get(nc, name_raw), error = function(e) NULL)
    # 2-D var: [N_LEVELS, N_PROF]
    if (!is.null(v_adj)) {
      col <- v_adj[, prof_idx]
      if (all(is.na(col)) && !is.null(v_raw)) col <- v_raw[, prof_idx]
      col
    } else if (!is.null(v_raw)) {
      v_raw[, prof_idx]
    } else NA_real_
  }
  list(
    pres = pick("PRES_ADJUSTED", "PRES"),
    chl  = pick("CHLA_ADJUSTED", "CHLA"),
    bbp  = pick("BBP700_ADJUSTED", "BBP700"),
    temp = pick("TEMP_ADJUSTED", "TEMP")
  )
}

# Read the float's near-surface DOWNWELLING_PAR for a given profile.
# Returns the shallowest valid sample within the upper 10 m of the profile,
# units umol m^-2 s^-1 (Argo BGC convention). NA if none available.
read_float_surface_par <- function(float_path, prof_idx) {
  nc <- nc_open(float_path)
  on.exit(nc_close(nc))
  pick <- function(name_adj, name_raw) {
    v_adj <- tryCatch(ncvar_get(nc, name_adj), error = function(e) NULL)
    v_raw <- tryCatch(ncvar_get(nc, name_raw), error = function(e) NULL)
    if (!is.null(v_adj)) {
      col <- v_adj[, prof_idx]
      if (all(is.na(col)) && !is.null(v_raw)) col <- v_raw[, prof_idx]
      col
    } else if (!is.null(v_raw)) {
      v_raw[, prof_idx]
    } else NA_real_
  }
  pres <- as.numeric(pick("PRES_ADJUSTED", "PRES"))
  par  <- as.numeric(pick("DOWNWELLING_PAR_ADJUSTED", "DOWNWELLING_PAR"))
  if (length(par) == 0 || all(is.na(par))) return(NA_real_)
  ok <- is.finite(pres) & is.finite(par) & pres <= 10
  if (!any(ok)) return(NA_real_)
  par[ok][which.min(pres[ok])]
}

# Batch clear-sky daily-PAR computation via the Python helper.
# `rows` is a tibble with columns id, datetime_utc, lat, lon, par_meas.
# Returns a named numeric vector keyed by id (mol photons m^-2 d^-1).
clearsky_daily_par <- function(rows) {
  rows <- dplyr::filter(rows, is.finite(par_meas), is.finite(lat), is.finite(lon))
  if (nrow(rows) == 0) return(setNames(numeric(0), character(0)))
  in_csv  <- tempfile(fileext = ".csv")
  out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(in_csv, out_csv)), add = TRUE)
  rows |>
    mutate(datetime_utc = format(datetime_utc, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) |>
    write_csv(in_csv)
  status <- system2(PYTHON_BIN,
                    c(CLEARSKY_SCRIPT, "--in", in_csv, "--out", out_csv),
                    stdout = "", stderr = "")
  if (status != 0 || !file.exists(out_csv))
    stop("clearsky_daily_par.py failed (status ", status, ")")
  out <- read_csv(out_csv, show_col_types = FALSE)
  setNames(as.numeric(out$daily_par_mol_m2_d), as.character(out$id))
}

build_cbpm_inputs <- function(prof_vars) {
  d <- tibble(pres = as.numeric(prof_vars$pres),
              chl  = as.numeric(prof_vars$chl),
              bbp  = as.numeric(prof_vars$bbp),
              temp = as.numeric(prof_vars$temp)) |>
    filter(!is.na(pres)) |>
    arrange(pres)

  if (nrow(d) < 5 || all(is.na(d$chl)) || all(is.na(d$bbp)))
    return(NULL)

  # Despike CHL and BBP (Briggs 2011, 7-point median baseline)
  if (sum(!is.na(d$chl)) >= 7)
    d$chl_clean <- despike(d$chl, window_size = 7)$baseline
  else
    d$chl_clean <- d$chl
  if (sum(!is.na(d$bbp)) >= 7)
    d$bbp_clean <- despike(d$bbp, window_size = 7)$baseline
  else
    d$bbp_clean <- d$bbp

  # Derive Cphyto (Graff et al. 2015, also used in 3d_products.ipynb cell 42)
  d$bbp470  <- d$bbp_clean / (470 / 400)
  d$cphyto  <- 12128 * d$bbp470 + 0.59

  # Surface SST = shallowest valid temperature in upper 10 m
  sst <- d |> filter(pres <= 10, !is.na(temp)) |> pull(temp) |> head(1)
  if (length(sst) == 0) sst <- NA_real_

  # Interpolate Chl, bbp, Cphyto onto 0..199 m grid (CbPM contract)
  grid <- 0:199
  chl_ok  <- !is.na(d$pres) & !is.na(d$chl_clean)
  bbp_ok  <- !is.na(d$pres) & !is.na(d$bbp_clean)
  cph_ok  <- !is.na(d$pres) & !is.na(d$cphyto)
  if (sum(chl_ok) < 2 || sum(cph_ok) < 2) return(NULL)
  chl_z    <- approx(d$pres[chl_ok], d$chl_clean[chl_ok], xout = grid, rule = 2)$y
  bbp_z    <- approx(d$pres[bbp_ok], d$bbp_clean[bbp_ok], xout = grid, rule = 2)$y
  cphyto_z <- approx(d$pres[cph_ok], d$cphyto[cph_ok],     xout = grid, rule = 2)$y
  list(chl_z = chl_z, bbp_z = bbp_z, cphyto_z = cphyto_z,
       sst = sst, depth_grid = grid)
}

float_file_for <- function(wmo) file.path(FLOAT_DIR, sprintf("%d_Sprof.nc", wmo))

inputs <- vector("list", nrow(matchup))
for (i in seq_len(nrow(matchup))) {
  pv <- read_float_profile_vars(float_file_for(matchup$float_wmo[i]),
                                matchup$prof_idx[i])
  inputs[[i]] <- build_cbpm_inputs(pv)
}
n_inputs_ok <- sum(!vapply(inputs, is.null, logical(1)))
cat("== Step 4: CbPM-ready inputs built ==\n")
cat("  matchups with usable Chl+bbp profiles:", n_inputs_ok, "/", nrow(matchup), "\n")

# Step 4b — Read calibrated CTD profiles and build CbPM inputs ---------------
# CSVs at Data/Processed/DY180_CTD_Calib/CTD_NNN_dn_calibrated.csv:
#   Depth_Center, latitude, longitude, t090C, calib_chl_qc_flC_X12_zero, bbp
read_calib_ctd <- function(stn_num) {
  fp <- file.path(CTD_CALIB_DIR, sprintf("CTD_%03d_dn_calibrated.csv", stn_num))
  if (!file.exists(fp)) return(NULL)
  d <- suppressWarnings(read_csv(fp, show_col_types = FALSE))
  tibble(pres = as.numeric(d$Depth_Center),
         chl  = as.numeric(d$calib_chl_qc_flC_X12_zero),
         bbp  = as.numeric(d$bbp),
         temp = as.numeric(d$t090C))
}

ctd_inputs <- vector("list", nrow(target_stations))
for (i in seq_len(nrow(target_stations))) {
  pv <- read_calib_ctd(target_stations$stn_num[i])
  ctd_inputs[[i]] <- if (is.null(pv)) NULL else build_cbpm_inputs(pv)
}
n_ctd_ok <- sum(!vapply(ctd_inputs, is.null, logical(1)))
cat("== Step 4b: CbPM inputs from calibrated CTD ==\n")
cat("  stations with usable Chl+bbp profiles:", n_ctd_ok, "/",
    nrow(target_stations), "\n")

# Step 5 — Surface PAR helper -------------------------------------------------
# Reads MODIS-Aqua daily L3m PAR (AQUA_MODIS.YYYYMMDD.L3m.DAY.PAR.x_par.nc)
# and returns the value at the nearest grid cell (E m^-2 d^-1).
get_daily_par <- function(date, lat, lon) {
  ymd_str <- format(as.Date(date), "%Y%m%d")
  candidates <- c(
    sprintf("AQUA_MODIS.%s.L3m.DAY.PAR.x_par.nc", ymd_str),
    sprintf("AQUA_MODIS.%s.L3m.DAY.PAR.par.nc",  ymd_str)
  )
  for (fn in candidates) {
    fp <- file.path(PAR_DIR, fn)
    if (file.exists(fp)) {
      nc <- nc_open(fp)
      on.exit(nc_close(nc), add = TRUE)
      par_name <- if ("par" %in% names(nc$var)) "par" else
                  if ("x_par" %in% names(nc$var)) "x_par" else names(nc$var)[1]
      lon_g <- ncvar_get(nc, "lon"); lat_g <- ncvar_get(nc, "lat")
      iy <- which.min(abs(lat_g - lat))
      ix <- which.min(abs(lon_g - lon))
      val <- ncvar_get(nc, par_name, start = c(ix, iy), count = c(1, 1))
      val <- as.numeric(val)
      if (!is.finite(val) || val < 0) val <- NA_real_
      return(val)
    }
  }
  return(NA_real_)
}

# Step 6 — Run CbPM and VGPM --------------------------------------------------
# Pre-compute clear-sky daily PAR (cbpm_ipar pathway): one batched call to the
# pvlib helper for every matchup that has a finite float surface PAR.
cs_inputs <- tibble(
  id            = as.character(seq_len(nrow(matchup))),
  datetime_utc = matchup$prof_time,
  lat           = matchup$prof_lat,
  lon           = matchup$prof_lon,
  par_meas      = vapply(seq_len(nrow(matchup)), function(i) {
    read_float_surface_par(float_file_for(matchup$float_wmo[i]),
                           matchup$prof_idx[i])
  }, numeric(1))
)
cat("== Step 6 prep: float surface PAR for cbpm_ipar ==\n")
cat("  matchups with surface PAR:",
    sum(is.finite(cs_inputs$par_meas)), "/", nrow(cs_inputs), "\n")
daily_par_lookup <- if (any(is.finite(cs_inputs$par_meas))) {
  clearsky_daily_par(cs_inputs)
} else setNames(numeric(0), character(0))
cat("  matchups with finite clear-sky daily PAR:",
    sum(is.finite(daily_par_lookup)), "/", nrow(cs_inputs), "\n")

results <- vector("list", nrow(matchup))
for (i in seq_len(nrow(matchup))) {
  inp <- inputs[[i]]
  if (is.null(inp)) next
  dt <- as.POSIXlt(matchup$ctd_time[i], tz = "UTC")
  yr <- dt$year + 1900; mo <- dt$mon + 1; dy <- dt$mday
  lat_p <- matchup$prof_lat[i]; lon_p <- matchup$prof_lon[i]

  irr <- get_daily_par(matchup$ctd_time[i], lat_p, lon_p)
  if (!is.finite(irr)) {
    cat(sprintf("  [%s float %d] PAR missing for %s -> skip CbPM/VGPM\n",
                matchup$ctd_id[i], matchup$float_wmo[i],
                format(as.Date(matchup$ctd_time[i]))))
    next
  }
  dayL <- daylength(yr, mo, dy, lat_p)

  cbpm <- tryCatch(
    cbpm_argo(inp$chl_z, inp$cphyto_z, irr, yr, mo, dy, lat_p),
    error = function(e) {
      cat("  CbPM error:", conditionMessage(e), "\n"); NULL
    })
  vgpm <- opp_befa(inp$chl_z[1], irr, inp$sst, dayL)

  if (is.null(cbpm)) next
  pp_z <- cbpm$pp_z
  z_eu <- if (is.na(cbpm$mzeu)) length(pp_z) - 1 else cbpm$mzeu
  cbpm_int <- sum(pp_z[1:(z_eu + 1)], na.rm = TRUE)         # mg C m^-2 d^-1

  # Parallel CbPM run with pvlib clear-sky daily PAR scaled from float surface PAR.
  irr_ipar      <- unname(daily_par_lookup[as.character(i)])
  par_meas      <- cs_inputs$par_meas[i]
  pp_z_ipar     <- rep(NA_real_, length(pp_z))
  z_eu_ipar     <- NA_integer_
  cbpm_ipar_int <- NA_real_
  if (length(irr_ipar) == 1 && is.finite(irr_ipar)) {
    cbpm_ipar_run <- tryCatch(
      cbpm_argo(inp$chl_z, inp$cphyto_z, irr_ipar, yr, mo, dy, lat_p),
      error = function(e) {
        cat("  CbPM (ipar) error:", conditionMessage(e), "\n"); NULL
      })
    if (!is.null(cbpm_ipar_run)) {
      pp_z_ipar     <- cbpm_ipar_run$pp_z
      z_eu_ipar     <- if (is.na(cbpm_ipar_run$mzeu)) length(pp_z_ipar) - 1
                       else cbpm_ipar_run$mzeu
      cbpm_ipar_int <- sum(pp_z_ipar[1:(z_eu_ipar + 1)], na.rm = TRUE)
    }
  }

  results[[i]] <- tibble(
    matchup_idx   = i,
    ctd_id        = matchup$ctd_id[i],
    rep_label     = matchup$rep_label[i],
    stn_label     = matchup$stn_label[i],
    float_wmo     = matchup$float_wmo[i],
    irr_par       = irr,
    par_meas_float = par_meas,
    irr_ipar      = if (length(irr_ipar) == 1) irr_ipar else NA_real_,
    sst           = inp$sst,
    dayL          = dayL,
    z_eu_cbpm     = z_eu,
    z_eu_cbpm_ipar = z_eu_ipar,
    cbpm_int      = cbpm_int,
    cbpm_ipar_int = cbpm_ipar_int,
    vgpm_int      = vgpm$npp,
    vgpm_pb_opt   = vgpm$pb_opt,
    vgpm_z_eu     = vgpm$z_eu,
    pp_z          = list(pp_z),
    pp_z_ipar     = list(pp_z_ipar),
    chl_z         = list(inp$chl_z),
    bbp_z         = list(inp$bbp_z)
  )
}
results_df <- bind_rows(results)
cat("== Step 6: CbPM + VGPM computed ==\n")
cat("  successful runs:", nrow(results_df), "/", nrow(matchup), "\n")
if (nrow(results_df) > 0) print(
  results_df |> select(-pp_z)
)

# Step 6b — Run CbPM/VGPM on calibrated CTD profiles --------------------------
ctd_results <- vector("list", nrow(target_stations))
for (i in seq_len(nrow(target_stations))) {
  inp <- ctd_inputs[[i]]
  if (is.null(inp)) next
  dt <- as.POSIXlt(target_stations$ctd_time[i], tz = "UTC")
  yr <- dt$year + 1900; mo <- dt$mon + 1; dy <- dt$mday
  lat_p <- target_stations$lat[i]; lon_p <- target_stations$lon[i]
  irr <- get_daily_par(target_stations$ctd_time[i], lat_p, lon_p)
  if (!is.finite(irr)) next
  dayL <- daylength(yr, mo, dy, lat_p)
  cbpm <- tryCatch(cbpm_argo(inp$chl_z, inp$cphyto_z, irr, yr, mo, dy, lat_p),
                   error = function(e) NULL)
  if (is.null(cbpm)) next
  vgpm <- opp_befa(inp$chl_z[1], irr, inp$sst, dayL)
  pp_z <- cbpm$pp_z
  z_eu <- if (is.na(cbpm$mzeu)) length(pp_z) - 1 else cbpm$mzeu

  # Sensitivity: divide bbp by N before the Cphyto conversion. Cphyto is linear
  # in bbp (Graff et al. 2015: Cphyto = 12128 * bbp470 + 0.59 with bbp470 = bbp/(470/400)).
  run_cbpm_bbp_div <- function(div) {
    cphyto_div <- 12128 * (inp$bbp_z / (470 / 400) / div) + 0.59
    out <- tryCatch(
      cbpm_argo(inp$chl_z, cphyto_div, irr, yr, mo, dy, lat_p),
      error = function(e) NULL)
    pp <- if (is.null(out)) rep(NA_real_, length(pp_z)) else out$pp_z
    zeu_div <- if (is.null(out) || is.na(out$mzeu))
                 length(pp) - 1 else out$mzeu
    int <- if (is.null(out)) NA_real_
           else sum(pp[1:(zeu_div + 1)], na.rm = TRUE)
    list(pp_z = pp, z_eu = zeu_div, cbpm_int = int)
  }
  bbp3  <- run_cbpm_bbp_div(3)
  bbp10 <- run_cbpm_bbp_div(10)

  ctd_results[[i]] <- tibble(
    ctd_id      = target_stations$ctd_id[i],
    rep_label   = target_stations$rep_label[i],
    stn_label   = target_stations$stn_label[i],
    irr_par     = irr,
    sst         = inp$sst,
    dayL        = dayL,
    z_eu_cbpm   = z_eu,
    z_eu_cbpm_bbp3  = bbp3$z_eu,
    z_eu_cbpm_bbp10 = bbp10$z_eu,
    cbpm_int_ctd       = sum(pp_z[1:(z_eu + 1)], na.rm = TRUE),
    cbpm_int_ctd_bbp3  = bbp3$cbpm_int,
    cbpm_int_ctd_bbp10 = bbp10$cbpm_int,
    vgpm_int_ctd = vgpm$npp,
    pp_z        = list(pp_z),
    pp_z_bbp3   = list(bbp3$pp_z),
    pp_z_bbp10  = list(bbp10$pp_z),
    chl_z       = list(inp$chl_z),
    bbp_z       = list(inp$bbp_z)
  )
}
ctd_results_df <- bind_rows(ctd_results)
cat("== Step 6b: CbPM + VGPM on calibrated CTDs ==\n")
cat("  successful runs:", nrow(ctd_results_df), "/", nrow(target_stations), "\n")
if (nrow(ctd_results_df) > 0) print(
  ctd_results_df |> select(ctd_id, rep_label, irr_par, sst, z_eu_cbpm,
                           cbpm_int_ctd, cbpm_int_ctd_bbp3,
                           cbpm_int_ctd_bbp10, vgpm_int_ctd)
)

# Step 7 — Compare to CTD NPP and save outputs --------------------------------
# Two CTD references per station (per Excel block):
#   PE-curve NPP    -> depth-resolved profile  -> depth-integrated for scatter
#   Incubation C-uptake -> 5 surface depths    -> sampled CbPM at those depths
trapz <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 2) NA_real_ else sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

if (nrow(results_df) > 0) {
  # CTD PE-curve NPP, depth-integrated 0..max-depth (mg C m^-2 d^-1)
  ctd_pe_int <- excel_long |>
    filter(src == "pe", !is.na(ctd_id)) |>
    arrange(ctd_id, block_id, depth_m) |>
    group_by(ctd_id, block_id) |>
    summarise(ctd_pe_int = trapz(depth_m, value), .groups = "drop")

  # CTD incubation carbon uptake, integrated over the 5 standard depths
  ctd_inc_int <- excel_long |>
    filter(src == "incub", !is.na(ctd_id)) |>
    arrange(ctd_id, block_id, depth_m) |>
    group_by(ctd_id, block_id) |>
    summarise(ctd_inc_int = trapz(depth_m, value),
              max_inc_depth = max(depth_m, na.rm = TRUE), .groups = "drop")

  # CbPM integrated over the same surface depth range as the incubations,
  # for an apples-to-apples comparison with the bottle data.
  inc_depth_per_station <- ctd_inc_int |>
    group_by(ctd_id) |>
    summarise(max_inc_depth = max(max_inc_depth, na.rm = TRUE), .groups = "drop")

  cbpm_surface_int <- results_df |>
    select(matchup_idx, ctd_id, pp_z) |>
    left_join(inc_depth_per_station, by = "ctd_id") |>
    mutate(cbpm_inc_int = purrr::map2_dbl(pp_z, max_inc_depth, function(pp, zmax) {
      if (is.na(zmax)) zmax <- 10
      zmax <- min(ceiling(zmax), length(pp) - 1)
      trapz(0:zmax, pp[1:(zmax + 1)])
    })) |>
    select(matchup_idx, cbpm_inc_int)

  compare_df <- results_df |>
    select(-pp_z, -pp_z_ipar, -chl_z, -bbp_z) |>
    left_join(ctd_pe_int,  by = "ctd_id", relationship = "many-to-many") |>
    left_join(ctd_inc_int, by = c("ctd_id", "block_id")) |>
    left_join(cbpm_surface_int, by = "matchup_idx") |>
    left_join(ctd_results_df |>
                select(ctd_id, cbpm_int_ctd,
                       cbpm_int_ctd_bbp3, cbpm_int_ctd_bbp10,
                       vgpm_int_ctd),
              by = "ctd_id")
  
  

  write_csv(compare_df,
            file.path(OUT_DIR, "dy180_npp_matchup_table.csv"))
  cat("  -> wrote", file.path(OUT_DIR, "dy180_npp_matchup_table.csv"), "\n")

  # Per-station rep_label lookup for plot labelling, augmented with distance
  # to each matched float profile (one short line per float).
  dist_summary <- matchup |>
    arrange(ctd_id, float_wmo) |>
    group_by(ctd_id) |>
    summarise(dist_str = paste(sprintf("%d: %.0f km", float_wmo, dist_km),
                               collapse = "\n"),
              .groups = "drop")

  station_label_lookup <- target_stations |>
    select(ctd_id, stn_label, rep_label) |>
    left_join(dist_summary, by = "ctd_id") |>
    mutate(stn_label = ifelse(is.na(dist_str),
                              stn_label,
                              paste(stn_label, dist_str, sep = "\n")))

  # Scatter A: column-integrated CbPM/VGPM vs depth-integrated PE-curve NPP
  sa <- compare_df |>
    select(ctd_id, block_id, float_wmo, ctd_pe_int, cbpm_int, vgpm_int) |>
    pivot_longer(c(cbpm_int, vgpm_int), names_to = "method", values_to = "float_npp") |>
    mutate(method = recode(method, cbpm_int = "CbPM (0..Z_eu)", vgpm_int = "VGPM")) |>
    left_join(station_label_lookup, by = "ctd_id")
  rng_a <- range(c(sa$ctd_pe_int, sa$float_npp), na.rm = TRUE)
  n_stns_a <- length(unique(sa$ctd_id))
  n_matchups_a <- length(unique(paste(sa$ctd_id, sa$float_wmo)))
  ggplot(sa, aes(ctd_pe_int, float_npp,
                 colour = factor(float_wmo), shape = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = rep_label), nudge_y = 60, size = 2.8,
              show.legend = FALSE) +
    coord_equal(xlim = rng_a, ylim = rng_a) +
    labs(x = "CTD PE-curve depth-integrated NPP (mg C m^-2 d^-1)",
         y = "Float-derived NPP (mg C m^-2 d^-1)",
         colour = "Float WMO", shape = "Method",
         title = "DY180: float CbPM/VGPM vs CTD PE-curve NPP",
         subtitle = sprintf("%d float-station matchups across %d CTD stations (label = rep)",
                            n_matchups_a, n_stns_a)) +
    theme_bw()
  ggsave(file.path(OUT_DIR, "dy180_npp_scatter_cbpm_vs_pe.png"),
         width = 7.5, height = 6.5, dpi = 150)
  cat("  -> wrote", file.path(OUT_DIR, "dy180_npp_scatter_cbpm_vs_pe.png"), "\n")

  # Scatter B: surface-integrated CbPM vs incubation carbon uptake
  sb <- compare_df |>
    distinct(ctd_id, float_wmo, ctd_inc_int, cbpm_inc_int) |>
    filter(!is.na(ctd_inc_int), !is.na(cbpm_inc_int)) |>
    left_join(station_label_lookup, by = "ctd_id")
  rng_b <- range(c(sb$ctd_inc_int, sb$cbpm_inc_int), na.rm = TRUE)
  ggplot(sb, aes(ctd_inc_int, cbpm_inc_int, colour = factor(float_wmo))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = rep_label), nudge_y = 50, size = 2.8,
              show.legend = FALSE) +
    coord_equal(xlim = rng_b, ylim = rng_b) +
    labs(x = "CTD incubation C-uptake, integrated over 5 surface depths (mg C m^-2 d^-1)",
         y = "CbPM integrated over same depth range (mg C m^-2 d^-1)",
         colour = "Float WMO",
         title = "DY180: float CbPM vs CTD incubation carbon uptake",
         subtitle = sprintf("%d matchups (label = rep)", nrow(sb))) +
    theme_bw()
  ggsave(file.path(OUT_DIR, "dy180_npp_scatter_cbpm_vs_incubation.png"),
         width = 7.5, height = 6.5, dpi = 150)
  cat("  -> wrote", file.path(OUT_DIR, "dy180_npp_scatter_cbpm_vs_incubation.png"), "\n")

  # Profile panels: CTD PE-curve profile vs CbPM pp_z + incubation points.
  # Use the augmented stn_label (with per-float distances) from station_label_lookup.
  attach_label <- function(df) {
    df |> select(-any_of("stn_label")) |>
      left_join(station_label_lookup |> select(ctd_id, stn_label),
                by = "ctd_id")
  }

  unnest_profile <- function(df, col, value_name) {
    df |>
      mutate(prof = purrr::map(.data[[col]],
                               ~ tibble(depth_m = 0:(length(.x) - 1), value = .x))) |>
      select(ctd_id, float_wmo, prof) |>
      tidyr::unnest(prof) |>
      rename(!!value_name := value) |>
      attach_label()
  }
  pp_profiles  <- unnest_profile(results_df, "pp_z",  "pp")
  chl_profiles <- unnest_profile(results_df, "chl_z", "chl")
  bbp_profiles <- unnest_profile(results_df, "bbp_z", "bbp")

  unnest_ctd <- function(df, col, value_name) {
    df |>
      mutate(prof = purrr::map(.data[[col]],
                               ~ tibble(depth_m = 0:(length(.x) - 1), value = .x))) |>
      select(ctd_id, prof) |>
      tidyr::unnest(prof) |>
      rename(!!value_name := value) |>
      attach_label()
  }
  ctd_pp_profiles  <- if (nrow(ctd_results_df)) unnest_ctd(ctd_results_df, "pp_z",  "pp")  else NULL
  ctd_chl_profiles <- if (nrow(ctd_results_df)) unnest_ctd(ctd_results_df, "chl_z", "chl") else NULL
  ctd_bbp_profiles <- if (nrow(ctd_results_df)) unnest_ctd(ctd_results_df, "bbp_z", "bbp") else NULL

  pe_profiles <- excel_long |>
    filter(src == "pe", !is.na(ctd_id)) |>
    select(ctd_id, block_id, depth_m, value) |>
    left_join(station_label_lookup, by = "ctd_id")
  inc_points <- excel_long |>
    filter(src == "incub", !is.na(ctd_id)) |>
    select(ctd_id, block_id, depth_m, value) |>
    left_join(station_label_lookup, by = "ctd_id")

  npp_plot <- ggplot() +
    geom_path(data = pe_profiles,
              aes(value, depth_m, group = block_id),
              colour = "black", alpha = 0.6) +
    geom_point(data = inc_points,
               aes(value, depth_m, group = block_id),
               colour = "black", shape = 17, size = 2) +
    geom_path(data = pp_profiles,
              aes(pp, depth_m, colour = factor(float_wmo))) +
    scale_y_reverse(limits = c(60, 0)) +
    facet_wrap(~ stn_label, scales = "free_x") +
    labs(x = "NPP / C-uptake (mg C m^-3 d^-1)",
         y = "Depth (m)",
         colour = "Float WMO",
         title = "DY180 per-station NPP profiles",
         subtitle = "Black line = CTD PE-curve; black triangles = incubation; grey dashed = CbPM-on-CTD; coloured = float CbPM") +
    theme_bw()
  if (!is.null(ctd_pp_profiles)) {
    npp_plot <- npp_plot +
      geom_path(data = ctd_pp_profiles,
                aes(pp, depth_m), colour = "grey30", linetype = "22", linewidth = 0.7)
  }
  ggsave(file.path(OUT_DIR, "dy180_npp_profiles_per_station.png"),
         plot = npp_plot, width = 11, height = 8, dpi = 150)
  cat("  -> wrote", file.path(OUT_DIR, "dy180_npp_profiles_per_station.png"), "\n")

  # Variant: overlay the cbpm_ipar profile (CbPM with pvlib clear-sky daily PAR
  # scaled from the float's instantaneous surface DOWNWELLING_PAR). Solid coloured
  # = MODIS-PAR CbPM (original), dotted coloured = clear-sky-PAR CbPM.
  ipar_profiles <- unnest_profile(results_df, "pp_z_ipar", "pp") |>
    dplyr::filter(is.finite(pp))
  if (nrow(ipar_profiles) > 0) {
    npp_plot_ipar <- npp_plot +
      geom_path(data = ipar_profiles,
                aes(pp, depth_m, colour = factor(float_wmo)),
                linetype = "dotted", linewidth = 0.8) +
      labs(subtitle = paste(
        "Black line = CTD PE-curve; black triangles = incubation;",
        "grey dashed = CbPM-on-CTD;",
        "coloured solid = float CbPM (MODIS daily PAR);",
        "coloured dotted = float CbPM (pvlib clear-sky daily PAR from float DOWNWELLING_PAR)",
        sep = "\n"))
    ggsave(file.path(OUT_DIR, "dy180_npp_profiles_per_station_with_ipar.png"),
           plot = npp_plot_ipar, width = 11, height = 8, dpi = 150)
    cat("  -> wrote",
        file.path(OUT_DIR, "dy180_npp_profiles_per_station_with_ipar.png"), "\n")
  } else {
    cat("  -> no finite cbpm_ipar profiles; skipping ipar overlay plot\n")
  }

  # CTD-only NPP comparison with bbp sensitivity --------------------------------
  # PE curve (depth-resolved), incubation C-uptake (5 depths), CbPM-on-CTD,
  # and CbPM-on-CTD with bbp divided by 3 and by 10 (sensitivity to scattering).
  if (nrow(ctd_results_df) > 0) {
    ctd_pp_profiles_full  <- unnest_ctd(ctd_results_df, "pp_z",       "pp")
    ctd_pp_profiles_bbp3  <- unnest_ctd(ctd_results_df, "pp_z_bbp3",  "pp") |>
      dplyr::filter(is.finite(pp))
    ctd_pp_profiles_bbp10 <- unnest_ctd(ctd_results_df, "pp_z_bbp10", "pp") |>
      dplyr::filter(is.finite(pp))

    ctd_npp_plot <- ggplot() +
      geom_path(data = pe_profiles,
                aes(value, depth_m, group = block_id,
                    colour = "PE curve (CTD)")) +
      geom_point(data = inc_points,
                 aes(value, depth_m, group = block_id,
                     colour = "C14 incubation (CTD)"),
                 shape = 17, size = 2) +
      geom_path(data = ctd_pp_profiles_full,
                aes(pp, depth_m, colour = "CbPM on CTD (bbp)"),
                linewidth = 0.8) +
      geom_path(data = ctd_pp_profiles_bbp3,
                aes(pp, depth_m, colour = "CbPM on CTD (bbp / 3)"),
                linewidth = 0.8, linetype = "dashed") +
      geom_path(data = ctd_pp_profiles_bbp10,
                aes(pp, depth_m, colour = "CbPM on CTD (bbp / 10)"),
                linewidth = 0.8, linetype = "dotted") +
      scale_colour_manual(
        name = "",
        values = c("PE curve (CTD)"         = "black",
                   "C14 incubation (CTD)"   = "black",
                   "CbPM on CTD (bbp)"      = "#1f77b4",
                   "CbPM on CTD (bbp / 3)"  = "#d62728",
                   "CbPM on CTD (bbp / 10)" = "#ff7f0e")) +
      scale_y_reverse(limits = c(150, 0)) +
      facet_wrap(~ stn_label, scales = "free_x") +
      labs(x = "NPP / C-uptake (mg C m^-3 d^-1)",
           y = "Depth (m)",
           title = "DY180 CTD-only NPP profiles — CbPM bbp sensitivity",
           subtitle = "Solid CbPM uses bbp as-is; dashed = bbp / 3; dotted = bbp / 10 (Cphyto = 12128 * bbp470 + 0.59)") +
      theme_bw()
    ggsave(file.path(OUT_DIR, "dy180_ctd_npp_profiles_bbp_sensitivity.png"),
           plot = ctd_npp_plot, width = 11, height = 8, dpi = 150)
    cat("  -> wrote",
        file.path(OUT_DIR, "dy180_ctd_npp_profiles_bbp_sensitivity.png"), "\n")

    # Integrated comparison: PE-curve integral vs CbPM integrals (all bbp variants).
    # PE is integrated over its full depth range; CbPM is integrated 0..z_eu.
    # ctd_pe_int may have multiple blocks per ctd_id (e.g. CTD036) -> average
    # to one PE-integral per station for the scatter.
    ctd_pe_int_per_stn <- ctd_pe_int |>
      group_by(ctd_id) |>
      summarise(ctd_pe_int = mean(ctd_pe_int, na.rm = TRUE), .groups = "drop")

    integ_compare <- ctd_results_df |>
      select(ctd_id, stn_label, rep_label,
             cbpm_int_ctd, cbpm_int_ctd_bbp3, cbpm_int_ctd_bbp10) |>
      left_join(ctd_pe_int_per_stn, by = "ctd_id") |>
      tidyr::pivot_longer(c(cbpm_int_ctd, cbpm_int_ctd_bbp3, cbpm_int_ctd_bbp10),
                          names_to = "method", values_to = "cbpm_int") |>
      mutate(method = recode(method,
                             cbpm_int_ctd       = "CbPM (bbp)",
                             cbpm_int_ctd_bbp3  = "CbPM (bbp / 3)",
                             cbpm_int_ctd_bbp10 = "CbPM (bbp / 10)"),
             method = factor(method,
                             levels = c("CbPM (bbp)",
                                        "CbPM (bbp / 3)",
                                        "CbPM (bbp / 10)")))

    rng_c <- range(c(integ_compare$ctd_pe_int, integ_compare$cbpm_int),
                   na.rm = TRUE)
    ctd_integ_plot <- ggplot(integ_compare,
                             aes(ctd_pe_int, cbpm_int, colour = method)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "grey50") +
      geom_point(size = 3, alpha = 0.85) +
      geom_text(aes(label = rep_label), nudge_y = 60, size = 2.8,
                show.legend = FALSE) +
      scale_colour_manual(values = c("CbPM (bbp)"      = "#1f77b4",
                                     "CbPM (bbp / 3)"  = "#d62728",
                                     "CbPM (bbp / 10)" = "#ff7f0e")) +
      coord_equal(xlim = rng_c, ylim = rng_c) +
      labs(x = "CTD PE-curve depth-integrated NPP (mg C m^-2 d^-1)",
           y = "CbPM-on-CTD depth-integrated NPP (mg C m^-2 d^-1)",
           colour = "",
           title = "DY180 CTD: integrated PE-curve NPP vs CbPM-on-CTD",
           subtitle = "Sensitivity of CbPM to scattering: full bbp vs bbp / 3 vs bbp / 10") +
      theme_bw()
    ggsave(file.path(OUT_DIR, "dy180_ctd_npp_integrated_bbp_sensitivity.png"),
           plot = ctd_integ_plot, width = 7.5, height = 6.5, dpi = 150)
    cat("  -> wrote",
        file.path(OUT_DIR, "dy180_ctd_npp_integrated_bbp_sensitivity.png"), "\n")

    # Baseline CTD-only profile (no bbp sensitivity overlays): PE curve, C14
    # incubation, original CbPM on CTD only.
    ctd_baseline_plot <- ggplot() +
      geom_path(data = pe_profiles,
                aes(value, depth_m, group = block_id,
                    colour = "PE curve (CTD)")) +
      geom_point(data = inc_points,
                 aes(value, depth_m, group = block_id,
                     colour = "C14 incubation (CTD)"),
                 shape = 17, size = 2) +
      geom_path(data = ctd_pp_profiles_full,
                aes(pp, depth_m, colour = "CbPM on CTD"),
                linewidth = 0.8) +
      scale_colour_manual(
        name = "",
        values = c("PE curve (CTD)"       = "black",
                   "C14 incubation (CTD)" = "black",
                   "CbPM on CTD"          = "#1f77b4")) +
      scale_y_reverse(limits = c(150, 0)) +
      facet_wrap(~ stn_label) +
      labs(x = "NPP / C-uptake (mg C m^-3 d^-1)",
           y = "Depth (m)",
           title = "DY180 CTD-only NPP profiles") +
      theme_bw()
    ggsave(file.path(OUT_DIR, "dy180_ctd_npp_profiles_baseline.png"),
           plot = ctd_baseline_plot, width = 11, height = 8, dpi = 150)
    cat("  -> wrote",
        file.path(OUT_DIR, "dy180_ctd_npp_ profiles_baseline.png"), "\n")
  }

  # Chl per-station: float-derived chl + CTD calibrated chl
  chl_plot <- ggplot() +
    geom_path(data = chl_profiles,
              aes(chl, depth_m, colour = factor(float_wmo))) +
    scale_y_reverse(limits = c(150, 0)) +
    facet_wrap(~ stn_label, scales = "free_x") +
    labs(x = "Chlorophyll-a (mg m^-3)", y = "Depth (m)",
         colour = "Float WMO",
         title = "DY180 per-station chlorophyll-a profiles",
         subtitle = "Black solid = calibrated CTD CHL; coloured = float CHLA_ADJUSTED (despiked)") +
    theme_bw()
  if (!is.null(ctd_chl_profiles)) {
    chl_plot <- chl_plot +
      geom_path(data = ctd_chl_profiles,
                aes(chl, depth_m), colour = "black", linewidth = 0.7)
  }
  ggsave(file.path(OUT_DIR, "dy180_chla_profiles_per_station.png"),
         plot = chl_plot, width = 11, height = 8, dpi = 150)
  cat("  -> wrote", file.path(OUT_DIR, "dy180_chla_profiles_per_station.png"), "\n")

  # bbp per-station: float-derived bbp + CTD calibrated bbp
  bbp_plot <- ggplot() +
    geom_path(data = bbp_profiles,
              aes(bbp, depth_m, colour = factor(float_wmo))) +
    scale_y_reverse(limits = c(150, 0)) +
    facet_wrap(~ stn_label, scales = "free_x") +
    labs(x = "Particulate backscatter b_bp(700) (m^-1)", y = "Depth (m)",
         colour = "Float WMO",
         title = "DY180 per-station b_bp(700) profiles",
         subtitle = "Black solid = calibrated CTD bbp; coloured = float BBP700_ADJUSTED (despiked)") +
    theme_bw()
  if (!is.null(ctd_bbp_profiles)) {
    bbp_plot <- bbp_plot +
      geom_path(data = ctd_bbp_profiles,
                aes(bbp, depth_m), colour = "black", linewidth = 0.7)
  }
  ggsave(file.path(OUT_DIR, "dy180_bbp_profiles_per_station.png"),
         plot = bbp_plot, width = 11, height = 8, dpi = 150)
  cat("  -> wrote", file.path(OUT_DIR, "dy180_bbp_profiles_per_station.png"), "\n")
}

cat("\n=== DONE ===\n")

