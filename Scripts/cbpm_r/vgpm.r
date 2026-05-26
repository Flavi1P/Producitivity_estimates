#' VGPM (Behrenfeld & Falkowski, 1997) — Vertically Generalized Production Model
#'
#' Port of `opp_befa()` from Scripts/3d_products.ipynb (cell 43).
#' Returns water-column-integrated NPP (mg C m^-2 d^-1).
#'
#' @param chl Surface chlorophyll-a concentration (mg m^-3)
#' @param irr Surface PAR (E m^-2 d^-1)
#' @param sst Sea surface temperature (degC)
#' @param dayL Day length (hours)
#'
#' @return Named list:
#'   - npp:     integrated NPP (mg C m^-2 d^-1)
#'   - pb_opt:  temperature-dependent maximum carbon fixation rate
#'   - z_eu:    euphotic depth (m, 1% surface light)
#'   - irrFunc: light limitation factor
opp_befa <- function(chl, irr, sst, dayL) {
  # Total chl in the water column (Morel & Berthon 1989)
  chl_tot <- if (chl < 1.0) 38.0 * chl ^ 0.425 else 40.2 * chl ^ 0.507

  # Euphotic depth (Morel 1988)
  z_eu <- 200.0 * chl_tot ^ (-0.293)
  if (z_eu <= 102.0) z_eu <- 568.2 * chl_tot ^ (-0.746)

  # Pb_opt — temperature-dependent maximum NPP (Behrenfeld & Falkowski 1997 7th-order polynomial)
  pb_opt <- if (sst < -10.0) {
    0.0
  } else if (sst < -1.0) {
    1.13
  } else if (sst > 28.5) {
    4.0
  } else {
    1.2956 + 2.749e-1 * sst + 6.17e-2 * sst^2 -
      2.05e-2 * sst^3 + 2.462e-3 * sst^4 -
      1.348e-4 * sst^5 + 3.4132e-6 * sst^6 -
      3.27e-8 * sst^7
  }

  irrFunc <- 0.66125 * irr / (irr + 4.1)
  npp <- pb_opt * chl * dayL * irrFunc * z_eu

  list(npp = npp, pb_opt = pb_opt, z_eu = z_eu, irrFunc = irrFunc)
}
