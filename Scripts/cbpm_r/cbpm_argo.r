#' CbPM-Argo: Depth-resolved phytoplankton growth and NPP
#'
#' Adapted from Westberry et al. (2008), following Lionel A. Arteaga’s Python implementation.
#' Requires AustinPetzold_1986() and daylength() equivalents in R.
#'
#' @param chl_z Chlorophyll profile (mg m^-3), length 200
#' @param Cphyto_z Phytoplankton carbon profile (mg m^-3), length 200
#' @param irr Surface irradiance (E m^-2 d^-1)
#' @param year Four digit year
#' @param month Month number (1–12)
#' @param day Day of month
#' @param lat Latitude (decimal degrees)
#'
#' @return A list with:
#' \describe{
#'   \item{pp_z}{Primary production profile (mg C m^-2 d^-1)}
#'   \item{mu_z}{Phytoplankton growth rate (d^-1)}
#'   \item{par_z}{PAR profile}
#'   \item{prcnt_z}{Fraction of surface irradiance}
#'   \item{nutTempFunc_z}{Nutrient limitation index}
#'   \item{IgFunc_z}{Light limitation index}
#'   \item{mzeu}{1% light depth (m)}
#' }
#'
cbpm_argo <- function(chl_z, Cphyto_z, irr, year, month, day, lat) {
  source("Scripts/cbpm_r/austinpetzold.r")
  source("Scripts/cbpm_r/daylength.r")
  library(matrixStats)
  library(purrr)
  
  # Spectral variables
  Lambda <- c(400, 412, 443, 490, 510, 555, 625, 670, 700)
  parFraction <- c(0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024)
  X <- c(0.11748, 0.122858, 0.107212, 0.07242, 0.05943, 0.03996, 0.04000, 0.05150, 0.03000)
  e <- c(0.64358, 0.653270, 0.673358, 0.68955, 0.68567, 0.64204, 0.64700, 0.69500, 0.60000)
  Kw <- c(0.01042, 0.007932, 0.009480, 0.01660, 0.03385, 0.06053, 0.28400, 0.43946, 0.62438)
  
  # Model constants
  y0 <- 0.0003
  umax <- 2.0
  
  ndepth <- length(chl_z)
  
  chlC_z <- rep(NA_real_, ndepth)
  nutTempFunc_z <- rep(NA_real_, ndepth)
  IgFunc_z <- rep(NA_real_, ndepth)
  mu_z <- rep(NA_real_, ndepth)
  prcnt_z <- rep(NA_real_, ndepth)
  pp_z <- rep(NA_real_, ndepth)
  Ezlambda <- matrix(NA_real_, nrow = ndepth, ncol = length(Lambda))
  par_z <- rep(NA_real_, ndepth)
  mzeu <- NA_real_
  
  # Attenuation coefficient at 490nm
  chl_surf <- chl_z[1]
  k490 <- 0.0166 + 0.0773 * chl_surf^0.6715
  
  # Daylength (use insol::daylength)
  Daylength <- daylength(year, month, day, lat)
  
  klambda <- map_dbl(Lambda, ~ austinpetzold(.x, k490))
  E0 <- irr * parFraction
  
  # Loop over depth
  for (z in seq_len(ndepth)) {
    if (z == 1) {
      Ezlambda[z, ] <- E0 * 0.975 * exp(-klambda * (z - 1))
    } else {
      for (n in seq_along(Lambda)) {
        chl_prev <- max(chl_z[z - 1], 0)
        kbio <- X[n] * chl_prev^e[n]
        kd <- Kw[n] + kbio
        Ezlambda[z, n] <- Ezlambda[z - 1, n] * exp(-kd * 1)
      }
    }
    
    # PAR as trapezoidal integration
    par_z[z] <- sum((diff(Lambda) * (Ezlambda[z, -1] + Ezlambda[z, -length(Lambda)])) / 2)
    
    chlC_z[z] <- chl_z[z] / Cphyto_z[z]
    chlCarbonMax_z <- 0.022 + (0.045 - 0.022) * exp(-3.0 * par_z[z] / Daylength)
    
    nutTempFunc_z[z] <- (chlC_z[z] - y0) / (chlCarbonMax_z - y0)
    nutTempFunc_z[z] <- min(nutTempFunc_z[z], 1)
    
    IgFunc_z[z] <- 1 - exp(-5.0 * par_z[z] / Daylength)
    IgFunc_z[z] <- min(IgFunc_z[z], 1)
    
    mu_z[z] <- umax * nutTempFunc_z[z] * IgFunc_z[z]
    mu_z[z] <- min(mu_z[z], umax)
    
    prcnt_z[z] <- par_z[z] / (irr * 0.975)
    
    if (!is.na(prcnt_z[z]) && prcnt_z[z] >= 0.01) {
      mzeu <- z
    }
    
    pp_z[z] <- mu_z[z] * Cphyto_z[z]
  }
  
  return(list(
    pp_z = pp_z,
    mu_z = mu_z,
    par_z = par_z,
    prcnt_z = prcnt_z,
    nutTempFunc_z = nutTempFunc_z,
    IgFunc_z = IgFunc_z,
    mzeu = mzeu
  ))
}
