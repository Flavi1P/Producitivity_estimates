#' Daylength in hours (Brock model, Forsythe et al. 1995)
#'
#' Computes daylight duration from latitude and calendar date.
#'
#' @param year Integer, e.g. 2025
#' @param month Integer (1–12)
#' @param day Integer (1–31)
#' @param lat Latitude in degrees (north positive, south negative)
#'
#' @return Daylength in hours
#'
daylength <- function(year, month, day, lat) {
  library(lubridate)
  
  # Day of year
  date <- ymd(sprintf("%04d-%02d-%02d", year, month, day))
  dayOfYear <- yday(date)
  
  # Latitude in radians
  latRad <- lat * pi / 180
  
  # Declination of Earth (degrees)
  decl <- 23.45 * sin((2 * pi / 365) * (283 + dayOfYear))
  
  # Convert declination to radians
  declRad <- decl * pi / 180
  
  # Edge cases: polar day / night
  x <- -tan(latRad) * tan(declRad)
  if (x <= -1) {
    return(24.0)
  } else if (x >= 1) {
    return(0.0)
  } else {
    # Hour angle in degrees
    hourAngle <- acos(x) * 180 / pi
    return(2.0 * hourAngle / 15.0)  # convert degrees to hours
  }
}
