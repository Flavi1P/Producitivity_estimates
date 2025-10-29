#' Despike a time series
#'
#' Compute a smooth baseline and identify anomalous spikes
#' following Briggs et al. (2011).
#'
#' @param var Numeric vector or ts
#' @param window_size Integer, rolling window length
#' @param spike_method "median" (default) or "minmax"
#'
#' @return A list with:
#' \describe{
#'   \item{baseline}{Smoothed baseline vector}
#'   \item{spikes}{Residuals = var - baseline}
#' }
#'
despike <- function(var, window_size, spike_method = "median") {
  library(zoo)
  
  arr <- as.numeric(var)
  n <- length(arr)
  
  baseline <- rep(NA_real_, n)
  mask <- !is.na(arr)
  
  if (startsWith(spike_method, "min")) {
    base_min <- rollapply(arr[mask], width = window_size, FUN = min, 
                          fill = NA, align = "center", na.rm = TRUE)
    base <- rollapply(base_min, width = window_size, FUN = max, 
                      fill = NA, align = "center", na.rm = TRUE)
  } else {
    base <- rollapply(arr[mask], width = window_size, FUN = median, 
                      fill = NA, align = "center", na.rm = TRUE)
  }
  
  baseline[mask] <- base
  spikes <- arr - baseline
  
  # Placeholder for NetCDF attribute transfer
  # baseline <- transfer_nc_attrs(var, baseline, suffix = "_baseline")
  # spikes   <- transfer_nc_attrs(var, spikes,   suffix = "_spikes")
  
  return(list(
    baseline = baseline,
    spikes = spikes
  ))
}
