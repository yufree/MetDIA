#' find peak index of matched precursor ions with given ppm
#' the precursor m/z below 400 is considered as 400 due to heavy mass shift in low mass range
#' @param peaktable provide the precursor ion m/z 
#' @param premz.lib the precursor ion m/z in library
#' @param ms1ppm the tolerance precursor ion match
FindPeaksIdx <- function(peaktable, premz.lib, ms1ppm){
  dev <- ms1ppm * 1e-6
  featureidx <- c(peaktable[, "mz"] >= premz.lib - max(dev * 400, dev * premz.lib) / 2 & peaktable[, "mz"] <= premz.lib + max(dev * 400,dev * premz.lib) / 2)  # feature
}