#' Extract extract m/z and intensity from scan 
#' applied on the extraction of both precursor ions and fragment ions
#' @param scan.all all scans of ms1 or ms2
#' @param scan.index scan index selected previously
#' @param mz.lib mz of targeted ions in library
#' @param ppm the given tolerance

ExtractMzI <- function(scan.all, scan.index, mz.lib, ppm){
  sample.num <- as.numeric(names(scan.index)[1])
  dev <- ppm * 1e-6
  record.all <- list()
  j <- 1
  for(i in scan.index){
    scan.single <- scan.all[[sample.num]][[i]]
    record.all[[j]] <- scan.single[(scan.single[, 1] <= mz.lib[1] + max(dev * 400 / 2, dev * mz.lib[1] / 2)) 
                                   & (scan.single[, 1] >= mz.lib[1] - max(dev * 400 / 2, dev * mz.lib[1] / 2)), , drop = FALSE]
    j <- j + 1
  }
  return(record.all)
}