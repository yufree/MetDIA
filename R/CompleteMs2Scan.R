#' Complete MS2 scans
#' MS2 scans could be lost, find missing MS2 scan and complete with 0
#' @param scantime.ms1 MS1 scan time
#' @param scantime.ms2 MS2 scan time
#' @param scan.ms2 MS2 scans 

CompleteMs2Scan <- function(scantime.ms1, scantime.ms2, scan.ms2){
  
  ms2.intervial <- median(diff(scantime.ms2))
  ms2.fold <- ceiling(length(scantime.ms2) / length(scantime.ms1))
  
  scantime.ms2.theoritical.list <- lapply(1:length(scantime.ms1), function(y){
    scantime.ms12mix.1group <- seq(from = scantime.ms1[y], length.out = ms2.fold + 1, by = ms2.intervial)
    scantime.ms2.1group <- scantime.ms12mix.1group[-1]
  })
  scantime.ms2.theoritical <- do.call(c, scantime.ms2.theoritical.list)
  
  names(scantime.ms2.theoritical) <- c(1 : length(scantime.ms2.theoritical))
  idx.ms2.theoritical <- sapply(scantime.ms2, function(y){
    tmp <- y - scantime.ms2.theoritical
    which.min(tmp[tmp >= 0])
  })
  
  names(scan.ms2) <- names(scantime.ms2)[idx.ms2.theoritical]
  scan.ms2.0 <- list()
  scan.ms2.0 <- sapply(1:length(scantime.ms2.theoritical), function(y){
    scan.ms2.0 <- list(matrix(0, ncol = 2))
  })
  
  scan.ms2.0[idx.ms2.theoritical] <- scan.ms2
  return(scan.ms2.0)
}
