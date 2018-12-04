#' Complete MS1 scans
#' MS1 scans could be lost, find missing MS1 scan and complete with 0
#' @param scantime.ms1.raw raw MS1 scan time
#' @param scan.ms1.raw raw MS1 scans, could contain missing scans

CompleteMs1Scan <- function(scantime.ms1.raw, scan.ms1.raw){
  # theoritical scan inf
  num.ms1scan <- sapply(scantime.ms1.raw, length)
  idx.maxscan <- which.max(num.ms1scan)
  scantime.ms1.full <- scantime.ms1.raw[[idx.maxscan]]
  scantime.ms1.theoritical <- scantime.ms1.full
  names(scantime.ms1.theoritical) <- seq(scantime.ms1.theoritical)

  lapply(seq(scan.ms1.raw), function(idx.s){
    if(length(scan.ms1.raw[[idx.s]]) == max(num.ms1scan)){
      scan.ms1.raw[[idx.s]]
    }else{
      scantime.ms1.missing <- scantime.ms1.raw[[idx.s]]
      scan.ms1.raw.s <- scan.ms1.raw[[idx.s]]

      idx.ms1.theoritical <- sapply(scantime.ms1.missing, function(y) {
        tmp <- y - scantime.ms1.theoritical
        which.min(abs(tmp))
      })
      names(scan.ms1.raw.s) <- names(scantime.ms1.theoritical)[idx.ms1.theoritical]

      scan.ms1.0 <- list()
      scan.ms1.0 <- sapply(seq(scantime.ms1.theoritical),function(y) {
        scan.ms2.0 <- list(matrix(0, ncol = 2))
      })
      scan.ms1.0[idx.ms1.theoritical] <- scan.ms1.raw.s
      return(scan.ms1.0)
    }
  })
}
