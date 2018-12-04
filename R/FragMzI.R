#'  Extract m/z and intensity of fragment ion
#'  Select SWATH window automatically according to precursor ion m/z, the extract fragment ion information 
#'  @param scan.ms2 all ms2 scan
#'  @param scan.index scan index selected previously
#'  @param premz.lib precursor ion mz in library
#'  @param fragmz.lib fragment ion mz in library
#'  @param ms2ppm m/z tolerance of fragment ions
#'  @param scantime.ms1 ms1 scantime
#'  @param massrange.ms1 ms1 mass range

FragMzI <- function(scan.ms2, scan.idx, premz.lib, fragmz.lib, ms2ppm, scantime.ms1, massrange.ms1, windows = file.windows){
  ms1.length <- length(scantime.ms1[[1]])
  idx <- which.max(premz.lib > windows[, "begin"] & premz.lib <= windows[, "end"])
  num.windows <- nrow(windows)
  ms2.windows.idx <- seq(0, ms1.length - 1) *num.windows + idx 
  
#   ms1.length <- length(scantime.ms1[[1]])
#   ms2.length <- length(scan.ms2[[1]])
#   ms2.fold <- ms2.length / ms1.length  # num.window
#   ms2.window <- (massrange.ms1[[1]][2] - massrange.ms1[[1]][1]) / ms2.fold  # window size
#   
#   ms2.window.idx <- seq(0, ms1.length - 1) * ms2.fold + ceiling((premz.lib - massrange.ms1[[1]][1]) / ms2.window)  # window index of all windows containing the precursor ion (the length is the same with the number of ms1 scans)

  ms2.scan.want.idx <- ms2.windows.idx[scan.idx]
  names(ms2.scan.want.idx) <- names(scan.idx)[1]
  record.ms2 <- ExtractMzI(scan.ms2, ms2.scan.want.idx, fragmz.lib, ms2ppm)
  record.ms2.c <- t(sapply(record.ms2, function(y) CompleteI(y, fragmz.lib)))
  return(record.ms2.c)
}
