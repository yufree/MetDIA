#'  Extract m/z and intensity of precursor ion
#'  Extract m/z and intensity of precursor ion and complete missing intensity
#'  @param scan.ms1 all ms1 scan
#'  @param scan.index scan index selected previously
#'  @param premz.lib precursor ion mz in library
#'  @param ms1ppm m/z tolerance of precursor ions

PreMzI <- function(scan.ms1, scan.index, premz.lib, ms1ppm){
  record.ms1 <- ExtractMzI(scan.ms1, scan.index, premz.lib, ms1ppm)
  
  record.ms1.c <- t(sapply(record.ms1, function(y) CompleteI(y, premz.lib)))
  return(record.ms1.c)
}