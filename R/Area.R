#' Extract peak information and calculate peak area
#' Use precursor ion m/z in library to match the detected peaks. Then targeted extract m/z and intensity of precursor and fragment ions. Calculate their peak areas
#' @param peaktable peaktable with raw retention time
#' @param peaktable.corrected peaktable with corrected retention time
#' @param scantime MS1 scan time
#' @param speclib.single metabolite information in library, contains PrecursorMZ, ProductMz, LibraryIntensity and Name
#' @param scan.ms1 MS1 scan, contain m/z and intensity
#' @param scan.ms2 MS2 scan, contain m/z and intensity
#' @param ms1ppm tolerance for precursor ion match
#' @param ms2ppm tolerance for precursor ion match
#' @param massrange.ms1 MS1 mass range
#' @param peakgroup.num for plot
#' @param mcicutoff the MCI score cutoff


Area <- function(peaktable, peaktable.corrected, scantime, speclib.single, scan.ms1, scan.ms2, ms1ppm, ms2ppm, massrange.ms1, peakgroup.num, mcicutoff = cutoff, windows = file.windows){
  
  # find peak index of matched precursor ions with given ppm
  ms1.mzrt.idx <- FindPeaksIdx(peaktable, premz.lib=speclib.single[1, 1], ms1ppm=ms1ppm)
  ms1.mzrt <- peaktable[ms1.mzrt.idx, , drop = FALSE]
  ms1.mzrt.corrected <- peaktable.corrected[ms1.mzrt.idx, , drop = FALSE]
  
  # extract scan index with rt of matched ms1 feature
  scan.idx <- ScanIdx(scantime, ms1.mzrt)
  if(all(is.na(scan.idx))){
    return(NA)
  }else{
    # browser("area")
  }
  
  # extract m/z and intensity of precursor ions and fragment ions according to scan indexϢ
  pre.mzi.raw <- lapply(scan.idx, function(y) PreMzI(scan.ms1, y, speclib.single[1, 1], ms1ppm = ms1ppm))
  frag.mzi.raw <- lapply(scan.idx, function(y){
    lapply(speclib.single[1 : nrow(speclib.single), "ProductMz"], function(z) 
      FragMzI(scan.ms2, y, speclib.single[1, 1], z, ms2ppm = ms2ppm, scantime.ms1 = scantime, massrange.ms1 = massrange.ms1, windows = windows))
  })
  
  # predict intensity for smooth
  pre.mzi <- lapply(pre.mzi.raw, function(pre.mzi.single){
    if(nrow(pre.mzi.single) > 3){
      pre.inf <- data.frame(idx = seq(length(pre.mzi.single[, 2])), intensity = pre.mzi.single[, 2])
      pre.fit <- loess(intensity~idx, data=pre.inf, span = 0.3, degree = 1)
      pre.predict <-  predict(pre.fit, data.frame(idx=seq(nrow(pre.mzi.single))))
      cbind(pre.mzi.single[, 1], pre.predict)
    }else{
      pre.mzi.single
    }
  })    
  
  frag.mzi <- lapply(frag.mzi.raw, function(frag.mzi.1st){
    lapply(frag.mzi.1st, function(frag.mzi.2nd){
      if(nrow(frag.mzi.2nd) > 3){
        frag.inf <- data.frame(idx = seq(length(frag.mzi.2nd[, 2])), intensity = frag.mzi.2nd[, 2])
        frag.fit <- loess(intensity~idx, data=frag.inf, span = 0.3, degree = 1)
        frag.predict <-  predict(frag.fit, data.frame(idx=seq(nrow(frag.mzi.2nd))))
        cbind(frag.mzi.2nd[, 1], frag.predict)
      }else{
        frag.mzi.2nd
      }
    })
  })  
  
  # peak area
  pre.area <- lapply(pre.mzi.raw, function(pre.mzi.single){
    sum(pre.mzi.single[, 2])
  })
  frag.area <- lapply(frag.mzi.raw, function(frag.mzi.s1){
    sapply(frag.mzi.s1, function(frag.mzi.s2){
      sum(frag.mzi.s2[, 2])
    })
  })
  
  # calculate correlation between the intensity of precursor ion and fragment ion
  cor.result <- CorPreFrag(pre.mzi, frag.mzi)
  cor.result.no0 <- lapply(cor.result, function(cor.s){
    ifelse(is.na(cor.s), 0, round(cor.s, 2))
  })
  speclib.single
  area.complex <- sapply(seq(length(pre.area)), function(idx){ 
    frag <- paste(round(speclib.single[, "ProductMz"], 5), 
                  round(cor.result.no0[[idx]], 2), 
                  round(frag.area[[idx]], 1), 
                  sep = ",", collapse = ";")
    pre <- paste(round(speclib.single[1, "PrecursorMz"], 5),
                 1, 
                 round(pre.area[[idx]], 1), 
                 sep = ",")
    paste(pre, frag, sep = ";")
  })
  # pre mz, pre cor, pre area;frag1 mz, frag1 cor, frag1 area;...

  
  information.table <- cbind(area = area.complex, ms1.mzrt.corrected)
  
  return(information.table)
}