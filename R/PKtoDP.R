#' Extract peak information and calculate score
#' Use precursor ion m/z in library to match the detected peaks. Then targeted extract m/z and intensity of precursor and fragment ions. Calculate dot product and correlation with these data
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

PKtoDP <- function(peaktable, peaktable.corrected, scantime, speclib.single, scan.ms1, scan.ms2, ms1ppm, ms2ppm, massrange.ms1, peakgroup.num, mcicutoff = cutoff, windows = file.windows){
  
  # find peak index of matched precursor ions with given ppm
  ms1.mzrt.idx <- FindPeaksIdx(peaktable, premz.lib=speclib.single[1, 1], ms1ppm=ms1ppm)
  ms1.mzrt <- peaktable[ms1.mzrt.idx, , drop = FALSE]
  ms1.mzrt.corrected <- peaktable.corrected[ms1.mzrt.idx, , drop = FALSE]
  
  # extract scan index with rt of matched ms1 feature
  scan.idx <- ScanIdx(scantime, ms1.mzrt)
  if(all(is.na(scan.idx))){
    return(NA)
  }else{
    
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
  
  # calculate correlation between the intensity of precursor ion and fragment ion
  cor.result <- CorPreFrag(pre.mzi, frag.mzi)
  
  # use top n/2 correlation to calculate cor score
  cor.mean <- sapply(cor.result, function(cor.single) {
    cor.single.order <- cor.single[order(cor.single, decreasing = TRUE)]
    cor.single.order.nato0 <- ifelse(is.na(cor.single.order), 0, cor.single.order)  # replace na with 0
    round(mean(cor.single.order.nato0[1 : ceiling(length(cor.single) / 2)]), 3) 
  })
  # use all correlation to calculate cor score
  cor.mean.all <- sapply(cor.result, function(cor.single) {
    cor.single.order <- cor.single[order(cor.single, decreasing = TRUE)]
    cor.single.order.nato0 <- ifelse(is.na(cor.single.order), 0, cor.single.order)  # replace na with 0
    round(mean(cor.single.order.nato0), 3) 
  })
  
  # Plot12(pre.mzi, frag.mzi, speclib.single, cor.result, peakgroup.num)
  # Plot12.cutoff(pre.mzi, frag.mzi, speclib.single, cor.result, peakgroup.num)
  
  # Extract intensity at summit
  ms2.maxi <- Ms2MaxI(pre.mzi, frag.mzi)
  # remove fragment ions with correlation below 0.6
  ms2.maxi.cutoff <- Ms2MaxIFilter(ms2.maxi, cor.result)
  # PlotDP(ms2.maxi, ms2.maxi.cutoff, speclib.single, peakgroup.num) 
  
  # calculate the dot product between library spectrum and filtered experiment spectrum
  dp.result <- sapply(1:length(ms2.maxi.cutoff), function(y) {
    score <- (sum(ms2.maxi.cutoff[[y]] * speclib.single[, "LibraryIntensity"]) / sqrt(sum(ms2.maxi.cutoff[[y]] ^ 2) * sum(speclib.single[, "LibraryIntensity"] ^ 2)))
    round(score, 3)
  })
  
  # calculate the dot product between library spectrum and filtered experiment spectrum
  dp.result.c <- sapply(1:length(ms2.maxi), function(y) {
    score <- (sum(ms2.maxi[[y]] * speclib.single[, "LibraryIntensity"]) / sqrt(sum(ms2.maxi[[y]] ^ 2) * sum(speclib.single[, "LibraryIntensity"] ^ 2)))
    round(score, 3)
  })
  
  dp.result0 <- ifelse(is.na(dp.result), 0, dp.result)  # change NA to 0
  dp.result.c0 <- ifelse(is.na(dp.result.c), 0, dp.result.c)
  
  # plot raw spectra
  PlotAllRaw(ms1 = pre.mzi.raw, ms2 = frag.mzi.raw, speclib.single, corvalue = cor.result, ms2.maxi = ms2.maxi, ms2.maxi.cutoff = ms2.maxi.cutoff, dpvalue = dp.result.c0, peakgroup.num)
  PlotAllRawFilter(ms1 = pre.mzi.raw, ms2 = frag.mzi.raw, speclib.single, corvalue = cor.result, ms2.maxi = ms2.maxi, ms2.maxi.cutoff = ms2.maxi.cutoff, dpvalue = dp.result.c0, corvalue.mean = cor.mean, peakgroup.num = peakgroup.num, mcicutoff = mcicutoff)
  
  # score.combine <- (dp.result0 + cor.mean) / 2
  # score.fall <- (dp.result0 + cor.mean.all) / 2  # !!!!!!!!!!!!!!
  score.cpart <- round((dp.result.c0 + cor.mean) / 2, 3)  # 
  # score.call <- (dp.result.c0 + cor.mean.all) / 2  # !!!!!!!!!!!
  
  # dp.idx <- (! is.na(dp.result)) & (dp.result >= 0.8)
  # information.table <- cbind(dp.result[dp.idx], ms1.mzrt[dp.idx, , drop = FALSE])
  # return(information.table)
  
  information.table <- cbind(score.cpart, dp.result.c0, cor.mean, ms1.mzrt.corrected)
  colnames(information.table) <- c("score.cpart", "dpc", "cormean.part", colnames(ms1.mzrt.corrected)) 
  
  return(information.table)
}