#' Generate the peakgroup with raw rentiontime
#' the retention time of peakgroup generated from xset3 (filled peaks) is aligned, need the raw retention time for scan extraction
#' 
#' @param xset.fillpeaks the xset with filled gaps
#' @param groupindex.single one sample only contributes one peak to a peak group 

PeakGroupRaw <- function(xset.fillpeaks=xset3, groupindex.single=groupindex.single){  
  peaks <- xset.fillpeaks @ peaks
  peakgroup.corrected <- lapply(groupindex.single, function(y) peaks[y, ])
  
  rt.corrected <- xset.fillpeaks @ rt $ corrected
  rt.raw <- xset.fillpeaks @ rt $ raw
  # peakgroup.raw <- peakgroup.corrected
  peakgroup.raw <- lapply(peakgroup.corrected, function(peaktable) {
    peaktable.raw.t <- sapply(1:nrow(peaktable), function(y) {
      rt.idx <- which.min(abs(peaktable[y, "rt"] - rt.corrected[[peaktable[y, "sample"]]]))
      peaktable[y, "rt"] <- rt.raw[[peaktable[y, "sample"]]][rt.idx]
      rtmin.idx <- which.min(abs(peaktable[y, "rtmin"] - rt.corrected[[peaktable[y, "sample"]]]))
      peaktable[y, "rtmin"] <- rt.raw[[peaktable[y, "sample"]]][rtmin.idx]
      rtmax.idx <- which.min(abs(peaktable[y, "rtmax"] - rt.corrected[[peaktable[y, "sample"]]]))
      peaktable[y, "rtmax"] <- rt.raw[[peaktable[y, "sample"]]][rtmax.idx]
      peaktable[y, ]
    })
    peaktable.raw <- t(peaktable.raw.t)
  })
}
