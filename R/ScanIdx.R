#' extract scan index with rt of matched ms1 feature
#' extract scan index with rt of matched ms1 feature. If without any match, return NA
#' @param scantime all ms1 scantime
#' @param ms1.mzrt mz and rt information of ms1matched ms1 feature 
ScanIdx <- function(scantime, ms1.mzrt){
  if(nrow(ms1.mzrt) == 0) return(NA)
  scan.idx <- lapply(1:nrow(ms1.mzrt), function(rnum){
    idx <- which(scantime[[ms1.mzrt[rnum, "sample"]]] >= ms1.mzrt[rnum, "rtmin"] & 
                   scantime[[ms1.mzrt[rnum, "sample"]]] <= ms1.mzrt[rnum, "rtmax"])
    names(idx) <- ms1.mzrt[rnum, "sample"]  # indicate these scan from which sample
    idx
  })
}  
