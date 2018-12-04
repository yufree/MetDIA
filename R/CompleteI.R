#' complete intensity information
#' the extract intensity may be null or more than one. Need to give mz (mz.lib) and int (0) to NULL. Selet the higher intensity if more than 1
#' @param record a matrix of m/z and intensity
#' @param the mz in library 
CompleteI <- function(record, mz.lib){
  if(nrow(record) == 0){
    record <- matrix(c(mz.lib, 0), 1, 2)
    colnames(record) <- c("mz", "intensity")
  }else{
    if(nrow(record) > 1){
      record <- record[which.max(record[, 2]), , drop = FALSE]
    }
  }
  return(record)
}
