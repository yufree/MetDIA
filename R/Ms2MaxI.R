#'  Extract fragment ion intensity at summit
#'  Extract fragment ion intensity at the MS1 EIC summit
#'  @param premzi a list contains a dataframe of precursor ion m/z and intensity
#'  @param fragmzi a list contains a dataframe of precursor ion m/z and intensity 

Ms2MaxI <- function(premzi, fragmzi){
  result <- list()
  for(i in 1:length(premzi)){
    idx <- which.max(premzi[[i]][, 2])
    fragmzi.max.i <- sapply(fragmzi[[i]], function(y) y[idx, 2])
    result[[i]] <- fragmzi.max.i
  }
  return(result)
}