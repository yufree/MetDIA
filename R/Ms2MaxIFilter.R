#'  Remove fragment ions with correlation below 0.6
#'  replace intensity of fragment ion whose correlation is below 0.6 with 0
#'  @param ms2maxi a list contains fragment ion intensity at the summit 
#'  @param corresult a list contains correaltion between each fragment ion and precursor ion

Ms2MaxIFilter <- function(ms2maxi, corresult){
  for(i in 1:length(ms2maxi)){
    ms2maxi[[i]][corresult[[i]] < 0.6] <- 0
  }
  return(ms2maxi)
}
