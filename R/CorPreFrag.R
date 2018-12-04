#'  Correlation between precursor and fragment ion
#'  Calculate correlation between precursor and fragment ions with Pearsom correlation
#'  @param premzi a list contains a dataframe of precursor ion m/z and intensity
#'  @param fragmzi a list contains a dataframe of precursor ion m/z and intensity 


CorPreFrag <- function(premzi, fragmzi){
  result.cor <- list()
  for(i in 1:length(premzi)){
    result.cor.2 <- c()  
    premzi.single <- premzi[[i]][, 2, drop = FALSE]
    for(j in 1:length(fragmzi[[i]])){
      fragmzi.single <- fragmzi[[i]][[j]][, 2, drop = FALSE]
      if(length(premzi.single) < 3 ){
        tmp.cor <- NA
        result.cor.2 <- c(result.cor.2, tmp.cor)
      }else{
        tmp.cor <- cor.test(premzi.single, fragmzi.single)  
        result.cor.2 <- c(result.cor.2, tmp.cor $ estimate)
      }
    }
    result.cor[[i]] <- result.cor.2
  }
  return(result.cor)
}