#'  Mirror plot and peak group plot with raw intensity data for all peaks
#'  Use m/z ang intensity of precursor and fragment ions for mirror plot and peak group plot. The left is peak group plot with all ions, the black line stands for precursor ion, and the colorful lines stand for fragment ions. The right is the mirror plot of all fragmet ions.
#'  @param ms1 precursor ion m/z and raw intensity
#'  @param ms2 fragment ion m/z and raw intensity
#'  @param feature spectrum information in library
#'  @param corvalue correlation value
#'  @param ms2.maxi MS2 intensity at summit
#'  @param ms2.maxi.cutoff MS2 intensity filtered by correlation (0.6) at summit 
#'  @param dpvalue dot product value
#'  @param peakgroup.num differentiate from different samples

PlotAllRaw <- function(ms1, ms2, feature, corvalue, ms2.maxi, ms2.maxi.cutoff, dpvalue, peakgroup.num){
  dir.create("EIC and spectra") 
  setwd("EIC and spectra")
  for(i in 1:length(ms1)){
    filename <- paste(feature[1, "Name"], peakgroup.num, i, ".png", sep = " ")
    png(file = filename, width = 960, height = 480)
    par(mfrow=c(1,2))
    plot(ms1[[i]][, 2], ylim = c(0, max(ms1[[i]][, 2])), type = "b", pch = 19, cex = 1, main = "EIC", ylab = "intensity", xlab = "index")
    for(j in 1:length(ms2[[i]])){
      points(ms2[[i]][[j]][, 2], cex = 0.8, col = j + 1)
      lines(ms2[[i]][[j]][, 2], cex = 0.8, col = j + 1)
      # in case that the max intensity is 0
    }

    quary.before <- cbind(feature[, "ProductMz"], ms2.maxi[[i]] / max(ms2.maxi[[i]]))
    lib.spec <- cbind(feature[, "ProductMz"], - feature[, "LibraryIntensity"] / max(feature[, "LibraryIntensity"]))
    ql.before <- rbind(quary.before, lib.spec)
    plot(ql.before, type = "h", pch = 19, cex = 1, ylim = c(-1, 1), main = "Spectrum to spectrum comparison", xlab = "m/z", ylab = "relative intensity")
    abline(h = 0)

    dev.off()
  }
  setwd("..")
}