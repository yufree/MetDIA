#' RefineReuslt.ms2: Refine output result of MS2 area
#' RefineReuslt.ms2: Add detailed information, such as KEGG, HMDB, to the result
#' @param result.raw raw result of MS2 area
#' @param inf.detail detail inforamtion containing metabolite name, lab ID, KEGG ID, HMDB ID
#' @param name.sample sample name

RefineReuslt.ms2 <- function(result.raw, inf.detail, name.sample){
  
  id.idx <- ! is.na(result.raw[, "id.add"])
  r.id <- result.raw[id.idx, , drop = FALSE]  # result.id
  r.id.order <- r.id[order(r.id[, "MCI"], decreasing = TRUE), , drop = FALSE]
  r.id.unique <- r.id.order[! duplicated(r.id.order[, "id.add"]), , drop = FALSE]
  
  # direction <- ifelse(is.na(r.id.unique[, "tstat"]), NA, ifelse(r.id.unique[, "tstat"] > 0, "down", "up"))
  id.id <- sapply(r.id.unique[, "id.add"], function(y) strsplit(y, split = "_")[[1]][1])  # identified id
  idx.metname.id <- match(id.id, inf.detail[, "LABID"])
  inf.detail.id <- inf.detail[idx.metname.id, c("name", "KEGG", "HMDB"), drop = FALSE]
  
  idx.area.begin <- which(colnames(result.raw) == name.sample[1])[2]
  idx.area.end <- which(colnames(result.raw) == name.sample[length(name.sample)])[2]
  
  r.refine <- data.frame(Metname = inf.detail.id[, "name"], 
                         MCI = r.id.unique[, "MCI"],
                         KEGG = inf.detail.id[, "KEGG"], 
                         HMDB = inf.detail.id[, "HMDB"], 
                         Labid = id.id, 
                         Name = r.id.unique[, "name"], 
                         mz = r.id.unique[, "mzmed"],
                         #Foldchange = r.id.unique[, "fold"], 
                         #tstat = r.id.unique[, "tstat"], 
                         #direction = direction, 
                         #pvaule = r.id.unique[, "pvalue"], 
                         RTmed = r.id.unique[, "rtmed"],
                         RTmin = r.id.unique[, "rtmin"],
                         RTmax = r.id.unique[, "rtmax"], 
                         r.id.unique[, idx.area.begin:idx.area.end])
  row.names(r.refine) <- NULL
  r.refine
}