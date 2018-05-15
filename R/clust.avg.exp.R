#' Scale data
#'
#' This function takes an object of class scSeqR and scales the normalized data.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' scale.data(my.obj)
#' }
#' @export
clust.avg.exp <- function (x = NULL,
                           clust.type = "tsne",
                           clust.dim = 2) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clust.type == "tsne") {
    if (clust.dim == 2) {
      DATA <- x@tsne.data
    }
    if (clust.dim == 3) {
      DATA <- x@tsne.data.3d
    }
  }
  if (clust.type == "pca") {
    if (clust.dim == 2) {
      DATA <- x@pca.data
    }
    if (clust.dim == 3) {
      DATA <- x@pca.data.3d
    }
  }
  sampleCondition <- DATA$clusters
  conditions <- unique(sampleCondition)
  DATA1 <- DATA
  Table = data.matrix(x@main.data)
  for(i in conditions){
    DATA <- Table[,row.names(subset(DATA1, sampleCondition == i))]
    DATA <- apply(DATA, 1, function(DATA) {mean(DATA)})
    DATA <- as.data.frame(DATA)
    Name=paste("meanExp_cluster",i,".txt",sep="_")
    NameCol=paste("cluster",i,sep="_")
    colnames(DATA) <- NameCol
    DATA <- cbind(gene = rownames(DATA), DATA)
    head(DATA)
    write.table((DATA),file=Name,sep="\t", row.names =F)
  }
  multmerge = function(mypath){
    filenames=list.files(pattern="meanExp")
    datalist = lapply(filenames, function(x){read.table(file=x,header=T)})
    Reduce(function(x,y) {merge(x,y)}, datalist)
  }
  MeanExpForClusters <- multmerge()
  file.remove(list.files(pattern="meanExp"))
  attributes(x)$clust.avg <- MeanExpForClusters
  return(x)
}
