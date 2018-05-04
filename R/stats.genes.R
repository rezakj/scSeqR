#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class scSeqR and creats plot.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' plot.stats(my.obj)
#' }
#' @export
stats.genes <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@main.data
  NumbCellPerGene <- as.data.frame(rowSums(DATA > 0))
  MeanExp <- as.data.frame(rowMeans(DATA))
  dat=as.matrix(DATA)
  mySD <- as.data.frame(apply(dat, 1, function(dat) sd(dat)))
  colnames(mySD) <- "SD"
  colnames(MeanExp) <- "MeanExp"
  colnames(NumbCellPerGene) <- "numberOfCells"
  NumbCellPerGene <- cbind(gene = rownames(NumbCellPerGene),NumbCellPerGene, MeanExp = MeanExp, SD = mySD)
  Name=paste("geneStatFile_",dim(DATA)[1],"_genes_",dim(DATA)[2],"_cells.tsv", sep="")
  write.table((NumbCellPerGene),file=Name,sep="\t",row.names =F)
  return(print("geneStatFile is created"))
}
