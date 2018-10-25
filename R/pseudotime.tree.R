#' Pseudotime Tree
#'
#' This function takes an object of class scSeqR and marker genes for clusters and performs pseudotime for differentiation or time course analysis.
#' @param x An object of class scSeqR.
#' @param clust.method Choose from "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid", defult is "complete".
#' @param dist.method Choose from "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski", defult is "euclidean".
#' @param clust.names A list of names for clusters.
#' @param marker.genes A list of marker genes for clusters.
#' @param label.offset Space between names and tree, defult = 0.5.
#' @param type Choose from "classic", "unrooted", "fan", "cladogram", "radial", defult = "classic".
#' @param cex Text size, defult = 1.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' my.obj <- run.pca(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
#' }
#' @import ape
#' @export
pseudotime.tree <- function (x = NULL,
                             marker.genes = "NULL",
                             clust.names = "NULL",
                             dist.method = "euclidean",
                             clust.method = "complete",
                             label.offset = 0.5,
                             type = "classic",
                             hang = 1,
                             cex = 1) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  # geth the genes and scale them based on model
  DATA <- x@clust.avg
  row.names(DATA) <- DATA$gene
  DATA <- DATA[,-1]
  if (clust.names[1] != "NULL") {
    colnames(DATA) <- clust.names
  }
  if (marker.genes[1] == "NULL") {
    stop("provide marker genes for clusters (e.g. top 10 for each cluster)")
  }
  MyGenes <- marker.genes
  topGenes <- subset(DATA, rownames(DATA) %in% MyGenes)
  DATA <- log(topGenes + 0.1)
  DATA <- dist(scale(t(DATA)), method = dist.method)
  hc <- hclust(DATA, method = clust.method)
  if (type == "classic") {
    return(plot(hc, hang = hang, ylab = "Height", xlab = "Clusters", sub=""))
  } else {
    par("mar")
    par(mar=c(0,0,0,0))
    return(plot(as.phylo(hc),
         type = type,
         cex = cex,
         label.offset = label.offset))
  }
}
