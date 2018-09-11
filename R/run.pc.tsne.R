#' Run tSNE on the main data
#'
#' This function takes an object of class scSeqR and runs tSNE on the main data.
#' @param x An object of class scSeqR.
#' @param clust.method Choose from "base.mean.rank" or "gene.model", defult is "base.mean.rank".
#' @param top.rank A number taking the top genes ranked by base mean, defult = 500.
#' @param gene.list A list of genes to be used for tSNE analysis. If "clust.method" is set to "gene.model", defult = "my_model_genes.txt".
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' my.obj <- run.tsne(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
#' }
#' @import Rtsne
#' @export
run.pc.tsne <- function (x = NULL,
                      dims = 1:20) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
# get PCA data
  DATA <- x@pca.data
  TransPosed <- DATA[dims]
  #  2 dimention
  #  if (clust.dim == 2) {
  # TransPosed <- t(TopNormLogScale)
  tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 2)
  tsne.data = as.data.frame(tsne$Y)
  tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
  rownames(tsne.data) <- tsne.data$cells
  tsne.data <- tsne.data[,-1]
  attributes(x)$tsne.data <- tsne.data
  #  }
  # choose 3 demention
  # tSNE
  # TransPosed <- t(TopNormLogScale)
  tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 3)
  tsne.data = as.data.frame(tsne$Y)
  tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
  rownames(tsne.data) <- tsne.data$cells
  tsne.data <- tsne.data[,-1]
  attributes(x)$tsne.data.3d <- tsne.data
  return(x)
}
