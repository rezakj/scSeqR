#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class scSeqR and creats plot.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' tsne.plot(my.obj)
#' }
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @export
cluster.plot <- function (x = NULL,
                          cell.size = 1,
                          plot.type = "tsne",
                          cell.color = "black",
                          clustered = TRUE) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clustered == FALSE) {
  if (plot.type == "tsne") {
    DATA <- x@tsne.data
    myPLOT <- ggplot(DATA, aes(x = V1, y = V2,
                               text = row.names(DATA))) +
      geom_point(size = cell.size, col = cell.color, alpha = 0.5) +
      xlab("Dim1") +
      ylab("Dim2") +
      ggtitle("t-SNE plot") + theme_bw()
    return(myPLOT)
  }
  if (plot.type == "pca") {
    DATA <- x@pca.data
    myPLOT <- ggplot(DATA, aes(x = PC1, y = PC2,
                               text = row.names(DATA))) +
      geom_point(size = cell.size, col = cell.color, alpha = 0.5) +
      xlab("Dim1") +
      ylab("Dim2") +
      ggtitle("PCA plot") + theme_bw()
    return(myPLOT)
  }
  if (plot.type == "distance") {
    DATA <- as.matrix(x@dist.data)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    return(pheatmap(DATA,col=colors))
  }
  }
  if (clustered == TRUE) {
    if (plot.type == "tsne") {
      DATA <- x@tsne.data
      myPLOT <- ggplot(DATA, aes(x = V1, y = V2,
                                 text = row.names(DATA),color = clusters)) +
        geom_point(size = cell.size, alpha = 0.5) +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle("t-SNE plot") + theme_bw()
      return(myPLOT)
    }
    if (plot.type == "pca") {
      DATA <- x@pca.data
      myPLOT <- ggplot(DATA, aes(x = PC1, y = PC2,
                                 text = row.names(DATA), col = clusters)) +
        geom_point(size = cell.size, alpha = 0.5) +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle("PCA plot") + theme_bw()
      return(myPLOT)
    }
    if (plot.type == "distance") {
      DATA <- as.matrix(x@dist.data)
      colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
      return(pheatmap(DATA,col=colors))
    }
  }
}
