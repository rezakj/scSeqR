#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class scSeqR and creats plot.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' tsne.plot(my.obj)
#' }
#' @export
tsne.plot <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@tsne.data
  myPLOT <- ggplot(DATA, aes(x = V1, y = V2,
                             text = cells ,
                             color = cl_hierarchical,
                             alpha = 0.5)) +
    geom_point(size=1) +
    xlab("Dim1") +
    ylab("Dim2") +
    ggtitle("t-SNE plot") + theme_bw()
  return(myPLOT)
}
