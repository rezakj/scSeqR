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
stats.plot <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@stats
  mito.percent.plot <- ggplot(DATA,aes(y=mito.percent,x=".", alpha = 0.5)) +
    theme_bw() +
    geom_jitter(color = "red") +
    geom_boxplot() + xlab("mito.percent") + ylab("percent of mito genes per cell")
  nGenes.plot <- ggplot(DATA,aes(y=nGenes,x=".", alpha = 0.5)) +
    theme_bw() +
    geom_jitter(color = "red") +
    geom_boxplot() + xlab("nGenes") + ylab("number of genes per cell")
  UMIsplot <- ggplot(DATA,aes(y=UMIs,x=".", alpha = 0.5)) +
    theme_bw() +
    geom_jitter(color = "red") +
    geom_boxplot() + xlab("UMIs") + ylab("number of UMIs per cell")
  return(grid.arrange(nGenes.plot, UMIsplot,mito.percent.plot, ncol = 3))
}
