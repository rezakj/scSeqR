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
stats.plot <- function (x = NULL,
                        cell.color = "slategray3",
                        box.color = "red",
                        cell.size = 1) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@stats
  mito.percent.plot <- ggplot(DATA,aes(y=mito.percent,x=".")) +
    theme_bw() +
    geom_jitter(color = cell.color, size = cell.size, alpha = 0.5) +
    geom_boxplot(fill = box.color, col = "green", notch = T, outlier.shape = NA, alpha = 0.2) +
    xlab("mito.percent") + ylab("percent of mito genes per cell")
  nGenes.plot <- ggplot(DATA,aes(y=nGenes,x=".")) +
    theme_bw() +
    geom_jitter(color = cell.color, size = cell.size, alpha = 0.5) +
    geom_boxplot(fill = box.color, col = "green", notch = T, outlier.shape = NA, alpha = 0.2) +
    xlab("nGenes") + ylab("number of genes per cell")
  UMIsplot <- ggplot(DATA,aes(y=UMIs,x=".")) +
    theme_bw() +
    geom_jitter(color = cell.color, size = cell.size, alpha = 0.5) +
    geom_boxplot(fill = box.color, col = "green", notch = T, outlier.shape = NA, alpha = 0.2) +
    xlab("UMIs") + ylab("number of UMIs per cell")
  return(grid.arrange(nGenes.plot, UMIsplot,mito.percent.plot, ncol = 3))
}
