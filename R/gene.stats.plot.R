#' Scale data
#'
#' This function takes an object of class scSeqR and scales the normalized data.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' gene.stats.plot(my.obj)
#' }
#' @export
gene.stats.plot <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@gene.data
  x = c(1:100)
  y = log1p(1:100)
  myPlot <- ggplot(DATA,aes(x=log2(numberOfCells),y=SDs, color = log2(SDs))) +
    theme_bw() +
    geom_point(size = 1, alpha = 0.5) +
    scale_y_continuous(trans = "log1p")
  return(myPlot)
}
