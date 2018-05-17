#' Scale data
#'
#' This function takes an object of class scSeqR and scales the normalized data.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' gene.stats(my.obj)
#' }
#' @export
gene.stats <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@main.data
  mymat = as.matrix(DATA)
  SDs <- apply(mymat, 1, function(mymat) {sd(mymat)})
  Table <- list(row.names(DATA),
             as.numeric(rowSums(DATA > 0)),
             as.numeric(rowMeans(DATA)),
             as.numeric(SDs))
  names(Table) <- c("genes","numberOfCells","meanExp","SDs")
  Table <- as.data.frame(Table)
  attributes(x)$gene.data <- Table
  return(x)
}
