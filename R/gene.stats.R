#' Make statistical information for each gene across all the cells (SD, mean, expression, etc.)
#'
#' This function takes an object of class scSeqR and provides some statistical information for the genes.
#' @param x An object of class scSeqR.
#' @param which.data Choose from "raw.data" or "main.data", defult = "raw.data".
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' my.obj <- gene.stats(my.obj, which.data = "main.data")
#' }
#' @export
gene.stats <- function (x = NULL,
                        which.data = "raw.data") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
# get data
  if (which.data == "raw.data") {
    DATA <- x@raw.data
  }
  if (which.data == "main.data") {
    DATA <- x@main.data
  }
# calculate
  mymat = as.matrix(DATA)
  SDs <- apply(mymat, 1, function(mymat) {sd(mymat)})
  Table <- list(row.names(DATA),
             as.numeric(rowSums(DATA > 0)),
             rep(dim(DATA)[2], dim(DATA)[1]),
             (as.numeric(rowSums(DATA > 0)) / dim(DATA)[2])*100,
             as.numeric(rowMeans(DATA)),
             as.numeric(SDs))
  names(Table) <- c("genes","numberOfCells","totalNumberOfCells","percentOfCells","meanExp","SDs")
  Table <- as.data.frame(Table)
  attributes(x)$gene.data <- Table
  return(x)
}
