#' Normalize data
#'
#' This function takes an object of class scSeqR and normalized the data based on "global.glsf", "ranked.glsf" or "spike.in" methods.
#' @param x An object of class scSeqR.
#' @param norm.method Choose a normalization method, there are three option currently.
#' Choose from "global.glsf", "ranked.glsf" or "spike.in", defult = "ranked.glsf".
#' @param top.rank If the method is set to "ranked.glsf", you need to set top number of genes sorted based on global base mean, defult = 500.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' norm(my.obj, "ranked.glsf", top.rank = 500)
#' }
#' @import gmp
#' @export
norm.data <- function (x = NULL, norm.method = "ranked.glsf", top.rank = 500) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (norm.method == "global.glsf") {
    DATA <- x@main.data
    libSiz <- colSums(DATA)
    factr <- max(factorize(libSiz))
    GLSF <- as.numeric(libSiz) / as.numeric(factr)
    decim <- nchar(as.character(round(max(GLSF))))-1
    FGLSF <- as.numeric(GLSF) / 10^decim
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, FGLSF, `/`))
  }
  if (norm.method == "ranked.glsf") {
    DATA <- x@main.data
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes = head(raw.data.order,top.rank)
    libSiz <- colSums(topGenes)
    factr <- max(factorize(libSiz))
    GLSF <- as.numeric(libSiz) / as.numeric(factr)
    decim <- nchar(as.character(round(max(GLSF))))-1
    topFGLSF <- as.numeric(GLSF) / 10^decim
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, topFGLSF, `/`))
  }
  attributes(x)$main.data <- normalized
  return(x)
}

