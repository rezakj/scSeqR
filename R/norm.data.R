#' Normalize data
#'
#' This function takes an object of class scSeqR and normalized the data based on "global.glsf", "ranked.glsf" or "spike.in" methods.
#' @param x An object of class scSeqR.
#' @param norm.method Choose a normalization method, there are three option currently.
#' Choose from "global.glsf", "ranked.glsf", "rpm","spike.in" or no.norm, defult = "ranked.glsf".
#' @param top.rank If the method is set to "ranked.glsf", you need to set top number of genes sorted based on global base mean, defult = 500.
#' @param rpm.factor If the norm.method is set to "rpm" the library sizes would be diveded by this number, defults = 1000 (higher numbers recomanded for bulk RNA-Seq).
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' norm(my.obj, "ranked.glsf", top.rank = 500)
#' }
#' @export
norm.data <- function (x = NULL,
                       norm.method = "ranked.glsf",
                       top.rank = 500,
                       spike.in.factors = NULL,
                       rpm.factor = 1000) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@main.data
  if (norm.method == "global.glsf") {
    libSiz <- colSums(DATA)
    norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
  if (norm.method == "ranked.glsf") {
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes = head(raw.data.order,top.rank)
    libSiz <- colSums(topGenes)
    norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
  if (norm.method == "spike.in") {
    norm.facts <- read.table(spike.in.factors,sep="\t")[2]
    norm.facts <- as.numeric(as.matrix(norm.facts))
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
  if (norm.method == "no.norm") {
    norm.facts <- colnames(DATA)
    norm.facts <- norm.facts == 1
    norm.facts[ norm.facts == "FALSE" ] <- 1
    normalized <- DATA
  }
  if (norm.method == "rpm") {
    libSiz <- colSums(DATA)
    norm.facts <- as.numeric(libSiz) / rpm.factor
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
  SizeFactors <- as.numeric(norm.facts)
  names(SizeFactors) <- c(colnames(DATA))
  SizeFactors <- as.data.frame(SizeFactors)
  attributes(x)$main.data <- normalized
  attributes(x)$norm.factors <- norm.facts
  return(x)
}

