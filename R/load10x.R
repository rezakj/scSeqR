#' Load 10X data as dgTMatrix
#'
#' This function takes 10X data files barcodes.tsv, genes.tsv and matrix.mtx and converts them to proper matrix file for scSeqR.
#' @param dir.10x A directory that includes the 10X barcodes.tsv, genes.tsv and matrix.mtx files.
#' @return The dgTMatrix matrix.
#' @examples
#' \dontrun{
#' load10x('/hg19')
#' }
#' @export
load10x <- function (dir.10x = NULL) {
  if (!dir.exists(dir.10x)) {
    stop("Directory is not provided. Please provide a standard 10x matrix directory
         that includes; barcodes.tsv, genes.tsv and matrix.mtx files")
  }
  if (dir.exists(dir.10x)) {
  Standard10xInput <- paste(c("barcodes.tsv","genes.tsv","matrix.mtx"),collapse="")
  InputFiles <- paste(list.files(dir.10x),collapse="")
  }
  if (Standard10xInput != InputFiles) {
    stop("Provided directory does not have a standard 10x matrix directory
         that includes; barcodes.tsv, genes.tsv and matrix.mtx files")
  }
  if (Standard10xInput == InputFiles) {
    MTX10x <- list.files(dir.10x,full.names=T,"matrix.mtx")
    cell.barcodes <- list.files(dir.10x,full.names=T,"barcodes.tsv")
    gene.names.ids <- list.files(dir.10x,full.names=T,"genes.tsv")
    MTX10x <- readMM(MTX10x)
    cell.barcodes <- readLines(cell.barcodes)
    gene.names.ids <- readLines(gene.names.ids)
  }
  colnames(x = MTX10x) <- cell.barcodes
  rownames(x = MTX10x) <- gene.names.ids
  data.10x <- append(x = MTX10x, values = MTX10x)
  return(data.10x)
}
