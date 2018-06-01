#' Scale data
#'
#' This function takes an object of class scSeqR and scales the normalized data.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' gene.stats.plot(my.obj)
#' }
#'
#' @import ggrepel
#'
#' @export
make.gene.model <- function (x = NULL,
                             dispersion.limit = 1.5,
                             base.mean.rank = 500,
                             non.sig.col = "darkgray",
                             right.sig.col = "chartreuse3",
                             left.sig.col = "cadetblue3",
                             disp.line.col = "black",
                             rank.line.col = "red",
                             cell.size = 1.75,
                             cell.transparency = 0.5,
                             no.mito.model = T,
                             interactive = TRUE,
                             out.name = "plot") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  data <- x@gene.data
#  x = c(1:100)
#  y = log1p(1:100)
# variables
  cellCountLimit = as.numeric(tail(head(data[order(data$meanExp, decreasing = T),],base.mean.rank)[2],1))
  top.rank.line = as.numeric(log2(cellCountLimit))
  SDlimit = dispersion.limit
# add data colors
  data <- data %>%
    mutate(color = ifelse(data$SDs > SDlimit & data$numberOfCells < cellCountLimit,
                          yes = "righLimit",
                          no = ifelse(data$numberOfCells > cellCountLimit,
                                      yes = "leftLimit",
                                      no = "none")))
# get mito genes
  mito.genes <- grep(pattern = "^mt\\-", x = data$genes, value = TRUE, ignore.case = TRUE)
  if ( length(mito.genes) == 0 ) {
    mito.genes <- grep(pattern = "^mt\\.", x = data$genes, value = TRUE, ignore.case = TRUE)
  }
  top_labelled <- subset(data, data$gene %in% mito.genes)
#  plot
  myPlot <- ggplot(data,aes(x=log2(numberOfCells),y=SDs, text = genes, label = genes)) +
    geom_point(aes(color = factor(color)),size = cell.size, alpha = cell.transparency) +
    theme_bw(base_size = 16) +
    theme(legend.position = "none") +
    ggtitle(label = "Dispersion Plot", subtitle = "modeled by basemean, number of cells and dispersion") +  # add title
    xlab("log2(number of cells per gene)") +
    ylab("dispersion") +
    geom_vline(xintercept = top.rank.line, colour = rank.line.col) +
    geom_hline(yintercept = SDlimit, colour = disp.line.col) +
    scale_color_manual(values = c("righLimit" = right.sig.col,
                                  "leftLimit" = left.sig.col,
                                  "none" = non.sig.col)) +
    scale_y_continuous(trans = "log1p")
  # geom_text(aes(label=ifelse(genes  %in% mito.genes ,as.character(mito.genes),'')))
  myPlot <- myPlot + geom_text_repel(data = top_labelled,
                                       mapping = aes(label = genes),
                                       size = 3,
                                       fontface = 'bold',
                                       color = 'black',
                                       box.padding = unit(0.5, "lines"),
                                       point.padding = unit(0.5, "lines"))
# get model genes
  my.clust.genes = subset(data, color != "none")[1]
# exclude mito genes from model genes
  if (no.mito.model == T) {
    my.clust.genes <- subset(my.clust.genes, !(genes %in% mito.genes))
  }
if (interactive == T) {
  OUT.PUT <- paste(out.name, ".html", sep="")
  htmlwidgets::saveWidget(ggplotly(myPlot),OUT.PUT)
}
if (interactive == F) {
  return(myPlot)
}
  # retrun
  # write out gene model
   write.table((my.clust.genes),file="my_model_genes.txt", row.names =F, quote = F, col.names = F)
   print("my_model_genes.txt file is generated, which can be use for clustering")
#  attributes(x)$gene.model <- my.clust.genes
}
