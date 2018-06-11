#' Find optimal number of clusters.
#'
#' This function takes an object of class scSeqR and finds optimal number of clusters based on three methods.
#' @param x An object of class scSeqR.
#' @param clust.num Number of clusters
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' gene.plot(my.obj, gene = "NULL", box.to.test = 0, box.pval = "sig.signs")
#' }
#' @import ggpubr
#' @import gridExtra
#' @export
gene.plot <- function (x = NULL,
                       gene = "NULL",
                       box.to.test = 0,
                       box.pval = "sig.signs") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (box.to.test == 0) {
    stop("provide cluster number to perform test")
  }
  if (gene == "NULL") {
    stop("There is no gene name provided. Please provide a gene name")
  }
  DATA <- x@main.data
  AllGenes = row.names(DATA)
  gene.availability = gene %in% AllGenes
  if(gene.availability != TRUE)
  {
    stop("Your gene name is not in the main data. To see the gene names issue this command:
         row.names(YOURobject@main.data)")
  }
  # get tSNE dementions
  datatsne <- x@tsne.data
  row.names(datatsne) <- datatsne$cells
  clusters <- as.character(datatsne$cl_hierarchical)
  # get the gene from the main data
  sub.data <- subset(DATA,rownames(DATA) == gene)
  data.t <- t(sub.data)
  data.expr <- as.data.frame(data.t)
  data.binary <- data.t > 0
  data.binary <- as.data.frame(data.binary)
  my.gene = data.binary[[gene]]
  # plot tSNE
  tSNEplot.heat <- ggplot(datatsne, aes(x = V1, y = V2, col = my.gene)) +
    geom_point(size=1.75,alpha = 0.5, show.legend = F) +
    xlab("Dim1") +
    ylab("Dim2") +
    ggtitle(gene) +
    theme_bw()
  # plot box plot
  Yaxis1 = data.expr[[gene]]
  Yaxis = Yaxis1 + 1
  Yaxis = log2(Yaxis)
  Yaxis[mapply(is.infinite, Yaxis)] <- 0

  bp2 <- ggplot(data.binary, aes(x=clusters,y=Yaxis ,fill = clusters)) +
    theme_bw() +
    geom_jitter(col = "slategray3" ,alpha = 0.5,show.legend = F) +
    xlab("clusters") +
    ggtitle(gene) +
    geom_boxplot(col = "green", notch = T, alpha = 0.6,show.legend = F,outlier.shape = NA) +
    ylab("scaled normalized expression")
  # add p-val
  if (box.pval == "sig.signs") {
    bp2 <- bp2 + stat_compare_means(label = "p.signif", ref.group = box.to.test)
  }
  if (box.pval == "sig.values") {
    bp2 <- bp2 + stat_compare_means(aes(label = paste0("p = ", ..p.format..)),ref.group = box.to.test)
  }
  # Bar plot
  BarP1 <- ggplot(data.expr, aes(x = clusters,y = Yaxis1 ,fill = clusters)) +
    stat_summary(fun.y="mean", geom="bar",alpha = 0.5,show.legend = F) +
    ylab("avraged normalized expression") +
    ggtitle(gene) +
    theme_bw()
  # return plots
  return(grid.arrange(arrangeGrob(tSNEplot.heat),
                      arrangeGrob(bp2,BarP1, ncol=1),
                      ncol=2, widths=c(1.5,0.9)))
}
