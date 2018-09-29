#' Create heatmaps for genes in clusters or conditions.
#'
#' This function takes an object of class scSeqR and genes and provides a heatmap.
#' @param x A data frame containing gene counts for cells.
#' @param gene A set of gene names to be heatmapped.
#' @param cluster.by Choose from "clusters" or "conditions", defult = "clusters".
#' @param cluster.rows If set to FALSE the genes would not be clustered, defult = TRUE.
#' @param scale Choose from "row" or "column", defult = "row".
#' @param heat.colors Colors for heatmap, defult = c("blue" ,"white", "red").
#' @return An object of class scSeqR
#' @examples
#' \dontrun{
#' MyGenes <- c("SOD1","CD7")
#' MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.8)
#' heatmap.plot (my.obj, gene = MyGenes, cluster.by = "clusters", cluster.rows = T)
#' }
#' @import pheatmap
#' @export
imm.gen <- function (immgen.data = "rna",
                     gene = "NULL",
                     plot.type = "heatmap",
                     heat.colors = c("blue","white", "red")) {
  ## get main data
#  my.data <- read.delim(x,header=TRUE)
#  rownames(my.data) <- toupper(my.data$gene)
if (immgen.data == "rna") {
  my.data <- immgen.rna
}
  MyRows <- toupper(row.names(my.data))
  MyRows <- gsub("-",".", MyRows)
  rownames(my.data) <- MyRows
#  rownames(my.data) <- make.names(MyRows, unique=TRUE)
#  my.data <- my.data[,-1]
  gene <- toupper(gene)
  my.data <- subset(my.data,row.names(my.data) %in% gene)
  # gg plot
  MYdf <- as.data.frame(sort(colSums(my.data),decreasing = F))
  colnames(MYdf) <- c("MyLev")
  MYdf <- cbind(cells = rownames(MYdf),MYdf)
  MYdf$cells <- factor(MYdf$cells, levels = row.names(MYdf))
  #
  if (plot.type == "point.plot") {
#    mySize = log2(MYdf$MyLev)
    return(ggplot(MYdf, aes(x = log2(MyLev),y=cells, col = log2(MyLev), size =log2(MyLev))) +
             geom_point() +
             xlab("log2(cumulative expression)") +
             ylab("ImmGen RNA-seq Cell Types") +
             theme_classic() +
             scale_size("") +
             scale_colour_gradient(low = "gray", high = "red", name="")
             )
  }
#  ggplot(MYdf, aes(x=cells, y=log2(MyLev))) +
#    geom_bar(stat = "identity") + coord_flip() +
#    theme_classic()
#  my.data <- my.data[ rowSums(counts(my.data)) == 0, ]
  # heatmap colors
  mycol <- colorRampPalette(heat.colors)(n = 100)
  # return
  if (plot.type == "heatmap") {
    return(pheatmap(as.matrix(my.data),
                    col = mycol,
                    show_colnames = T,
                    cluster_rows = T,
                    cluster_cols = T,
                    scale = "row"))
  }
}
