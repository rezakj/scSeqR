#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class scSeqR and creats plot.
#' @param x An object of class scSeqR.
#' @param plot.type Choose from "box.umi", "box.mito", "box.gene", "box.gene.umi.mito", "point.mito.umi", "point.gene.umi".
#' @param cell.color Choose a color for dots.
#' @param cell.size Choose a size for dots.
#' @param cell.transparency Transparency of colors.
#' @param box.color Choose color for box.
#' @param box.line.col Choose a color for line around box plot.
#' @param interactive Make html intractive plot.
#' @param out.name Name for intractive plot html file.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' plot.stats(my.obj)
#' }
#' @export
stats.plot <- function (x = NULL,
                        plot.type = "box.umi",
                        cell.color = "slategray3",
                        cell.size = 1,
                        cell.transparency = 0.5,
                        box.color = "red",
                        box.line.col = "green",
                        back.col = "white",
                        interactive = TRUE,
                        out.name = "plot")
{
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@stats
  w = log1p(1:100)
  # get conditions
  do <- data.frame(do.call('rbind', strsplit(as.character(head(DATA$CellIds,1)),'_',fixed=TRUE)))
  do <- dim(do)[2]
if (do == 2) {
  col.legend <- data.frame(do.call('rbind', strsplit(as.character(DATA$CellIds),'_',fixed=TRUE)))[1]
  col.legend <- factor(as.matrix(col.legend))
} else {
  col.legend = "."
}
  # Box plots
  # mito
  mito.percent.plot <- ggplot(DATA,aes(y=mito.percent,x=col.legend)) +
    theme_bw() +
    geom_jitter(color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = "green", notch = T, outlier.shape = NA, alpha = cell.transparency) +
    xlab("mito.percent") + ylab("percent of mito genes per cell")
    # nGenes
  nGenes.plot <- ggplot(DATA,aes(y=nGenes,x=col.legend)) +
    theme_bw() +
    geom_jitter(color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = box.line.col, notch = T, outlier.shape = NA, alpha = cell.transparency) +
    xlab("nGenes") + ylab("number of genes per cell")
    # UMIs
  UMIsplot <- ggplot(DATA,aes(y=UMIs,x=col.legend)) +
    theme_bw() +
    geom_jitter(color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = box.line.col, notch = T, outlier.shape = NA, alpha = cell.transparency) +
    xlab("UMIs") + ylab("number of UMIs per cell")
# scatter plots
  if (col.legend[1] == ".") {
    Mito.UMIs <- ggplot(DATA,aes(y=mito.percent,x=UMIs,
                                 text = paste("UMIs =",DATA$UMIs,",",DATA$CellIds,sep=" "))) +
      geom_point(color = cell.color, size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col))
    #
    Genes.UMIs <- ggplot(DATA,aes(y=nGenes,x=UMIs,
                                  text = paste("nGenes =",DATA$nGenes,",",DATA$CellIds,sep=" "))) +
      geom_point(color = cell.color, size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_y_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col))
  } else {
    Mito.UMIs <- ggplot(DATA,aes(y=mito.percent,x=UMIs, col = col.legend,
                                 text = paste("UMIs =",DATA$UMIs,",",DATA$CellIds,sep=" "))) +
      geom_point( size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col))
    #
    Genes.UMIs <- ggplot(DATA,aes(y=nGenes,x=UMIs, col = col.legend,
                                  text = paste("nGenes =",DATA$nGenes,",",DATA$CellIds,sep=" "))) +
      geom_point(size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_y_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col))
  }
  # out puts
  if (plot.type == "point.mito.umi") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
    return(Mito.UMIs)
  }
  #
  if (plot.type == "point.gene.umi") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
    return(Genes.UMIs)
  }
  #
  if (plot.type == "box.umi") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
    return(UMIsplot)
  }
  #
  if (plot.type == "box.mito") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
    return(mito.percent.plot)
  }
  #
  if (plot.type == "box.gene") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
    return(nGenes.plot)
  }
  #
  if (plot.type == "box.gene.umi.mito") {
    return(grid.arrange(nGenes.plot, UMIsplot,mito.percent.plot, ncol = 3))
  }
}
