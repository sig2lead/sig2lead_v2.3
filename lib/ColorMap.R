color_dend <- function()
{
  dist_mat <<- 1-simMA
  hc <<-hclust(as.dist(dist_mat), method="average")
  dend <<- as.dendrogram(hc)
#  P2 <- read.csv("OrderedP2.csv")
#  P4 <- read.csv("OrderedP4.csv")
  
  #added_compounds <- rep("Other", length(rownames(dist_mat)))
  #is_added <- grepl("ZINC", rownames(dist_mat))
  #added_compounds[is_added] <- "ZINC"
  #added_compounds <- factor(added_compounds)
  #n_added_compounds <- length(unique(added_compounds))
  #cols_4 <- colorspace::terrain_hcl(n_added_compounds, c=150, l=50)
  #cols_4 <- colorspace::rainbow_hcl(n_added_compounds, c=150, l=50)
  #col_added_compounds <-cols_4[added_compounds]
  
  heatmap_colors <- rep('black', nrow(dist_mat))
  heatmap_colors[rownames(dist_mat) %in% cid(adds_SMI)] <- 'magenta'
  #heatmap_colors[rownames(dist_mat) %in% P4[[1]]] <- 'blue'
  #heatmap_colors[is_added == TRUE] <- 'red'
  
  #colors <- as.numeric(dist[,5])
  #colors <- colors[order.dendrogram(dend)]
  
  #labels_colors(dend)[] <- 'black'
  #labels_colors(dend)[rownames(dist_mat) %in% P2[[1]]] <- 'red'
  #labels_colors(dend)[rownames(dist_mat) %in% P4[[1]]] <- 'blue'
  
  #labels_colors(dend) <- col_added_compounds[order.dendrogram(dend)]
  #par(cex=0.2)
  #plot(dend,edgePar=list, lwd=.2, horiz=TRUE)
  
  #test2 <<-heatmap.2(dist_mat, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), colRow = heatmap_colors, colCol = heatmap_colors, col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5)
  output$distPlot<<- renderPlot(heatmap.2(dist_mat, Rowv=dend, Colv=dend, colRow = heatmap_colors, colCol = heatmap_colors, col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
}