MDS_plot <- function(centroid_list)
{
  x <- 0
  indices <- list()
  MDS_dists <- matrix(0.0, nrow = length(centroid_list), ncol = length(centroid_list))
  
  for (i in 1:length(centroid_list))
  {
    x<- x+1
    indices[[x]] <- grep(centroid_list[[i]], row.names(dist_mat))
    #    print(indices)
    if (length(indices[[i]])>1){  
      indices[[i]] <- indices[[i]][[1]]
    }
  }
  for (j in 1:(length(indices)-1)){
    for (k in (j+1):length(indices)){
      MDS_dists[j, k] <- dist_mat[indices[[j]],indices[[k]]]
      MDS_dists[k, j] <- dist_mat[indices[[j]],indices[[k]]]
      
      #  print(indices)
    }
  }
  fit <- cmdscale(MDS_dists,eig=TRUE, k=2)
  #  print(fit$points[,1])
  #  print(fit$points[,2])
  x1 <- fit$points[,1]
  y1 <- fit$points[,2] 
  #  print(class(x1))
  #  fit <<-fit
  cluster_info$x1 <- x1 
  cluster_info$y1 <- y1 
  #cluster_info$P2Count <- as.numeric(unlist(cluster_info$P2Count, recursive = TRUE))
  #cluster_info$P4Count <- as.numeric(unlist(cluster_info$P4Count, recursive = TRUE))
  cluster_info$LINCS <- as.numeric(unlist(cluster_info$LINCS, recursive = TRUE))
  #cluster_info$input$AddedLabel <- as.numeric(unlist(cluster_info$input$AddedLabel), recursive = TRUE)
  
  
  cluster_info$Added <- as.numeric(unlist(cluster_info$Added, recursive = TRUE))
  #test <- cluster_info
  
#  cluster_info$Added <- 0
  #cluster_info <<- test
  cluster_info <<-cluster_info
  #  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  #       main="Metric MDS", type="p")
  
  #  plot(x1, y1, xlab="Coordinate 1", ylab="Coordinate 2",
  #       main="Metric MDS", type="n")
  #  text(x1, y1, labels = centroid_list[1:length(centroid_list)], cex=.7) 
  
  #  test <- plot(x1, y1, xlab="Coordinate 1", ylab="Coordinate 2",
  #       main="Metric MDS", type="p") 
  #  test <- plot(x1, y1, xlab="Coordinate 1", ylab="Coordinate 2",
  #       main="Metric MDS", type="t") 
  #  ggplot() + geom_scatterpie(mapping = aes(x=x1, y=y1), data=cluster_info, cols=c("P2Count", "P4Count", "LINCS")) + coord_fixed()
  cluster_info$radius <- (as.numeric(cluster_info$ClusterSize) / 500)
  colnames(cluster_info) <- c("Centroid", "ClusterSize", "My Candidates", "LINCS", "x1", "y1", "radius")
  #ggplot() + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("P2Count", "P4Count", "LINCS")) + coord_fixed()
  #output$MDSPlot<<- renderPlot(ggplot() + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS"))) #"P2Count", "P4Count", "LINCS"))
  cluster_info <<- cluster_info
  x <<- cluster_info$x1
  y <<- cluster_info$y1
}
