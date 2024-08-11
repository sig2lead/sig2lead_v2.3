cut_tree <- function(threshold)
{
  cut<-cutree(hc, h=threshold)
  cut[hc$order]
  #print(cut)
  #ClusterMembers <- as.matrix(cut)
  #clusterlabs <<- labels(cut)
  #labs <- row.names(ClusterMembers)
  ClusterMembers <<- cbind(clusterID=cut)
  
  
  #ClusterMembers <- ClusterMembers[mixedsort(ClusterMembers),]
  
  #ClusterMembers[hc$order,]
  #ClusterMembers <<- ClusterMembers
}



find_centroid <- function(cut_tree, cluster_size)
{
  dist_mat<<-1-simMA
  x <- 0
  centroid <- list()
  cluster_counts <- list()
  #P2_counts <- list()
  #P4_counts <- list()
  h <- 0
  
  #P2 <- read.csv("OrderedP2.csv")
  #P2 <- read.csv("SOS1.csv")
  #P4 <- read.csv("OrderedP4.csv")
  
  for (i in 1:length(cut_tree))
  {
    x <- x+1
    y<-0
    z <- 0
    
    newcluster <- list()
    for (j in 1:length(cut_tree)){
      if (cut_tree[[j]]==x)
      {
        y <- y+1
        newcluster[[y]] <- labels(cut_tree)[[1]][[j]]
        #newcluster[[y]] <- cut_tree[[j]][[1]]
      }
    }
    #    z<-0
    if (length(newcluster)>0 && length(newcluster)>cluster_size){
      #print(newcluster)
      ####Do stuff here, find all members of this cluster and measure distances
      centroid_index <- list()
      #print(length(newcluster))
      
      
      for (k in 1:length(newcluster))
      {
        z <- z+1
        centroid_index[[z]]<-grep(newcluster[[k]], row.names(simMA))
        #      print(test)
        if (length(centroid_index[[z]])>1){  
          centroid_index[[z]] <- centroid_index[[z]][[1]]
        }
        #      print(test)
      }
      testfinal <<- centroid_index
      
      #      print(test)
      min_matrix <- rep(0.0, length(centroid_index))
      #      print(test)
      #      print(length(test))
      for (p in 1:(length(centroid_index)-1))
      {
        #        print(test)
        for (q in (p+1):length(centroid_index))
        {
          min_matrix[p] <- min_matrix[p]+dist_mat[centroid_index[[p]], centroid_index[[q]]]
          min_matrix[q] <- min_matrix[q]+dist_mat[centroid_index[[p]], centroid_index[[q]]]
        }
      }
      minimum <<- which.min(min_matrix)
      #      print(newcluster[[minimum]])
      centroid[[x]]<-newcluster[[minimum]]
      #      list.append(centroid, newcluster[[minimum]]
      newcluster<<-newcluster
      
      h <- h+1
      #      print(length(newcluster))
      cluster_counts[[h]] <- length(newcluster)
      #P2_counts[[h]] <- sum(newcluster %in% P2[[1]], na.rm = TRUE)
      #P4_counts[[h]] <- sum(newcluster %in% P4[[1]], na.rm = TRUE)
      #      print(sum(newcluster %in% P2[[1]]), na.rm = TRUE)
      #      print(sum(newcluster %in% P4[[1]]), na.rm = TRUE)
      #      test <<- test
      #print(test)
      centroid<<-centroid
    }
  }
}
  