cut_tree <- function(threshold)
{
  cut<-cutree(hc, h=threshold)
  cut[hc$order]
  #print(cut)
  #ClusterMembers <- as.matrix(cut)
  #clusterlabs <<- labels(cut)
  #labs <- row.names(ClusterMembers)
  #ClusterMembers <<- cbind(clusterID=cut)
  
  ClusterMembers <<- cbind(clusterID=cut)
  test4 <- as.data.frame(ClusterMembers)
  test7 <- max(test4$clusterID, na.rm = TRUE)
  test5 <- data.frame(matrix(0,nrow = test7, ncol=2))
  test6 <- c("cluster_number", "cluster_size")
  colnames(test5) <- test6
  test8 <- data.frame(ncol = 2)
  
  for (i in 1:max(test4$clusterID, na.rm = TRUE)){
    test5$cluster_number[i] <- i
    test5$cluster_size[i] <- nrow(filter(test4, clusterID == i))
  }
  test9 <- test5[order(-as.integer(test5$cluster_size)),]
  test10 <- test9
  for (j in 1:test7){
    test10$new_cluster_number[j] <- j
  }
  
  test11 <- ClusterMembers
  for (k in 1:length(ClusterMembers)){
      h <- match(k, test10$cluster_number)
      test11[[k]] <- test10$new_cluster_number
  }
  ClusterMembers <- test11
  ClusterMembers <<- as.matrix(ClusterMembers[order(ClusterMembers),])
  
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
  
  h <- 0
  
  
  
  for (i in 1:length(cut_tree))
    #  for (i in 1:length(cut_tree$new_cluster_number))
  {
    x <- x+1
    y<-0
    z <- 0
    
    newcluster <- list()
    for (j in 1:length(cut_tree)){
      #    for (j in 1:length(cut_tree$new_cluster_number)){
      if (cut_tree[[j]]==x)
        #      if (cut_tree$new_cluster_number[[j]]==x)
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
      #Added_Counts[[h]] <- sum(newcluster %in% cid(adds_SMI), na.rm = TRUE)
      #P2_counts[[h]] <- sum(newcluster %in% P2[[1]], na.rm = TRUE)
      #P4_counts[[h]] <- sum(newcluster %in% P4[[1]], na.rm = TRUE)
      #      print(sum(newcluster %in% P2[[1]]), na.rm = TRUE)
      #      print(sum(newcluster %in% P4[[1]]), na.rm = TRUE)
      #      test <<- test
      #print(test)
      centroid <<- centroid[!sapply(centroid, is.null)]
      cluster_counts <<- cluster_counts[!sapply(cluster_counts, is.null)]
      #P2_counts <<- P2_counts
      #P4_counts <<- P4_counts
      cluster_info <- as.data.frame(cbind(Centroid=unlist(centroid), ClusterSize=as.numeric(unlist(cluster_counts))), stringsAsFactors=FALSE) #P2Count=as.numeric(unlist(P2_counts)), P4Count=as.numeric(unlist(P4_counts))), stringsAsFactors = FALSE)
      LINCSFreq <- as.numeric(cluster_info$ClusterSize)# - as.numeric(cluster_info$P2Count) - as.numeric(cluster_info$P4Count)
      cluster_info$Added <- 0
      cluster_info$LINCS <- LINCSFreq
      cluster_info <<- cluster_info
    }
  }
}

find_centroid_adds <- function(cut_tree, cluster_size)
{
  dist_mat<<-1-simMA
  x <- 0
  centroid <- list()
  cluster_counts <- list()
  Added_Counts <- list()
  #P2_counts <- list()
  #P4_counts <- list()
  h <- 0
  
  adds<-input$AddedCompounds
  adds_SMI<-read.SMIset(paste(getwd(), adds$name, sep="/"))
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
      Added_Counts[[h]] <- sum(newcluster %in% cid(adds_SMI), na.rm = TRUE)
      #P2_counts[[h]] <- sum(newcluster %in% P2[[1]], na.rm = TRUE)
      #P4_counts[[h]] <- sum(newcluster %in% P4[[1]], na.rm = TRUE)
      #      print(sum(newcluster %in% P2[[1]]), na.rm = TRUE)
      #      print(sum(newcluster %in% P4[[1]]), na.rm = TRUE)
      #      test <<- test
      #print(test)
      centroid <<- centroid[!sapply(centroid, is.null)]
      cluster_counts <<- cluster_counts[!sapply(cluster_counts, is.null)]
      Added_Counts <<- Added_Counts
      #P2_counts <<- P2_counts
      #P4_counts <<- P4_counts
      cluster_info <- as.data.frame(cbind(Centroid=unlist(centroid), ClusterSize=as.numeric(unlist(cluster_counts)), Added=as.numeric(unlist(Added_Counts))), stringsAsFactors=FALSE) #P2Count=as.numeric(unlist(P2_counts)), P4Count=as.numeric(unlist(P4_counts))), stringsAsFactors = FALSE)
      LINCSFreq <- as.numeric(cluster_info$ClusterSize)-as.numeric(cluster_info$Added)# - as.numeric(cluster_info$P2Count) - as.numeric(cluster_info$P4Count)
      cluster_info$LINCS <- LINCSFreq
      cluster_info <<- cluster_info
    }
  }
}

centroid_clusters <- function(centroid_list, cluster_members){
  ClusterCentroid <- list()
  for (i in 1:length(centroid_list)){
    ind <- grep(centroid_list[[i]], row.names(cluster_members))
    if (length(ind>1)){
      ind <- ind[1]
    }
    ClusterCentroid[[i]] <- cluster_members[[ind]]
  }
  ClusterCentroid <<- ClusterCentroid
}


#test <- sort.list(cluster_counts)
#test2 <- list()
#for (i in 1:length(test)){
#  test2[[i]] <- centroid[[test[[i]]]]
#}
#test3 <- cluster_info[order(-as.integer(cluster_info$ClusterSize)),]




#test4 <- as.data.frame(ClusterMembers)
#test7 <- max(test4$clusterID, na.rm = TRUE)
#test5 <- data.frame(matrix(0,nrow = test7, ncol=2))
#test6 <- c("cluster_number", "cluster_size")
#colnames(test5) <- test6
#test8 <- data.frame(ncol = 2)

#for (i in 1:max(test4$clusterID, na.rm = TRUE)){
#  test5$cluster_number[i] <- i
#  test5$cluster_size[i] <- nrow(filter(test4, clusterID == i))
#}
#test9 <- test5[order(-as.integer(test5$cluster_size)),]
#test10 <- test9
#for (j in 1:test7){
#  test10$new_cluster_number[j] <- j
#}
