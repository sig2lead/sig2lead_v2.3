GetClusteredAddedCompounds <- function(InFile){
  indices <- list()
  h <- 1
  LINCSCluster <- (InFile$LSM.Compounds > 0)
  for (i in 1:length(LINCSCluster)){
    if (LINCSCluster[[i]] == TRUE && InFile$Non.LSM.compounds[[i]] > 0){
      indices[[h]] <- InFile$Cluster[[i]]
      h <- h+1
    }
  }
  indices <<- indices
AddedInClusters <- list()
g <- 1
  for (j in 1:length(test$clusterID)){
          if (test$clusterID[[j]] %in% indices && !grepl("LSM-", test$Compound[[j]])){
           AddedInClusters[[g]] <- as.character(test$Compound[[j]])
           g<- g+1
          }
    }
AddedInClusters <<- unlist(AddedInClusters)
}

ClusteredAdds_SMILES <- function(){
  h <- 1
  SMILESClusterAdd <- list()
  IDClusterAdd <- list()
  SMILESwithIDAdds <- data.frame()
  AddsInClusterBinary <- adds_csv$LSM_ID %in% AddedInClusters
  for (i in 1:length(AddsInClusterBinary)){
  if (AddsInClusterBinary[[i]] == TRUE){
    SMILESClusterAdd[[h]] <- as.character(adds_csv$SMILES[[i]])
    IDClusterAdd[[h]] <- as.character(adds_csv$LSM_ID[[i]])
    h <- h+1
  }
  }
  SMILESClusterAdd <<- SMILESClusterAdd
  IDClusterAdd <<- IDClusterAdd
  SMILESwithIDAdds <- cbind("SMILES" = unlist(SMILESClusterAdd), "Compound" = unlist(IDClusterAdd))
  SMILESwithIDAdds <<- SMILESwithIDAdds
  write.table(SMILESwithIDAdds, file = "SMILESforAdds.smi", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
