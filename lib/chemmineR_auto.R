###############################
# ChemmineR Clustering of Analog
##############################
library(ChemmineR)
#library(help="ChemmineR")
library("gplots")
library("cluster")
library("ggplot2")
#source("lib/minSim.R", local=TRUE)

cluster_compounds <- function(lsm_rows, lsm_labels_heatmap, fpset_added_lsm, input_lsm, m_chemmineR)
  {
  
  print("Starting to Calculate Similarity Matrix")
  #begin_sim_calc <- Sys.time()
  
  
  
  begin_sim_calc <- Sys.time()
  
################ LSM to LSM #################################################### 
  
  
  tanimoto_plusplus_lsm <- list()
  simMA_lsm <- replicate(2594, 0)
  tanimoto_final_lsm <- list()
  
  
  for(l in 1:length(lsm_rows)){
    #l <- 696
    query_1_lsm <- (which(input[(lsm_rows[l] * 2),c(3:1026)] == 1))
    #query_1_lsm <- (which(input[(lsm_rows[l] * 2),] == 1))
    sim_lsm <- vector("integer", length=42764)
    #sim_lsm <- 0
    m_query_lsm <- length(query_1_lsm)
    
    
    for(h in 1:length(query_1_lsm)){
      #sim[sets_minority_cols[[query_1[]]][l]] <- sim[sets_minority_cols[[j]][l]] + 1
      #sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
      
      #sim_lsm[sets_minority_cols[[query_1_lsm[h]]]] <- sim_lsm[sets_minority_cols[[query_1_lsm[h]]]] + 1
      sim_lsm[sets_minority_cols[[query_1_lsm[h]]]] <- sim_lsm[sets_minority_cols[[query_1_lsm[h]]]] + 1
    }
    
    #end_sim_calc <- Sys.time()
    
    #begin_sim_calc <- Sys.time()
    tanimoto_plusplus_lsm[[l]] = sim_lsm/ (m + m_query_lsm - sim_lsm)
    
    #end_sim_calc <- Sys.time()
    
    tanimoto_final_lsm[[l]] <- tanimoto_plusplus_lsm[[l]][lsm_rows]
    simMA_lsm <- cbind(simMA_lsm, tanimoto_final_lsm[[l]])  
  } 

  
################################ Added Compounds to LSM #########################################  
tanimoto_plusplus_added_lsm <- list()
simMA_added_lsm <- replicate(2594, 0)
tanimoto_final_added_lsm <- list() 

fpset_added_lsm <- fpset_added
 for(k in 1:dim(fpset_added_lsm)[1]){
   
    query_1_a_l <- which(fpset_added_lsm[k,] == 1)
    sim_a_l <- vector("integer", length=42764)
    m_query_a_l <- length(query_1_a_l)
    
   
    for(j in 1:length(query_1_a_l)){
      #sim[sets_minority_cols[[query_1[]]][l]] <- sim[sets_minority_cols[[j]][l]] + 1
      #sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
      sim_a_l[sets_minority_cols[[query_1_a_l[j]]]] <- sim_a_l[sets_minority_cols[[query_1_a_l[j]]]] + 1
    }
    #end_sim_calc <- Sys.time()
    
    #begin_sim_calc <- Sys.time()
    tanimoto_plusplus_added_lsm[[k]] = sim_a_l / (m + m_query_a_l - sim_a_l)
    
    #end_sim_calc <- Sys.time()
    
    tanimoto_final_added_lsm[[k]] <- tanimoto_plusplus_added_lsm[[k]][lsm_rows]
    simMA_added_lsm <- cbind(simMA_added_lsm, tanimoto_final_added_lsm[[k]])  
  }  
######################################################################################## 

###################### Added Compounds to Added Compounds ############################################
fpset_added_added <- fpset_added
# simMA_cluster <-  apply(fpset_combined, MARGIN=1, FUN=function(x) {minSim(x, lsm_rows[1:2], m_chemmineR)})
# simMA_added_added <-  apply(fpset_added_added, MARGIN=1, FUN=function(x) {fpSim(fpset[x],fpset, method="Tanimoto")})

fpSim_added <- list()
for (i in 1:nrow(fpset_added_added)){ 
 fpSim_added[[i]] <- fpSim(x=fpset_added_added[i,] , y=fpset_added_added , method="Tanimoto", sorted=FALSE)
}

simMA_Added_added <- do.call(rbind,fpSim_added)
#simMA_Added_added <- matrix(unlist(fpSim_added), byrow=FALSE, nrow=length(fpSim_added) ) 
 # end_sim_calc <- Sys.time()
# end_sim_calc <- Sys.time()
#simMA_added_added <- apply(fpset_added_added, 1,  function(x) fpSim())

################################## Combine all Similarity Matrices ###################################################################
  
  simMA_added_added_lsm  <- rbind(simMA_Added_added, simMA_added_lsm[,-1])
  simMA_1 <- cbind(simMA_added_lsm[,-1], simMA_lsm[,-1])
  
  simMA <- rbind(t(simMA_added_added_lsm), simMA_1)
  
  added_labels <- colnames(simMA_added_added_lsm)
  all_labels <- c(added_labels,lsm_labels_heatmap)
  
  #write.csv(simMA, file="./example_simMA.csv")
#####################################################################################################################################3
  
  end_sim_calc <- Sys.time()
  #colnames(simMA) <- sdfid(sub  
  #simMA <<- fpSim(x=lsm_compounds, y=added_compounds, method="Tanimoto")
  #end_sim_calc <- Sys.time()
  total_sim_time <<- end_sim_calc - begin_sim_calc
  print(paste("Total Time to Calculate Similarity Matrix: ", total_sim_time, sep= " "))

########################### Generate Heatmap ###############################################################################################  
  hc <<-hclust(as.dist(1-simMA), method="average")
  #hc <<-hclust(as.dist(1-simMA), method="complete")
  #par(cex=0.3)
  #plot(as.dendrogram(hc),edgePar=list(col=4, lwd=2),horiz=TRUE)
  
  windows()
  heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=all_labels, labRow=all_labels, cexRow=0.5, cexCol=0.5)
  
  #cut.4<-cutree(hc, h=0.4)
  #cut.2<-cutree(hc, h=0.2)
  
#sdfset <- sdf_smiles
#sdfset <- read.SDFset("sdf_smiles.sdf")
#sdfset <- read.SDFset("lib/A1TestedWithiLINCSExtract")
#sdfset <- sdf_test
#valid <- validSDF(sdfset)
##unique_ids <- makeUnique(sdfid(sdfset))
#cid(sdfset) <- unique_ids
#rpropma <- data.frame(MF=MF(sdfset), MW=MW(sdfset))
#plot(sdfset[1:3], print=FALSE)



# apset_lsm <-sdf2ap(lsm_sdf)
# fpset_lsm <-as.matrix(desc2fp(apset_lsm))
# lsm_labels <- input[seq(from=1,to=85527,by=2),1]
# lsm_labels <- gsub(">","",lsm_labels)
# lsm_rows <- which(lsm_labels %in% rownames(fpset_lsm)) 
# 
# 
# apset <-sdf2ap(sdfset)
# fpset<<-desc2fp(apset)
# 
# apset_added <-sdf2ap(added_cmpd_sdf)
# fpset_added <<-as.matrix(desc2fp(apset_added))
# 
# fpset_combined <- rbind(fpset_lsm[1:2,], fpset_added[1:2,])
# print("Calculating Similarity Matrix")
# #begin_sim_time <- Sys.time()
# #simMA <<- sapply(cid(fpset), function(x) fpSim(fpset[x], fpset, sorted=FALSE))
# begin_sim_calc <- Sys.time()
# simMA_cluster <-  apply(fpset_combined, MARGIN=1, FUN=function(x) {minSim(x, lsm_rows[1:2], m_chemmineR)})
# end_sim_calc <- Sys.time()
# total_sim_time <<- end_sim_calc - begin_sim_calc
# 
# print(paste("Total Time to Calculate Similarity Matrix: ", total_sim_time, sep= " "))
# hc <<-hclust(as.dist(1-simMA), method="average")
# #hc <<-hclust(as.dist(1-simMA), method="complete")
# par(cex=0.3)
#plot(as.dendrogram(hc),edgePar=list(col=4, lwd=2),horiz=TRUE)

#heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5)

#cut.4<-cutree(hc, h=0.4)
#cut.2<-cutree(hc, h=0.2)

}


#lsm_sdf <- sdf_smiles
#lsm_labels <- input[seq(from=1,to=85527,by=2),1]
#lsm_labels <- gsub(">","",lsm_labels)
#lsm_rows <- which(lsm_labels %in% rownames(fpset_lsm)) 
# m_chemmineR <- fp_data_3$Sum
# added_cmpd_sdf <- subset_add[1]

#test_bypass <- bypass_clustering(lsm_sdf, added_cmpd_sdf)

bypass_clustering <- function(lsm_rows, fpset_added,m_chemmineR)
{
  #sdfset <- sdf_smiles_2
  #sdfset <- read.SDFset("sdf_smiles.sdf")
  #sdfset <- read.SDFset("lib/A1TestedWithiLINCSExtract")
  #sdfset <- sdf_test
  #valid <- validSDF(sdfset)
 # unique_ids <- makeUnique(sdfid(sdfset))
  #cid(sdfset) <- unique_ids
  #rpropma <- data.frame(MF=MF(sdfset), MW=MW(sdfset))
  #plot(sdfset[1:3], print=FALSE)
  
  
  #apset <-sdf2ap(sdfset)
  #fpset<<-desc2fp(apset)
  
  #apset_lsm <-sdf2ap(lsm_sdf)
  #fpset_lsm <-as.matrix(desc2fp(apset_lsm))
  #lsm_labels <- input[seq(from=1,to=85527,by=2),1]
  #lsm_labels <- gsub(">","",lsm_labels)
  #lsm_rows <- which(lsm_labels %in% rownames(fpset_lsm)) 
  
  #apset_added <-sdf2ap(added_cmpd_sdf)
  #fpset_added <<-as.matrix(desc2fp(apset_added))
  
  
  #fpset_added_2 <- fpset_added[1,]
 
 ######################### minSim ########################################################## 
  #simMA <-  apply(fpset_added, MARGIN=1, FUN=function(x) {minSim(x, lsm_rows, m_chemmineR)})
  
  #simMA <<- sapply(cid(fpset_added), function(x) minSim(fpset_added[x,], fpset_lsm, sorted=FALSE))
  print("Starting to Calculate Similarity Matrix")
  #begin_sim_calc <- Sys.time()
  
  
  tanimoto_plusplus <- list()
  #simMA <- replicate(length(lsm_rows), 0)
  simMA <- replicate(42764, 0)
  tanimoto_final <- list()
  
  begin_sim_calc <- Sys.time()
  for(k in 1:dim(fpset_added)[1]){
  #lsm_labels <- input[seq(from=1,to=85527,by=2),1]
  #lsm_labels <- gsub(">","",lsm_labels)
    query_1 <- which(fpset_added[k,] == 1)
    sim <- vector("integer", length=42764)
    m_query <- length(query_1)
    #print(paste("m_query", m_query))
    if(m_query == 0){
      #print(paste("k", k))
      next
    }
  
  # query_1 <- which(query[] == 1)
  # m_query <- sum(query_1)
  # 
  #query_labels <- rownames(query_mat)
  #lsm_labels <- rownames(lsm_mat)
  
    #begin_sim_calc <- Sys.time()
    for(j in 1:length(query_1)){
      #sim[sets_minority_cols[[query_1[]]][l]] <- sim[sets_minority_cols[[j]][l]] + 1
      #sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
      sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
    }
    #end_sim_calc <- Sys.time()
    
    #begin_sim_calc <- Sys.time()
    tanimoto_plusplus[[k]] = sim / (m + m_query - sim)
   
    #end_sim_calc <- Sys.time()

    
    # tanimoto_final[[k]] <- tanimoto_plusplus[[k]]
     
    tanimoto_final[[k]] <- tanimoto_plusplus[[k]][tanimoto_plusplus[[k]] > 0.8]
     
     
     #simMA <- cbind(simMA, tanimoto_final[[k]])  
     #simMA <- cbind(simMA, tanimoto_final[[k]])  
}  

end_sim_calc <- Sys.time()
  #colnames(simMA) <- sdfid(subset_add)
  simMA <<- simMA
tanimoto_final <<- tanimoto_final
  #simMA <<- fpSim(x=lsm_compounds, y=added_compounds, method="Tanimoto")
  #end_sim_calc <- Sys.time()
  total_sim_time <<- end_sim_calc - begin_sim_calc
  print(paste("Total Time to Calculate Similarity Search: ", total_sim_time, sep= " "))
  #hc <<-hclust(as.dist(1-simMA), method="average")
  #hc <<-hclust(as.dist(1-simMA), method="complete")
  #par(cex=0.3)
  #plot(as.dendrogram(hc),edgePar=list(col=4, lwd=2),horiz=TRUE)
  
  #heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5)
  
  #cut.4<-cutree(hc, h=0.4)
  #cut.2<-cutree(hc, h=0.2)
  
  ######################## End minSim ################################################################
  
}
####################################
#Get Clusters
####################################

#rect.hclust(hc, h=.4, border="red")

#cut.4[hc$order]
#ClusterMembers<-cbind(clusterID=cut.4)
#cluster1 <- 1-simMA[cut.4==49]
#cluster1

#ClusterMembers[hc$order,]

      
