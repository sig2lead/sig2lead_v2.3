fpsim_cluster <- function(fpset_added, lsm_rows, all_cmpd){
  load("./lincs_fps.RData")
  #fpSim_start_time <- Sys.time()
  #lsm_rows <- 1:100
  lincs <- lincs_fps_2[lsm_rows,]
  fpset_input <- rbind(fpset_added, lincs)
  
  if (nrow(fpset_input) < 500){
    fpset_input <- fpset_input
  }
  else if (((nrow(fpset_input) >= 500) & (all_cmpd==0))){
    fpset_input <- fpset_input[1:500,]
  }
  else if (((nrow(fpset_input) >= 500) & (all_cmpd==1))){
    fpset_input <- fpset_input
  }
  
 
  
 # if ((rownames(fpset_input) != ""))){
  simMA <<- sapply(rownames(fpset_input), function(x) fpSim(fpset_input[x,], fpset_input, sorted=FALSE))
#  }
#  else
#  {
#    simMA <<- fpSim(fpset_input, sorted=FALSE)
#  }  
  hc <<-hclust(as.dist(1-simMA), method="average")
  #hc <<-hclust(as.dist(1-simMA), method="complete")
  par(cex=0.3)
  #simmA_fpSim <-sapply(rownames(fpset_added), function(x) max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  #lincs_fpsim <-sapply(rownames(fpset_added), function(x) which.max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  #lincs_fpsim_2 <- lincs_fps_2[lincs_fps_2]
  #fpSim_end_time <- Sys.time()
  #total_fpSim_time[[j]] <- fpSim_end_time - fpSim_start_time
  #print(paste("Total time to run fpSim: ", total_fpSim_time[[j]], sep=""))
  
  #lincs_fp_ind <- lsm_rows[lincs_fpsim]
  #lincs_vec <- rownames(lincs_fps_2)[lincs_fp_ind]
  #simmA_fpSim_df <- data.frame(lincs_vec,simmA_fpSim)
  #return(simmA_fpSim_df)
  
  
 
}

fpsim_cluster_no_added <- function(lsm_rows, all_cmpd){
  load("./lincs_fps.RData")
  #fpSim_start_time <- Sys.time()
  #lsm_rows <- 1:100
  lincs <- lincs_fps_2[lsm_rows,]
  #fpset_input <- rbind(fpset_added, lincs)
  fpset_input <- lincs
  
  if (nrow(fpset_input) < 5000){
    fpset_input <- fpset_input
  }
  else if (((nrow(fpset_input) > 5000) & (all_cmpd==0))){
    fpset_input <- fpset_input[1:5000,]
  }
  else if (((nrow(fpset_input) > 5000) & (all_cmpd == 1))){
    fpset_input <- fpset_input
  }
  
  
  
  
  simMA <<- sapply(rownames(fpset_input), function(x) fpSim(fpset_input[x,], fpset_input, sorted=FALSE))
  hc <<-hclust(as.dist(1-simMA), method="average")
  #hc <<-hclust(as.dist(1-simMA), method="complete")
  par(cex=0.3)
  #simmA_fpSim <-sapply(rownames(fpset_added), function(x) max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  #lincs_fpsim <-sapply(rownames(fpset_added), function(x) which.max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  #lincs_fpsim_2 <- lincs_fps_2[lincs_fps_2]
  #fpSim_end_time <- Sys.time()
  #total_fpSim_time[[j]] <- fpSim_end_time - fpSim_start_time
  #print(paste("Total time to run fpSim: ", total_fpSim_time[[j]], sep=""))
  
  #lincs_fp_ind <- lsm_rows[lincs_fpsim]
  #lincs_vec <- rownames(lincs_fps_2)[lincs_fp_ind]
  #simmA_fpSim_df <- data.frame(lincs_vec,simmA_fpSim)
  #return(simmA_fpSim_df)
  
  
  
}


# cluster <- function()
# {
#   sdfset <- sdf_smiles
#   #sdfset <- read.SDFset("sdf_smiles.sdf")
#   #sdfset <- read.SDFset("lib/A1TestedWithiLINCSExtract")
#   #sdfset <- sdf_test
#   valid <- validSDF(sdfset)
#   unique_ids <- makeUnique(sdfid(sdfset))
#   cid(sdfset) <- unique_ids
#   rpropma <- data.frame(MF=MF(sdfset), MW=MW(sdfset))
#   #plot(sdfset[1:3], print=FALSE)
#   
#   
#   apset <<-sdf2ap(sdfset)
#   fpset<<-desc2fp(apset)
#   print("Calculating Similarity Matrix")
#   #begin_sim_time <- Sys.time()
#   simMA <<- sapply(cid(fpset), function(x) fpSim(fpset[x], fpset, sorted=FALSE))
#   
#   hc <<-hclust(as.dist(1-simMA), method="average")
#   #hc <<-hclust(as.dist(1-simMA), method="complete")
#   par(cex=0.3)
#   #plot(as.dendrogram(hc),edgePar=list(col=4, lwd=2),horiz=TRUE)
#   
#   #heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5)
#   
#   #cut.4<-cutree(hc, h=0.4)
#   #cut.2<-cutree(hc, h=0.2)
#   
# }