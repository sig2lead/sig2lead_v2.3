fpsim_sim_search <- function(fpset_added, lsm_rows, number_analogs){
  load("./lincs_fps.RData")
  #fpSim_start_time <- Sys.time()
  
  #simmA_fpSim <-sapply(rownames(fpset_added), function(x) max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  #lincs_fpsim <-sapply(rownames(fpset_added), function(x) which.max(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE)))
  #fpset_added <- sdfset_add
  #lsm_rows <- 1:41572
  #number_analogs <- 2
 
  lincs_fpsim <-sapply(rownames(fpset_added), function(x) which.maxn(fpSim(fpset_added[x,], lincs_fps_2[lsm_rows,],sorted=FALSE), n=number_analogs))
  #lincs_fpsim <- as.vector(lincs_fpsim)
  lincs_fpsim <<- lincs_fpsim
  
  
  #for(i in 1:length(fpset_added)){
  #  simma_fpSim_temp <- fpsim(fpset_added[i], lincs_fps_2[c(1:2),i])
    
  #}
  
  #simmA_fpSim <-sapply(rownames(fpset_added), function(x) (fpSim(fpset_added[x,], lincs_fps_2[lincs_fpsim,],sorted=FALSE)))
  simmA_fpSim <- 1:number_analogs
  #test <- 1:2
  for (i in 1:nrow(fpset_added)){
    tmp <- fpSim(fpset_added[i,], lincs_fps_2[lincs_fpsim[c(1:number_analogs),i],])
    #test <- test[lincs_fpsim[1:2,i]]
    #simmA_fpSim <- cbind(simmA_fpSim,test)
    simmA_fpSim <- c(simmA_fpSim,tmp)
  }
  simmA_fpSim <- simmA_fpSim[-c(1:number_analogs)]
  #sim
  #simmA_fpSim <-sapply(rownames(fpset_added), function(x,y) (fpSim(fpset_added[x,], lincs_fps_2,],sorted=FALSE)))
  
  #simmA_fpSim <<- simmA_fpSim
  #simmA_fpSim_vec <- as.vector(simmA_fpSim)
  
  #simmA_fpSim <<- simmA_fpSim
  #lincs_fpsim_2 <- lincs_fps_2[lincs_fps_2]
  #lincs_fpsim_2 <- lincs_fps_2[lincs_fps_2]
  #fpSim_end_time <- Sys.time()
  #total_fpSim_time[[j]] <- fpSim_end_time - fpSim_start_time
  #print(paste("Total time to run fpSim: ", total_fpSim_time[[j]], sep=""))
  
  #lincs_fp_ind <- lsm_rows[lincs_fpsim]
  #lincs_vec <- rownames(lincs_fps_2)[lincs_fp_ind]
  #simmA_fpSim_df <- data.frame(lincs_vec,simmA_fpSim)
  #return(simmA_fpSim_df)
  
  #lincs_fp_ind <- lsm_rows[lincs_fpsim]
  lincs_fp_ind <- as.vector(lincs_fpsim)
  lincs_vec <- rownames(lincs_fps_2)[lincs_fp_ind]
  simmA_fpSim_df <- data.frame(lincs_vec,simmA_fpSim)
  return(simmA_fpSim_df)
  
}