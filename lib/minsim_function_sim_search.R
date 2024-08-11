

minsim_sim_search <- function(fpset_added,lsm_rows,number_analogs){
  load("./minSim_apfp_RObjects.RData")
  load("./lincs_fps.RData")
  tanimoto_final <- list()
  lincs_final <- list()
  m <- m_lincs
  # 
  # begin_sim_calc <- Sys.time()
  for(k in 1:dim(fpset_added)[1]){
  #   
    query_1 <- which(fpset_added[k,] == 1)
    sim <- vector("integer", length=41572)
    m_query <- length(query_1)

    #if(m_query == 0){
    #  #print(paste("k", k))
    #  next
    #}
    
    if(m_query == 0){
      #print(paste("k", k))
      tanimoto_plusplus <- NA
      tanimoto_final[[k]] <- NA
      lincs_final[[k]] <- NA
      next
    }
    

    for(j in 1:length(query_1)){

      sim[minSim_list[[query_1[j]]]] <- sim[minSim_list[[query_1[j]]]] + 1
    }

  #   #tanimoto_plusplus[[k]] = sim / (m + m_query - sim)
     tanimoto_plusplus = sim / (m + m_query - sim)
  #   #print(paste("Tanimotoplusplus: ",tanimoto_plusplus))
  #   
  #   #tanimoto_final[[k]] <- tanimoto_plusplus[[k]][tanimoto_plusplus[[k]] > 0.8]
    
     
     ####################### Edited 02/15/2021 ####################################
     #tanimoto_final[[k]] <- max(tanimoto_plusplus[lsm_rows])
     #lincs_final[[k]] <- which.max(tanimoto_plusplus[lsm_rows])
     
     #tanimoto_final[[k]] <- slice_max(tanimoto_plusplus[lsm_rows], n=2, with_ties=FALSE)
     lincs_final[[k]] <- which.maxn(tanimoto_plusplus[lsm_rows], n=number_analogs)
     tanimoto_final[[k]] <- tanimoto_plusplus[lincs_final[[k]]]
     
     ##############################################################################
  #   #print(paste("Max Tanimotoplusplus: ", max(tanimoto_plusplus)))
  #   #tanimoto_final[[k]] <- max(tanimoto_plusplus[[k]])
  #   
  #   
  #   #simMA <- cbind(simMA, tanimoto_final[[k]])  
  #   #simMA <- cbind(simMA, tanimoto_final[[k]])
  #   #print(paste("Finished Compound: ", k))
  }
  
  # 
  # end_sim_calc <- Sys.time()
  # #colnames(simMA) <- sdfid(subset_add)
  # #simMA <<- simMA
  # tanimoto_final <<- tanimoto_final
  # #simMA <<- fpSim(x=lsm_compounds, y=added_compounds, method="Tanimoto")
  # #end_sim_calc <- Sys.time()
  # total_sim_time <<- end_sim_calc - begin_sim_calc
  # print(paste("Total Time to Calculate Similarity Search: ", total_sim_time, sep= " "))
  #tanimoto_w_lincs <- data.frame(unlist(lincs_final), unlist(tanimoto_final))
  #return(tanimoto_w_lincs)
  
  #################### Edited 02/15/2021 ###############################################
  #lincs_fp_ind <- lsm_rows[unlist(lincs_final)]
  #lincs_vec <- rownames(lincs_fps_2)[lincs_fp_ind]
  #minSim_df <- data.frame(lincs_vec,unlist(tanimoto_final))
  
  lincs_fp_ind <- lsm_rows[unlist(lincs_final)]
  lincs_vec <- rownames(lincs_fps_2)[lincs_fp_ind]
  minSim_df <- data.frame(lincs_vec,unlist(tanimoto_final))
  
  lincs_final <<- lincs_final
  tanimoto_final <<- tanimoto_final
  lincs_fp_ind <<- lincs_fp_ind
  lincs_vec <<- lincs_vec
  minSim_df <<- minSim_df
  #######################################################################################
  
 
  return(minSim_df)
  
  # 
}