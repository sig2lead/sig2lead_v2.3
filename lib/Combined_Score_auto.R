combined_score <- function(){
  
  #######Extract Vectors of Similarites for Added Compounds from Similarity Matrix
  
      #adds_csv
  
#      if(!(grepl(".sdf", adds$name))){
#        added_sims <- simMA[(rownames(simMA) %in% adds_csv[,2]),]
#      }
#      else if(grepl(".sdf", adds$name)){
        #added_ids <<- unlist(lapply(1:length(sdfset_add),function(x) sdfset_add@SDF[[x]]@header[1]))
        #added_sims <<- simMA[(rownames(simMA) %in% added_ids),]
        #added_ids_vec <- paste(added_ids,collapse="|")
        #added_sims <<- simMA[which(grepl(added_ids_vec, rownames(simMA))),]
#      }
      
  
  #######Extract Concordances for LSM-IDs in Similarity Matrix
      
      
  
  ########Calculate Concordance Matrix###############
  # Added Compounds X LINCS Compounds
  ##################################################
      
      # added_sims for LSM x Con for LSM
      added_sims <<- simMA
      #combined_score <-  max_cons_2[which(max_cons_2$`LSM-ID`%in% colnames(added_sims)),]
    
      combined_score <-  max_cons_4[which(max_cons_4$`LSM-ID`%in% rownames(added_sims)),]
    
      #added_sims_2 <- added_sims[,(colnames(added_sims) %in% combined_score[,1])]
      
      #added_sims_o <- added_sims_2[,order(colnames(added_sims_2))]
      #added_sims_o <- added_sims[,order(colnames(added_sims))]
      #added_sims_t <<- t(added_sims_o)
      #print(paste("Dimensions of added_sims_t:", dim(added_sims_t), sep= " "))
            
      
  ####################################
  # Order Concordances
  ####################################
      
      combined_score_o <<- combined_score[order(as.character(combined_score[,1])),]
      concordance <<- combined_score_o
      
      
  ################################
  # Multiply concordance x Similarity
  ################################    
      
  score_df <-  matrix(0, nrow=nrow(added_sims), ncol=ncol(added_sims))
  rownames(score_df) <- rownames(added_sims)
  colnames(score_df) <- colnames(added_sims)
  concordance_vector <<- as.vector(concordance[,4])
  
  concordance_df <- concordance[, c(1,4)]
  
  print("Made it this far...")
#  for (i in 1:nrow(added_sims_o)){    
#      for (j in 1:length(concordance[,1])){
#      score_df[i,j] <- added_sims_o[i,j] * concordance[j,3]
#      }
#    print(paste("Compound Score Interation: ", i, sep=""))
#  }
  
  print("Going to run vectorized solution for combined score")
  start_combined_score <<- Sys.time()
  
  #score_df <- added_sims * concordance_vector
  #score_df <- added_sims * concordance_vector
  
#-032919------------------------------------------------------------------------------  
  added_sims_df <- as.data.frame(added_sims)
  added_sims_df$'LSM-ID' <- rownames(added_sims_df)
  score_df2 <- added_sims_df
  library(dplyr)
  score_df2 <- left_join(added_sims_df, concordance_df, by = "LSM-ID")
  
score_df3 <- score_df2[ , 1:10]
 rownames(score_df3) <- score_df2$`LSM-ID`
 score_df4 <- score_df3 * as.numeric(score_df2$Concordance)
 score_df4$'LSM-ID' <- score_df2$`LSM-ID`
 
 score_df <- score_df4[, 1:10]
 
 x <- colnames(score_df)
 mycompound <- as.character(colnames(score_df))
 
 
 z <- list()
 for(i in 1:ncol(score_df)){
 ind[i] <- which(score_df[,i] == max(score_df[,i]), arr.ind = TRUE)
 z[[i]] <- c(mycompound[i], as.character(rownames(score_df[ind[i], ])),  score_df2$Concordance[ind[i]],
             score_df2[ind[i], i], score_df[ind[i], i])
 }
 
 
 z_df <<- data.frame(matrix(unlist(z), nrow = length(z), byrow =T)) 
 
 

 
 
 
#----------------------------------------------------------------------------------------- 
  
  end_combined_score <<- Sys.time()
  print(" Done calculating combined score")
  
  total_combined_score_time <<- end_combined_score - start_combined_score
  print(paste("Time to calculate combined score: ", total_combined_score_time, sep=""))
  #####################
#  combinations  <<- as.vector(score_df) 
#  added_sims_vec <<- as.vector(added_sims)
  #concordance_vec <<- rep(concordance[,3], each=(length(combinations)/nrow(concordance)))
#  concordance_vec <<- rep(concordance_vector, times=(length(combinations)/nrow(added_sims)))
  #added_sims_vec <- as.vector(as.matrix(added_sims_o))
  
  #times_row <- nrow
#  rows_rep <<- rep(rownames(score_df), times=(length(combinations)/nrow(score_df)))
#  cols_rep <<- rep(colnames(score_df), each = length(combinations)/ncol(score_df))

   
#  ranked_scores <- score_df
   #ranked_scores <- data.frame(rank(score_df)) 
   #return(score_df)
#  scores_combined <- as.data.frame(cbind(cols_rep, rows_rep, concordance_vec, added_sims_vec, combinations)) 
  #scores_combined <- as.data.frame(cbind(rows_rep, cols_rep, combinations))
#  scores_combined$combinations <- as.numeric(as.character(scores_combined$combinations))
#  ranked_scores <- scores_combined[order(scores_combined[,5], decreasing = TRUE),]
#  colnames(ranked_scores) <- c("My Compounds", "LINCS Compounds", "Concordance", "Similarity", "Combined Score")
#  ranked_scores_2 <<- ranked_scores[!duplicated(ranked_scores$'My Compounds'),]
  end_time <<- Sys.time()
#  ranked_scores_2[,-c(1,2)] <- apply(ranked_scores_2[,-c(1,2)], 2, as.numeric)
#  ranked_scores_2[,-c(1,2)] <- round(ranked_scores_2[,-c(1,2)],3)
#  return(ranked_scores_2) 
 
return(z_df) 
}

