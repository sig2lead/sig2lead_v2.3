######################################
# Min-Sim
#####################################

#fp_test <- fpset_added[1,]
#fp_test <- input[2, c(3:1026)]
#test_minSim <- minSim(fp_test, lsm_rows, m )

minSim <- function(query, lsm_rows, m_test){
  
  #load("./lib/minsim_objects.RData")
  #input <- read.table(file="C:\\Data\\LINCS_NCI\\phiClust_input_label_tab.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")

  sim <- vector("integer", length=42764)
  tanimoto_plusplus <- vector("numeric", length=42764)

  #lsm_labels <- input[seq(from=1,to=85527,by=2),1]
  #lsm_labels <- gsub(">","",lsm_labels)
  
  query_1 <- which(query == 1)
  m_query <- length(query_1)
  
  # query_1 <- which(query[] == 1)
  # m_query <- sum(query_1)
  # 
  #query_labels <- rownames(query_mat)
  #lsm_labels <- rownames(lsm_mat)


  for(j in 1:length(query_1)){
    #sim[sets_minority_cols[[query_1[]]][l]] <- sim[sets_minority_cols[[j]][l]] + 1
    #sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
    sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
  }

  for(p in 1:42764){
    tanimoto_plusplus = sim / (m_test + m_query - sim)}

  #return(unlist(tanimoto_plusplus[[lsm_rows]]))

  return(tanimoto_plusplus[lsm_rows])
  
  
}