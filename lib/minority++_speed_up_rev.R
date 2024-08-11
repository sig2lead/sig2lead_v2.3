########################################
# Tanimoto++
########################################

########################################
# Create lists of minority state vector 
########################################

m <- fp_data_3$Sum

fp_data_5 <- fp_data_3[,1:1024]

save(fp_data_3, sets_minority_cols, file="./minsim_objects.RData")
write.csv(fp_data_3, file="./lincs_fps_w_m.csv")

sets_minority_cols <- list()

for(i in 1:1024){
  sets_minority_cols[[i]] <- which(fp_data_5[,i] == 1)
  print(paste("Iteration: ", i))
}

increment <- function(x)
{
  eval.parent(substitute(x <= x + 1))
}

total_time_minoritylsh <- list()
#library(iterators)
tanimoto_plusplus <- list()

for(p in 1000:2000){
query_1 <- which(fp_data_5[p,] == 1)
length_query1 <- length(query_1)
sim <- vector("integer", length=42764)


#start_minoritylsh <- proc.time()
start_minoritylsh <- Sys.time()
 # for(j in 1:length_query1){
 #   for(l in 1:length(sets_minority_cols[[j]])){
 #     #sim[sets_minority_cols[[query_1[]]][l]] <- sim[sets_minority_cols[[j]][l]] + 1
 #     sim[sets_minority_cols[[query_1[j]]][l]] <- sim[sets_minority_cols[[query_1[j]]][l]] + 1
 #   }        
 # }

for(j in 1:length_query1){
    #sim[sets_minority_cols[[query_1[]]][l]] <- sim[sets_minority_cols[[j]][l]] + 1
    #sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
    sim[sets_minority_cols[[query_1[j]]]] <- sim[sets_minority_cols[[query_1[j]]]] + 1
}

# #for(j in 1:length_query1){
# 
 # for(k in 1:42764){
 # 
 #   tanimoto_plusplus[k] <- (sim[k]) / (m[p] + m[k] - sim[k])
 #  #   #print(paste("Iteration: ", k))
 # }
  # 
tanimoto_plusplus[[p]] = sim / (m[p] + m - sim)
  
#end_minoritylsh <- proc.time()
end_minoritylsh <- Sys.time()
#total_time_minoritylsh <- end_minoritylsh - start_minoritylsh
total_time_minoritylsh[[p]] <- end_minoritylsh - start_minoritylsh
}


total_time <- mean(unlist(total_time_minoritylsh))

# list_times <- vector("numeric", length=100)
# 
 for (z in 1000:2000){
   list_times[z] <- total_time_minoritylsh[[z]][3]
# }
 mean_time <-mean(list_times)

list_times <- mean(unlist(total_time_minoritylsh))

neighbors_9 <- which(tanimoto_plusplus > 0.90)


sim_real <- simil(x=input[2,c(3:1026)], y=input[6,c(3:1026)], method="Rogers", by_rows=TRUE, auto_convert_data_frames = TRUE, convert_distances=FALSE)

a1 = (paste(input[2,c(3:1026)], sep="",collapse=" "))
b1 = paste(input[6,c(3:1026)], sep="", collapse=" ")

l1 <- list(a1,b1)
sim_real <- 
  dist_real <- 1 - sim_real

x1 <- vector("integer", length=10)
a <- c(0,0,1,1,0,1,0,1,1,1,1,0)
b <- c(0,0,0,0,0,1,0,1,1,1,0,0)
l <- list(a,b)
sim_real <- simil(l, method="Jaccard", by_rows = TRUE)
test_jacc <-
  
  test <- rbind(input[2,c(3:1026)], input[4,c(3:1026)])
test_2 <- rbind(a,b)
sim_real_tan <- simil(test_2, method="Tanimoto", by_rows=TRUE, auto_convert_data_frames = TRUE, convert_distances=FALSE)  
#######################################################################################################
###################################################
# Calculate Exact Tanimoto Similarity/Distance 
# For Any Two Points
###################################################


library(proxy)

proxy_total_3 <- list()
sim_real_test_3 <- list()

for (z in 1000:2000){
  proxy_start <- Sys.time()
  sim_real_test_3[[z]] <- simil(x=input[(z*2),c(3:1026)], y=input[seq(from=2,to=85528,by=2), c(3:1026)], method="Jaccard", by_rows=TRUE, auto_convert_data_frames = TRUE, convert_distances=FALSE)
  proxy_end <- Sys.time()
  
  proxy_total_3[[z]] <- proxy_end - proxy_start
}

mean_proxy_3 <- mean(unlist(proxy_total_3))
neighbors_9 