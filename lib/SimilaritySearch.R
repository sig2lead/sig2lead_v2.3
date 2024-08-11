similarity_search <- function(LSM_ID){
lincs <- read.csv(file="./similarity_search/lincs_cluster_centroid.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE)  
lincs2nci <- read.csv(file="./similarity_search/lincs2nci_dataBase.csv", sep=",", stringsAsFactors = FALSE)

#############################################################
# Search Lincs Column 1 for LSM_ID and Get Centroid Column 3
############################################################

lincs_centroid <- lincs[which(lincs[,1] == LSM_ID),3]

#################################################################################
# Search for Centroid in lincs2nci file Column 1 and Output Matching NCI Compunds
#################################################################################

nci_compounds <- lincs2nci[which(lincs2nci[,1] == lincs_centroid), 2]

library(bazar)
#ifelse(!is.empty(nci_compounds), return(nci_compounds), return(0))
output_nci <- list()
ifelse(!is.empty(nci_compounds),
       output_nci <- nci_compounds, 
       output_nci <- 0)
output_nci <<- output_nci
#return(nci_compounds)                           
}
  

#lincs <- read.csv(file="./similarity_search/lincs_cluster_centroid.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE)  
#lincs2nci <- read.csv(file="./similarity_search/lincs2nci_dataBase.csv", sep=",", stringsAsFactors = FALSE)

#lincs_centroid <- lincs[which(lincs[,1] == "LSM-5741"),3]
#nci_compounds <- lincs2nci[which(lincs2nci[,1] == lincs_centroid), 2]
#library(bazar)
#ifelse(!is.empty(nci_compounds), return(nci_compounds), return(0))
