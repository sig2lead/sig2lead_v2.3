library(enrichR)
library(visNetwork)
library(XML)
library(RCurl)
library(bitops)
library(scrapeR)
library(readr)
library(dplyr)

mStitch <- function(chem = 'NULL', 
                    gene = 'NULL', 
                    confidence = '400', 
                    species = '9606',
                    limit = '10'){
  
  chemInput <- chem
  chemInput_clean <- chemInput
  for (i in 1:length(chemInput_clean)){
    if(grepl("\\(\\(", chemInput_clean[i])){ 
      print(paste("Format not recognized by STITCH:", " ",chemInput_clean[grep("\\(\\(", chemInput_clean[i])]))
      chemInput_clean <- chemInput_clean[!grepl("\\(\\(", chemInput_clean)] 
    } else if(grepl(",", chemInput_clean[i])){
      print(paste("Format not recognized by STITCH:", " ", chemInput_clean[i]))
      chemInput_clean <-  chemInput_clean[!grepl(",", chemInput_clean)]
    }
  }
  
  # Get API url to pull the chemical-target interaction data from STITCH (stitch.embl.de)
  url_stitch <- 'http://stitch.embl.de/api/psi-mi-tab/interactionsList?'  
  
  chems <- paste(chemInput_clean, collapse = '%0D')
  chems <- gsub(" ", "%20", chems)
  genes <<- gene 
  score <- confidence # score 400 = moderate confident (set as a default value), score 700 = strong confident, score >900 = very strong confident, score range = 0-999
  speciesID <- '9606' #By default, species ID = 9606 = H.sapians; Mus Musculus txid = 10090, Zebra fish (Danio rerio) txid= 7955
  limitInt <- (as.numeric(limit) - 1) # By default, maximum numbers of connection of the input nodes = 10; a range of 0 - 25 is suggestive
  
  api <- paste0(url_stitch, 'identifiers=', chems,'%0D', genes,
                '&required_score=', score, 
                '&species=', speciesID, 
                '&limit=', limitInt)
  
  doc <- suppressMessages(read_tsv(api, col_names = FALSE)) 
}
