#####################################
# StitchNetFn v9.0                  #
# Last update: Feb 14, 2019         #
# Author: Somchai Chutipongtanate   #
#####################################

library(visNetwork)
library(XML)
library(RCurl)
library(bitops)
library(scrapeR)
library(readr)
library(dplyr)

# Creating stitchNet function
  stitchNet <- function(chem = 'NULL', 
                        gene = 'NULL', 
                        confidence = '400', 
                        species = '9606',
                        limit = '10'){
    
# Disable this below code in order to run this function as part of Sig2Lead shinyApp
  chemInput <- chem
 
# Clean chemInput to solve Stitch API compatibility issues
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
                '&limit=', limitInt
  )
  print(api)
  
# Get interScore (edge) for visNetwork
  ####### Old script - a messy version of data collection and preprocessing for interScore object #####  
  # for the record: it's been 6 months since the first stitchNetFn finished (on Aug 2, 2018) - That was truly one of the best moment in my life
  ## this code chunk is parts of my precious memory, so this may be messy but always beautiful in my eyes

# Get interScore (edge) for visNetwork  
#  doc <<- readLines(api)
#  doc <- data.frame(id = doc)
#  newCol1 <- strsplit(as.character(doc$id), '\t',fixed=TRUE)
#  doc <- data.frame(doc, do.call(rbind, newCol1))
#  newCol2 <-strsplit(as.character(doc$X15),'|',fixed=TRUE)
#  doc <- data.frame(doc, do.call(rbind, newCol2))
#  newCol3 <- strsplit(as.character(doc$X1.1),':',fixed=TRUE)
#  doc <- data.frame(doc,do.call(rbind, newCol3))
#  #interScore <- data.frame(doc[, c(4, 5, 22)])
#  interScore <- data.frame(doc$X3, doc$X4, doc$X2.2)
#  names(interScore) <- c('from', 'to', 'score')
#  width <- (as.numeric(as.character(interScore$score))*2)^2
#  interScore <- data.frame(interScore, width)
#  names(interScore) <- c('from', 'to', 'score', 'width')
#  interScore <- data.frame(interScore, title = paste0("<p>", "Score:", "<br>", interScore$score, "</p>"))
#  interScore <<- interScore
  ##########################################################  
  doc <- suppressMessages(read_tsv(api, col_names = FALSE)) 
  tmp <- strsplit(doc$X15, "|", fixed = TRUE) %>% sapply("[", 1) %>% 
            strsplit(":", fixed = TRUE) %>% sapply("[", 2) %>% 
            as.numeric()
  interScore <- data.frame("from" = doc$X3, "to" = doc$X4, "score" = tmp, "width" = (tmp*2)^2, 
                           "title" = paste0("<p>", "Score:", "<br>", tmp, "</p>"), stringsAsFactors = FALSE)
  interScore <<- interScore
  
# Get molnodes (node) for visNetwork
  #### Old script: not that bad #############################
#  mol <- c(as.character(interScore$from), as.character(interScore$to))
#  mol <- unique(mol)
#  mol <<- mol
#  group <- replace(mol, toupper(mol) %in% toupper(chemInput), "Input chemicals present in STITCH")
#  group <- replace(group, toupper(group) %in% toupper(genes), "Gene KD")
#  group <- replace(group, !group %in% c("Input chemicals present in STITCH", "Gene KD"), "STITCH connectors")
#  molnodes <- data.frame(label = mol, id = mol, group = group, 
#                         title = paste0('<a target= "_blank"', '', 'href=', '', 
#                                        'https://en.wikipedia.org/wiki/', mol, '>', mol, '</a>')
#  )
  ###############################################
  mol <- unique(c(interScore$from, interScore$to)) 
  group <- toupper(mol) %>% replace(. %in% toupper(chemInput), "Candidates in STITCH") %>%
    #replace(. %in% toupper(genes), "Gene KD") %>%
    replace(. %in% toupper(genes), "Target Gene") %>%
    replace(!. %in% c("Candidates in STITCH", "Target Gene"), "STITCH gene/compound")

# Create landing pages for compounds and genes in molnodes
  link_modnodes <- character()
  for (i in seq_along(mol)) for (j in seq_along(chemInput)) {
    if(grepl(toupper(genes), toupper(mol[i])) ){
      link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=', mol[i], '>', mol[i], '</a>')
    } else if(grepl("CHEMBL", mol[i]) ){
      link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://www.ebi.ac.uk/chembl/beta/compound_report_card/', mol[i], '>', mol[i], '</a>')
    } else if(grepl("ZINC", mol[i]) ){
      link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'http://zinc15.docking.org/substances/',  mol[i], '>', mol[i], '</a>')
    } else if(grepl("LSM", mol[i]) ){
      link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'http://lincsportal.ccs.miami.edu/SmallMolecules/view/', mol[i], '>', mol[i], '</a>')
    } else if(toupper(mol[i]) %in% toupper(chemInput[j]) ){
      link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://pubchem.ncbi.nlm.nih.gov/search/#query=', mol[i], '>', mol[i], '</a>')
    } else if( is.na(link_modnodes[i]) ){
      link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://en.wikipedia.org/wiki/', mol[i], '>', mol[i], '</a>')
    }
  }

  molnodes <- data.frame(label = mol, id = mol, group = group, title = link_modnodes)
  
# get stitch undetected compounds (caused by whatever reasons i.e., name/id confusing between chemInput and stitch) and their landing pages 
  ind <- toupper(chemInput) %in% toupper(mol)
  stitch_undetect <- chemInput[!ind]
  stitch_undetect <- if(as.numeric(length(stitch_undetect))==0){"NULL"} else{stitch_undetect} 

# create landing pages for undetect_Node
  link_undetect_Node <- character()
  for (i in seq_along(stitch_undetect)){
    if(grepl("CHEMBL", stitch_undetect[i]) == 1){
      link_undetect_Node[i] <- paste0('<a target="_blank"href=https://www.ebi.ac.uk/chembl/beta/compound_report_card/', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
    } else if(grepl("ZINC", stitch_undetect[i]) == 1){
      link_undetect_Node[i] <- paste0('<a target="_blank"href=http://zinc15.docking.org/substances/', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
    } else if(grepl("LSM", stitch_undetect[i]) == 1){
      link_undetect_Node[i] <- paste0('<a target="_blank"href=http://lincsportal.ccs.miami.edu/SmallMolecules/view/', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
    } else {
      link_undetect_Node[i] <- paste0('<a target="_blank"href=https://pubchem.ncbi.nlm.nih.gov/search/#query=', stitch_undetect[i], '>', stitch_undetect[i], '</a>')
    }
  }
  
  undetect_Node <- data.frame(label = stitch_undetect, id = stitch_undetect, group = "Undetected by STITCH", title = link_undetect_Node) 

# Adding undetected compounds into molnodes
  molnodes <- rbind(molnodes, undetect_Node) %>% unique()
  molnodes <<- molnodes

# Generate a network!
  visNetwork(molnodes, interScore) %>%
      visNodes(physics = TRUE) %>%
      visEdges(physics= TRUE, smooth = TRUE, color = list(color = "#008080", highlight = "blue")) %>% 
      visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
                font = list(color = "660000"), 
                color = list(background = "ffa7a7", border = "black", 
                             highlight = list(background = "ffa7a7", 
                                              border = "red"))) %>%
      visGroups(groupname = "Gene KD", shape = "circle",
                font =  list(color = "blue"),
                color = list(background = "fdff00", border = "black",
                             highlight = list(background = "feff8f", 
                                              border = "red")) ) %>%
      visGroups(groupname = "STITCH gene/compound", shape = "dot", size = 10,
                font = list(color = "magenta", size = 10),
                color = list(background = "ddd" , border = "black", 
                             highlight = list(background = "caff44", 
                                              border = "red"))) %>%
      visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                color = list(background = "82a9ff", border = "253085",  
                             highlight = list(background = "4e91fc", 
                                              border = "red"))) %>%
      visIgraphLayout(randomSeed = 9, layout = "layout_with_mds", physics = TRUE, smooth = TRUE) %>% 
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
      visInteraction(navigationButtons = TRUE) %>%
      visLegend(width = 0.25) %>%
      visExport(type="png", name="export-network", float="left") # note: visExport only works under Shiny environment

  }

# Example
  # stitchNet(chem = c("nilotinib", "sunitinib", "curcumin", "MG-132", "P-5091"), gene = "mycn")
