library(XML)
library(RCurl)
library(bitops)
library(scrapeR)
library(visNetwork)

# Creating the stitchNet function
stitchNet <- function(chem = 'NULL', 
                    gene = 'NULL', 
                    confidence = '400', 
                    species = '9606',
                    limit = '10'){

# Get API url to pull the chemical-target interaction data from STITCH (stitch.embl.de)
url_stitch <- 'http://stitch.embl.de/api/psi-mi-tab/interactionsList?'
chemInput <<- chem
chems <- paste(chemInput, collapse = '%0D')
genes <<- gene 
score <- '400' #score 400 = moderate confident (set as a default value), score 700 = strong confident, score >900 = very strong confident, score range = 0-999
speciesID <- '9606' #By default, species ID = 9606 = H.sapians; Mus Musculus txid = 10090, Zebra fish (Danio rerio) txid= 7955
limitInt <- '10' #By default, maximum numbers of connected nodes= 10; a range of 0 - 25 is suggestive

api <- paste0(url_stitch, 'identifiers=', chems,'%0D', genes,
              '&required_score=', score, 
              '&species=', speciesID, 
              '&limit=', limitInt
              )

# Get interScore (edge) for visNetwork
doc <<- readLines(api)
doc <- data.frame(id = doc)
newCol1 <- strsplit(as.character(doc$id), '\t',fixed=TRUE)
doc <- data.frame(doc, do.call(rbind, newCol1))
newCol2 <-strsplit(as.character(doc$X15),'|',fixed=TRUE)
doc <- data.frame(doc, do.call(rbind, newCol2))
newCol3 <- strsplit(as.character(doc$X1.1),':',fixed=TRUE)
doc <- data.frame(doc,do.call(rbind, newCol3))
#interScore <- data.frame(doc[, c(4, 5, 22)])
interScore <- data.frame(doc$X3, doc$X4, doc$X2.2)
names(interScore) <- c('from', 'to', 'score')
width <- (as.numeric(as.character(interScore$score))*2)^2
interScore <- data.frame(interScore, width)
names(interScore) <- c('from', 'to', 'score', 'width')
interScore <- data.frame(interScore, title = paste0("<p>", "Score:", "<br>", interScore$score, "</p>"))
interScore <<- interScore

# Get molnodes (node) for visNetwork
mol <- c(as.character(interScore$from), as.character(interScore$to))
mol <- unique(mol)
mol <<- mol
group <- replace(mol, toupper(mol) %in% toupper(chemInput), "Input chemicals present in STITCH")
group <- replace(group, toupper(group) %in% toupper(genes), "Gene KD")
group <- replace(group, !group %in% c("Input chemicals present in STITCH", "Gene KD"), "STITCH connectors")
molnodes <- data.frame(label = mol, id = mol, group = group, 
                       title = paste0('<a target= "_blank"', '', 'href=', '', 
                                      'https://en.wikipedia.org/wiki/', mol, '>', mol, '</a>')
                      )

## get information of LSM compounds (undetected by stitch) by hyperlink to iLINCS
inputChem <- toupper(chemInput) %in% toupper(mol)
inputChem <<- inputChem
stitch_undetect <- chemInput[!inputChem]

## When stitch_undetect produces any missing value, replacing missing value with NULL
stitch_undetect <- if(as.numeric(length(stitch_undetect))==0){"NULL"} else{stitch_undetect} 
undetect_Node <- data.frame(label = stitch_undetect, id = stitch_undetect, group = "Undetected by STITCH", 
                            title = paste0('<a target= "_blank"', '', 'href=', '', 
                                    'http://lincsportal.ccs.miami.edu/SmallMolecules/view/', 
                                    stitch_undetect, '>', stitch_undetect, '</a>')
                            )

## Adding undetected compounds into molnodes
molnodes <- rbind(molnodes, undetect_Node)
molnodes <<- molnodes

# Generate a network!
visNetwork(molnodes, interScore) %>%
  visNodes(physics = FALSE) %>%
  visEdges(physics= FALSE, smooth = FALSE, color = list(color = "gray", highlight = "blue")) %>%
  visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
            color = list(background = "ffa7a7", border = "black", 
                         highlight = list(background = "ffa7a7", 
                                          border = "red"))) %>%
  visGroups(groupname = "Gene KD", shape = "circle",
            color = list(background = "fdff00", border = "black",
                         highlight = list(background = "feff8f", 
                                          border = "red")) ) %>%
  visGroups(groupname = "STITCH connectors", shape = "dot", size = 10,
            color = list(background = "ddd", border = "black", 
                         highlight = list(background = "caff44", 
                                          border = "red"))) %>%
  visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
            color = list(background = "4e91fc", border = "253085", 
                         highlight = list(background = "4e91fc", 
                                          border = "red"))) %>%
  visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = FALSE) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
  visInteraction(navigationButtons = TRUE) %>%
  visLegend(width = 0.25) %>%
  #visExport only works under Shiny environment, so enable a code below when running in Shiny environment
  visExport(type="png", name="export-network", float="left") 

}
