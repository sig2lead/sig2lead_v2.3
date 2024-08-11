library(XML)
library(RCurl)
library(bitops)
library(scrapeR)
library(visNetwork)

stitchNet <- function(chem = 'NULL', 
                    gene = 'NULL', 
                    confidence = '400', 
                    species = '9606',
                    limit = '2'){

url_stitch <- 'http://stitch.embl.de/api/psi-mi-tab/interactionsList?'

#chemInput <<- chem
#chems <- chem 
chems <- paste(chem, collapse = '%0D')
#for multiple chems, by using c('chem1', 'chem2',..., 'chem_n')

genes <<- gene 
#genes <- paste(gene, collapse = '%0D')
# the same input in the front page 'Sig2Lead' app 
##but can be applied for multiple genes, by using c('gene1', 'gene2',..., 'gene_n')

score <- confidence
#score 400 = moderate confident (set as a default value)
##score 700 = strong confident
###score >900 = very strong confident
####score range = 0-999

speciesID <- species
#By default, species ID = 9606 = H.sapians
##Mus Musculus txid = 10090
###Zebra fish (Danio rerio) txid= 7955

limitInt <- limit
#By default, maximum numbers of node to return = 10 
## lower the limit, less the interaction  
### limit range = 0 - 25 is suggestive

api <- paste0(url_stitch, 'identifiers=', chems,'%0D', genes,
              '&required_score=', score, 
              '&species=', speciesID, 
              '&limit=', limitInt
              )

# convert HTMLInternalDocument to character vector
doc <<- readLines(api)

# separate xxxxx\txxxx\t to multiple columns
doc <- data.frame(id = doc)
newCol1 <- strsplit(as.character(doc$id),'\t',fixed=TRUE)
doc <- data.frame(doc,do.call(rbind, newCol1))

# separate xxx|xxx|xxx to multiple columns
newCol2 <-strsplit(as.character(doc$X15),'|',fixed=TRUE)
doc <- data.frame(doc,do.call(rbind, newCol2))

# separate xxx:xxx to 2 columns 
newCol3 <- strsplit(as.character(doc$X1.1),':',fixed=TRUE)
doc <- data.frame(doc,do.call(rbind, newCol3))

# get molecular pairs with interactive scores by excluding unwanted data
## get edge for visNetwork!
interScore <- doc[ , -c(1:3, 6:(ncol(doc)-1))]

# change column names
names(interScore) <- c("from", "to", "value")
interScore <<- interScore

# pull all molecular names and convert to a vector
mol <- c(as.character(interScore$from), as.character(interScore$to))

# get unique names and labeling the input gene
mol <- unique(mol)
#group <- as.integer(mol==toupper(genes))
#shape <- as.integer(mol==chemInput)
#shape <- replace(shape, shape==0, 'circle')
#shape <- replace(shape, shape==1, 'triangle')
color <- as.integer(mol==toupper(genes))
color <- replace(color, color==1, '#fff400')
color <- replace(color, color==0, '#97C2FC')


# node labeling
## get node for visNetwork!
#molnodes <<- data.frame(label = mol, id = mol, group = group, shape = 'circle')
molnodes <<- data.frame(label = mol, id = mol, color = color, shape = 'circle')

# labeling the input gene
#group <- as.integer(molnodes$id==toupper(genes)) 
#molnodes <<- cbind(molnodes, group)
#}
# Generate a network!
#visNetwork(molnodes, interScore) %>% 
#  visEdges(color = "#efc0f2") %>%
#  visLayout(randomSeed = 1999) %>%
#  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
#  visInteraction(navigationButtons = TRUE)

}
