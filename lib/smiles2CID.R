library(XML)
library(magick)

smiles2CID <- function(smiles = 'NULL'){
# only the canonical smiles is accepted, since the isomeric smiles contain "/" backslash which is an issue in R

# generate API call
url <- 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'
smilesInput <- smiles
path <- '/cids/xml'
#path <- '/record/xml'
pathImg <- '/PNG'

# get API
pubChemAPI <- paste0(url, smilesInput, path)
pubChemAPIimg <- paste0(url, smilesInput, pathImg)

# read API XML response as character vector
cid <- readLines(pubChemAPI)

# get CID
cid <- gsub("<CID>", "", cid[7])
cid <- gsub("</CID>", "", cid)
cid <- gsub("  ", "CID", cid)
CID <<- cid
print(CID)

# also get a chemical structure
image_read(pubChemAPIimg)

}
