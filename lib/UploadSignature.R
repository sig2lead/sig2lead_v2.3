
sigUpload <- function(filename){
  #Do stuff here
  #signature <- read.table(file = paste("%22", filename, "%22", sep = ""))
  #ftemp = tempfile(pattern="filename", fileext=".xls", tmpdir = tempdir())
  #write.csv(filename, ftemp)
  print(paste("filename: ", filename))
  ############# Added 08/08/2024
  #rv <- reactiveValues(r = NULL)
  
  #observe({
  #  req(input$UploadSignature)
  #  rv$data <- read.csv(input$inFile$datapath)

  #r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(filename)))
  #})
  #################################
  ##################################################################################################################
  # Before Edit
  
  #r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(filename)))
  #r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file("TestSig1.txt")))
  
  r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(filename)), timeout(30))
  
  #r <<- r
  
  rjson <- toJSON(content(r)$status$concordanceTable)
  
  rjson <- fromJSON(rjson)
  
  rjson_2 <<- rjson
  
  #rjson <- rjson
  #rjson_table <- (data.frame(matrix(unlist(rjson), nrow=length(rjson), byrow=T)))
  #colnames(rjson_table) <- c("similarity", "pValue", "nGenes", "compound", "lincsPertID", "genes", "concentration", "time", "_row", "signatureid", "cellline")
  
  #rjson_table <<- rjson_table
  
  connected_sigs <- rjson[,c(4,5,10,1)]
  
  #connected_sigs <- rjson
  
  colnames(connected_sigs) <- c("Candidate Name","Candidate LSM ID", "Cell Line", "Concordance")
  
  connected_sigs$Concordance <- round(as.numeric(connected_sigs$Concordance),4)
  connected_sigs <<- connected_sigs
  ###################################################################################################################
  # After Edit
  #signatureFile  <- ("sample2.csv")
  #signatureFile$ID_geneid <- as.character(signatureFile$ID_geneid)
  #apiUrl <- "http://www.ilincs.org/api/ilincsR/findConcordances"
  #req <- (POST(apiUrl, body = list(file=signatureFile, lib="LIB_5"), encode = "form"))
  #req <- (POST(apiUrl, body = list(file=signatureFile, lib="LIB_5"), encode = "form"))
  #req <- (POST(apiUrl, body = list(file=signatureFile, lib="LIB_5"), encode="json"))
  #req <- (POST(apiUrl, body = list(file=(signatureFile)), encode="json"))
  #rjson <- data.table::rbindlist(httr::content(req)$concordanceTable, use.names = TRUE, fill = TRUE)
  #rjson <- toJSON(content(req)$status$concordanceTable)
  #rjson <<- fromJSON(rjson)
  
  #req <- POST("http://www.ilincs.org/api/ilincsR/findConcordances?lib=LIB_5", body = list(file = upload_file("sample2.txt")))
  #req <- POST("http://www.ilincs.org/api/ilincsR/findConcordances", body = list(lib="LIB_5", file = upload_file("./sample2.txt")))
  #http://www.ilincs.org/api/ilincsR/findConcordances?lib=LIB_5&file=upload_file("./sample2.txt")
  
  #head(output)
  
  #####Down the road status$data will likely become something else (status$concordanceTable)
  #l <- lapply(content(r)$status$data, function(x)unlist(x))
  #l <- lapply(content(r)$status$concordanceTable, function(x)unlist(x))
  #####
  
  
  #ilincsresult <<- data.frame(t(sapply(l, c)), stringsAsFactors = FALSE)
  #ilincsresult <- data.frame(fromJSON(l))
  #colnames(ilincsresult) <- c("similarity", "pValue", "nGenes", "compound", "lincsPertID", "concentration", "time", "_row", "signatureid", "cellline")
  #ilincsresult <<- ilincsresult
  #print(content(r))
  #r$status_code
}

