####################################
# User Input Variable (Gene KD)
####################################
#library(httr)
#library(jsonlite)
#library(RJSONIO)
#detach("package:RJSONIO", unload=TRUE)

define_knockdown <- function(gene)
{
  usr <- gene
  
  #usr <- "BCL2A1"
  if(input$ILINCS == "Legacy"){
    
      url1 <- paste("http://ilincs2018.ilincs.org/ilincs/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%2C%20%22libraryid%22%3A%22LIB_6%22%7D%7D", sep="")
      #url1 <- paste("http://ilincs.org/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%2C%20%22libraryid%22%3A%22LIB_6%22%7D%7D", sep="")
  
      #try - if status code != 200 OR
      
      raw.result <- GET(url = url1)
      #req <- GET(url = url1)
      
      shiny::validate(need(raw.result$status_code == 200, "Your query returned a server error."))
    
      sigid <<- list("character")
  
      if(length(content(raw.result)) != 0){
        for (i in 1:length(content(raw.result)))
        {
          sigid[i] <- content(raw.result)[[i]]$signatureid
      #catch    
      
        }
      } else{
        print("A knockdown signature for youe gene target is not in LINCS.  Would you like to try a different gene target?")
        sigid <- NULL
        }
  
    #sigid <<- sigid
  }
  else if(input$ILINCS == "Current"){
    
      #url1 <- paste("http://www.ilincs.org/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%7D%7D", sep="")
      #usr <- 'bcl2a1'
      
      url1 <- paste("http://www.ilincs.org/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%2C%20%22libraryid%22%3A%22LIB_6%22%7D%7D",sep="")
    
      #try status code != 200 
      req <- GET(url = url1)
      shiny::validate(need(req$status_code == 200, "Your query returned a server error."))
      json_test <- httr::content(req,type="text", encoding="UTF-8")
      shiny::validate(need(jsonlite::validate(json_test), "The output from your query is in the wrong format."))
      sigid <-fromJSON(httr::content(req,type="text", encoding="UTF-8"))$signatureid
    
      if (length(sigid) == 0){
        sigid <<- NULL
        print("A knockdown signature for youe gene target is not in LINCS.  Would you like to try a different gene target?")
        #return()
      }
      else {
        sigid <<- sigid
      }
    #catch
  }
  #return(sigid)
  sigid <<- sigid
}
