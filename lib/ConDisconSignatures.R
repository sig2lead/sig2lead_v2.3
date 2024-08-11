########################################
# Find Con/Dis/cordant Signtures
########################################

ConDisconSignatures <- function(signatureID)
{
    #signatureID <- sigid
    #input$ILINCS <- "Current"
  
    con_comma <<- list("character")
    #print(input$ConOrDiscon)
    for (i in 1:length(signatureID))
    {
      #try
      if(input$ILINCS == "Legacy"){
          #url2 <- paste("http://www.ilincs.org/api/SignatureMeta/findConcordantSignatures?sigID=%22",signatureID,"%22&lib=%22LIB_5%22", sep="")
          url2 <- (paste("http://ilincs2018.ilincs.org/ilincs/api/SignatureMeta/findConcordantSignatures?sigID=%22",signatureID[[i]],"%22&lib=%22LIB_5%22", sep=""))
          df <- fromJSON(url2)
      }
      #try
      else if (input$ILINCS == "Current"){
          #signatureID <- sigid
          url2 <- paste("http://www.ilincs.org/api/SignatureMeta/findConcordantSignatures?sigID=%22",signatureID[i],"%22&lib=%22LIB_5%22", sep="")
          #url2 <- (paste("http://ilincs2018.ilincs.org/ilincs/api/SignatureMeta/findConcordantSignatures?sigID=%22",signatureID[[i]],"%22&lib=%22LIB_5%22", sep=""))
          req <- GET(url = url2)
          shiny::validate(need(req$status_code == 200, "Your query returned a server error."))
          #df <-fromJSON(httr::content(req,type="text", encoding="UTF-8"))
          json_test_con <- httr::content(req,type="text", encoding="UTF-8")
          shiny::validate(need(jsonlite::validate(json_test_con), "The output from your query is in the wrong format."))
          
          df <-fromJSON(httr::content(req,type="text", encoding="UTF-8"))
      }
    
      if (length(df)>0)
          df2<<- df[order(df$similarity), ]
    
      if(input$ConOrDiscon == "Inhibit"){
          con <<- df2[df2$similarity > input$Concordance,]
          #con <<- df2[df2$similarity > .2,]
          #discon <- df2[df2$similarity < -0.235, ]
          discon <- df2[df2$similarity < input$Concordance, ]
      }
      else if(input$ConOrDiscon == "Activate"){
          con <<- df2[df2$similarity < input$Concordance,]
          #discon <- df2[df2$similarity < -0.235, ]
          discon <- df2[df2$similarity < input$Concordance, ]
      }
    
      concordant_list <- list("character")
      for (j in 1:length(con))
      {
          concordant_list[j] <- con[j]
      }
    
      con_comma[[i]] <- concordant_list
    
    }
  
  con_comma <<- con_comma
}

