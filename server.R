#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#




#########Load all required libraries here, they cannot be sourced
library(shiny)
library(shinyjs)
library(httr)
library(jsonlite)
library(DT)
#library(heatmaply)
#library(shinyHeatmaply)
library(visNetwork)
library(ChemmineOB)
library(ChemmineR)
library(plyr)
library("dendextend")
library("colorspace")
library(ggforce)
library(rlist)
library(scatterpie)
library(ggrepel)
library(bazar)
library(XML)
library(RCurl)
library(bitops)
#library(scrapeR)
library(igraph)
library(circlize)
library(dplyr)
library(doBy)
#library(ComplexHeatmap)

#source("http://bioconductor.org/biocLite.R")

#biocLite("ChemmineOB")
#biocLite("ChemmineR")

load("./minSim_apfp_RObjects.RData")
source("lib/chemmineR_2_option.R")
source("./lib/cluster_fpsim.R")
######Discordant or Concordant
######Benchmarking -- Behrouz

begin_time <<- Sys.time()
options(shiny.maxRequestSize=1000*1024^2)
#options(shiny.maxRequestSize=200*1024^2)

#Added 8/8/2024###############
#session$onSessionEnded(stopApp)
###############################
shinyServer(function(input, output, session) {
  #session$reload()
  #output$mytest <- renderUI(
  #  if(!is.null(input$AddedCompounds)) {
  #    div(textInput("AddedLabel", "Label for Added Compounds", value="My Candidates"),
  #        checkboxInput("Bypass", "Bypass Clustering", value=TRUE)
  #    )
  #  }
  #)

 output$logo <- renderImage({
   list(src='./www/sig2lead_rev_old.png',
        position='right')
 }, deleteFile=FALSE)
 
 output$logo2 <- renderImage({
   list(src='./www/sig2lead_rev.png',
        position='right')
 }, deleteFile=FALSE)
 

 output$wait_message <- renderText({ "Please note that concurrent use of the server by multiple users\ncan cause delays in processing your job. Consider resubmitting your\njob if you do not see the results after several minutes." })
 
 output$ExampleFile <- downloadHandler(
   filename <- function() {
     paste("example_candidates_sdf", "sdf", sep=".")
   },
   
   content <- function(filename) {
     file.copy("./www/EGFR.sdf", filename)
   }
   
 )
 output$ExampleFileSmiles <- downloadHandler(
   filename <- function() {
     paste("example_candidates_smiles", "smi", sep=".")
   },
   
   content <- function(filename) {
     file.copy("./www/example_candidates_smiles.smi", filename)
   }
   
 )
 output$ExampleSig <- downloadHandler(
   filename <- function() {
     paste("Example_Signature", "txt", sep=".")
   },
   
   content <- function(filename) {
     file.copy("./www/Example_Signature.txt", filename)
   }
   
 )
 

 ###############################################################################################################################################
 ###################### Signature Connectivity Analysis (Tab 1) ################################################################################
 ###############################################################################################################################################
 
 index_analogs <<- 1


###########################################################################################
 observeEvent(input$Go,
  
    {
     
      rm(list=ls())
      adds <- input$AddedCompounds
      print(input$AddedCompounds$name)
      shinyjs::hideElement("CANDIDATES")
      shinyjs::hideElement("SMILES")
      shinyjs::hideElement(id="logo")
      shinyjs::hideElement(id="logo2")
     
      Bypass = TRUE 
      
      withProgress(message = "Running Sig2Lead", value = 0, {
        
       
        
        source("lib/GeneKD.R", local = TRUE)
        source("lib/ConDisconSignatures.R", local = TRUE)
        
        
        lincs_compounds <- read.csv("./LINCSCompounds.csv", stringsAsFactors = FALSE)
        #lsm_rows <- 1:41572
      #1. Query iLINCS for user input gene
        if (input$Signature == "Define Target Gene"){
          
          shinyjs::hideElement("check_kd")
          shinyjs::hideElement("check_connection")
          
         
          
          ##########################################################################################################
          # Query iLINCS and Extract Gene KD Signatures
          #############################################
          define_knockdown(input$gene_knockdown)
                           
          if (length(sigid) == 0){
            shinyjs::hideElement("logo")
            #shinyjs::hideElement("check_connection")
            shinyjs::showElement("check_kd")
            output$check_kd <- renderUI({
            p("Please enter another gene target. A knockdown signature for that target is not in LINCS.")
            })
            return()
          }
          print("Found your knockdowns") 
          if(input$ConOrDiscon=="Concordant"){
          incProgress(1/7, message = "Finding Concordant Signatures")
         
          
          ###########################################################################################################  
          # Query iLINCS and Extract Small Molecule Signatures Concordant/Discordant to Gene KD Signatures
          ################################################################################################
          ConDisconSignatures(sigid)}
          
          ConDisconSignatures <<- ConDisconSignatures(sigid)
          
          if (length(con_comma) == 0){
            shinyjs::hideElement("logo")
            shinyjs::showElement("check_sm")
            output$check_sm <- renderUI({
              p("Your query returned no concordant/discordant small molecules.")
            })
            return()
          }
          
          ##############################################################
          ##############################################################
          # Redo with LINCS fingerprints loaded
          load("./lincs_fps.RData")
          #lincs_fps_2 <- read.csv("./lincs_fps.csv")
          #write.csv(lincs_fps_2, file="./lincs_fps.csv")
          
          source("lib/GetSMILES.R", local = TRUE)
          #3. Generate SMILES of compounds
             GetSMILES()
             #print("SMILES acquired")
             #incProgress(1/7, message = "Converting SMILES to SDF...Slow")
          # #4. Generate Data Table
             display <<- final_compounds
             colnames(display) <- c("LSM_ID", "SMILES")
          #   
            
          #   
            makeLink <- function(val) {
               #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
               paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
             }
         
          # Need to tweak with LINCS fingerprints pre-loaded
           df <- list()
           for (d in 1:length(con_comma)){
           
             df[[d]] <- as.data.frame(matrix(unlist(con_comma[[d]]), nrow=length(con_comma[[d]][[1]])))
           
           }
           df_concordances <<- ldply(df, data.frame)
           #df_concordances_2 <- df_concordances[,c(11,6,2,3)]
           
           if(input$ILINCS == "Legacy"){
              df_concordances_2 <- df_concordances[,c(11,6,8,2,3)]
           }
           else if (input$ILINCS == "Current"){
             df_concordances_2 <- df_concordances[,c(13,7,9,2,3)]
           }
           colnames(df_concordances_2) <- c("LSM-ID", "Compound", "Cell Line", "Concordance", "Significance")
           df_concordances_2$Concordance <- as.numeric(as.character(df_concordances_2$Concordance))
           df_concordances_o <<- df_concordances_2[order(df_concordances_2$Concordance, decreasing = TRUE ),]
          # 
          # 
           df_concordances_no_dups <<- df_concordances_o[!duplicated(df_concordances_o["LSM-ID"]),]
           display_no_dups <- display[!duplicated(display["LSM_ID"]),]
          # 
           df_concordances_o_smiles <- merge(df_concordances_no_dups, display_no_dups, by.x="LSM-ID", by.y="LSM_ID", all.x=FALSE, all.y=FALSE)
           df_concordances_o_smiles<<- df_concordances_o_smiles[order(df_concordances_o_smiles$Concordance, decreasing = TRUE ), c(1,5,2,3,4)]
          # 
          # #############Output to LINCS Compounds Table;
          # 
          # #################
          # # Max Concordance
          # #################
           max_cons <- df_concordances_o_smiles[order(df_concordances_o$`LSM-ID`, -abs(df_concordances_o$Concordance)), ]
           max_cons_2 <- max_cons[!duplicated(max_cons$`LSM-ID`),]
           max_cons_3 <<- max_cons_2[order(max_cons_2$Concordance, decreasing = TRUE), -5]
           max_cons_4 <<- max_cons_3[-which(is.na(max_cons_3$`LSM-ID`)),]
           con_df <- as.data.frame((matrix(unlist(con_comma),byrow=T)))
           max_cons_4$Concordance <- round(max_cons_4$Concordance,3)
           lincs_output <- max_cons_4
           colnames(lincs_output) <- c("Candidate LSM ID", "LINCS Candidate Name", "Cell Line", "Concordance", "Significance")
           lincs_output <- lincs_output[,-5]
           lincs_output <- lincs_output[,c(2,1,3,4)]
           shinyjs::hideElement("sim_search")
           shinyjs::hideElement("logo")
           shinyjs::hide("logo")
           shinyjs::showElement("SMILES")
           
           output$SMILES <<- DT::renderDataTable({
             lincs_output$`Candidate LSM ID`<- makeLink(lincs_output$`Candidate LSM ID`)
             lincs_output
             return(lincs_output)}, caption="LINCS Candidates Ranked", escape=FALSE, rownames=FALSE)
           lincs_output <<- lincs_output
           load("./lincs_fps.RData")
           lincs_max <<- max(lincs_output$Concordance)
          # 
           lsm_rows <<- which(rownames(lincs_fps_2) %in% max_cons_4$`LSM-ID`)
        
           display <<- display
          # 
           source("lib/ChemmineOB.R", local = TRUE)
          # 
          # ###This is the slow step, see if we can go straight from smiles to fingerprints and seperately convert to SDF
           #print("Converting SMILES to SDF")
           Get_SDF()
          index <<- 1
          index_analogs <<- 2
         
        }  
        else if (input$Signature == "Upload a Signature"){
          #source("lib/UploadSignature.R")
          
          ##### Added 8/8/2024
          
          d <- tempdir(check = FALSE)
          baseuuid <- paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")
          uuid <- paste(
            substr(baseuuid,1,8),
            "-",
            substr(baseuuid,9,12),
            "-",
            "4",
            substr(baseuuid,13,15),
            "-",
            sample(c("8","9","a","b"),1),
            substr(baseuuid,16,18),
            "-",
            substr(baseuuid,19,30),
            sep="",
            collapse=""
          )
          destfile <- file.path(d,uuid)
          file.rename(from=input$UploadSignature$datapath,to=destfile)
          
          r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(destfile)), timeout(30))
          
          
          print("Made it to before post call")
         
          print("Made it past Post")
          print(r)
          #rjson <- content(r)
          rjson <- parsed_content(r)
          rjson_p <- rjson$status$concordanceTable
          #rjson <- do.call(rbind, rjson)
          rjson_2 <- do.call(rbind, lapply(rjson_p, data.frame))
          #rjson <- toJSON(content(r)$status$concordanceTable)
          print(rjson_2)
          rjson_2 <<- rjson_2
          print("After first rjson")
          #rjson <- fromJSON(rjson)
          print("After second rjson")
          
          #rjson_2 <<- rjson
          connected_sigs <- rjson_2[,c(4,5,11,1)]
          
          #connected_sigs <- rjson
          
          colnames(connected_sigs) <- c("Candidate Name","Candidate LSM ID", "Cell Line", "Concordance")
          
          connected_sigs$Concordance <- round(as.numeric(connected_sigs$Concordance),4)
          connected_sigs <<- connected_sigs
          
          ##############################################################
          #sigUpload(input$UploadSignature$datapath)
          print(input$UploadSignature$datapath)
          print("Found your signature")
          incProgress(1/7, message = "Finding Concordant Signatures")
          shinyjs::hideElement("sim_search")
          shinyjs::hideElement("logo")
          shinyjs::showElement("SMILES")
          
          #########################################################################################
          # Added 08/13/2022
          
          
          if(input$ConOrDiscon_2 == "Concordant"){
            # For SM Signatures
            con <<- connected_sigs[connected_sigs$Concordance > input$Concordance_2,]
            
            connected_sigs_2 <- con
          }
          else if(input$ConOrDiscon_2 == "Discordant"){
            con <<- connected_sigs[connected_sigs$Concordance < input$Concordance_3,]
            
            connected_sigs_2 <- con
          }
          #########################################################################################
          shinyjs::hideElement("sim_search")
          shinyjs::hideElement("logo")
          shinyjs::showElement("SMILES")
          
          output$SMILES <<- DT::renderDataTable({
            connected_sigs_2$`Candidate LSM ID`<- makeLink5(connected_sigs_2$`Candidate LSM ID`)
            connected_sigs_2
            return(connected_sigs_2)}, caption="LINCS Candidates Ranked", escape=FALSE, rownames=FALSE)
          
          load("./lincs_fps.RData")
          
          lsm_rows <<- which(rownames(lincs_fps_2) %in% connected_sigs_2$`Candidate LSM ID`)
          display <- connected_sigs_2[,c(2,1)]
          
          ###########11/28/2023
          display <- as.data.frame(display)
          rownames(display) <- NULL
          display[,1] <- as.character(display[,1])
          display[,2] <- as.character(display[,2])
          display <<- display
          ########################
          max_cons_4 <- connected_sigs_2
          colnames(max_cons_4) <- c("Compound", "LSM-ID", "Cell Line", "Concordance")
          max_cons_4$'Compound' <- as.character(max_cons_4$'Compound')
          max_cons_4$'LSM-ID' <- as.character(max_cons_4$'LSM-ID')
          max_cons_4$'Cell Line' <- as.character(max_cons_4$'Cell Line')
          
          max_cons_4 <- max_cons_4[!duplicated(max_cons_4$`LSM-ID`),]
          
          ######## Added
          
          max_cons_4 <<- max_cons_4
          lincs_output <- connected_sigs_2
          lincs_output$ID <- 1:nrow(lincs_output)
          lincs_output_2 <- apply(lincs_output,2,as.character)
          lincs_output_2 <<- lincs_output_2
          ####### Added 04/17/2023####################
          lincs_max <<- max(lincs_output$Concordance)
          index <<- 2
          index_analogs <<- 2
          #############################################
          
          
          ##################################################################
          ############# Added 11/27/2023 ###################################
          
          
          #index <<- 2
          #index_analogs <<- 2
          
          makeLink5 <- function(val) {
           
            paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
          }
          
         
          if (input$Signature == "Upload a Signature"){
            output$downloadUpload <- renderUI({
              if((input$Signature == "Upload a Signature")){
                downloadButton("downloadSigs", label = "LINCS Signatures Ranked")}
            })
          }
          output$downloadSigs <- downloadHandler(
            filename = function() {
              'lincs_signatures_ranked.csv'
            },
            content = function(filename){
              write.csv(lincs_output_2[,-5], filename)
            }
          )
          }
        
        
        else if(input$Signature == "Find Analogs in LINCS")
        {
          lsm_rows <- 1:41572
          load('./lincs_fps.RData')
          
          ####### Added 11/15/2023####################
          #lincs_max <<- max(lincs_output$Concordance)
          index <<- 2
          index_analogs <<- 2
          ############################################# 
          
          
          adds<<-input$AddedCompounds
          
          if (is.null(input$AddedCompounds)){
            shinyjs::hideElement("logo")
            output$sim_search <- renderUI({
              p("No added compounds.  Similarity search cannot be performed.")
              
            })
            return()
          }
          if (!grepl(".sdf|.smi",adds$name)){
            shinyjs::hideElement("logo")
            output$check_format <- renderUI({
              p("File in Wrong Format.  Please input sdf or smi file.")
              
            })
            return()
          }
          else if(!(grepl(".sdf", adds$name))){
            adds_SMI <<- read.SMIset(adds$datapath)
            
            sdfset_add <<- smiles2sdf(adds_SMI)
            
            ####### Added 11/15/2023####################
            #lincs_max <<- max(lincs_output$Concordance)
            index <<- 2
            index_analogs <<- 2
            display <<- NULL
            ############################################# 
            
          }
          else if(grepl(".sdf", adds$name)){
            print("Added compounds in sdf format")
            max_cons_4 <<- NULL
            topNum <<- input$topN
            
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ####### Added 11/15/2023####################
            #lincs_max <<- max(lincs_output$Concordance)
            index <<- 2
            index_analogs <<- 2
            display <<- NULL
            sdfset_add <- read.SDFset(adds$datapath)
            added_labels <- sdfid(sdfset_add)
            
            ####################################
            # ChemmineR Atom Pair Fingerprints
            ##################################
            apset_added <-sdf2ap(sdfset_add)
            fpset_added <- as.matrix(desc2fp(apset_added))
            
            sdfset_add <<- as.matrix(fpset_added)
            
            ############################################# 
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          }
        }
        
        
        
      #2. Get a list of all concordant compounds
        
        
    
       makeLink5 <- function(val) {
           #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
           #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"</a>", sep="")
           paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
         }
     
      #5. Cluster compounds identified through LINCS, generate heatmaps
       

    if(!is.null(input$AddedCompounds)){
      adds<<-input$AddedCompounds
      if (!grepl(".sdf|.smi",adds$name)){
        shinyjs::hideElement("logo")
        output$check_format <- renderUI({
          p("File in Wrong Format.  Please input sdf or smi file.")
          
        })
        return()
      }
      
      #5.1. Include Added Compounds if applicable
        
         if (Bypass == FALSE)      {                     
                      #adds<<-input$AddedCompounds
                      
                      if(!(grepl(".sdf", adds$name))){
                        adds_SMI <<- read.SMIset(adds$datapath)
                       
                        sdfset_add <<- smiles2sdf(adds_SMI)
                        
                        cluster_compounds()
                        print("Your compounds have been clustered with LINCS compounds")
                        source("lib/ColorMap.R", local = TRUE)
                        color_dend()
                      }
                      else if(grepl(".sdf", adds$name)){
                        print("Added compounds in sdf format")
                        
                        lsm_smiles_2 <- lsm_smiles
                        sdfset_add <<- read.SDFset(adds$datapath)
                        added_smiles <<- sdf2smiles(sdfset_add)
                        added_labels <- cid(added_smiles)
                        added_smiles_2 <- as.character(added_smiles[1:length(added_smiles)])
                        added_smiles_3 <- cbind(added_labels, added_smiles_2)
                        
                       
                        colnames(added_smiles_3) <- c("Compound_ID", "SMILES")
                        
                        lsm_smiles_2_labels <- cid(lsm_smiles_2)
                        lsm_smiles_3 <- as.character(lsm_smiles_2[1:length(lsm_smiles_2)])
                        lsm_smiles_3b <- unname(lsm_smiles_3)
                        lsm_smiles_4 <- cbind(lsm_smiles_2_labels, lsm_smiles_3b)
                        colnames(lsm_smiles_4) <- c("Compound_ID", "SMILES")
                        colnames(display) <-  c("Compound_ID", "SMILES")
                        
                        test2 <<- rbind(added_smiles_3, display)
                        
                        colnames(test2) <<- c("Compound_ID", "SMILES")
                        
                        
                        sdf_smiles <<- c(sdf_smiles, sdfset_add)
                        cluster_compounds()
                        print("Your compounds have been clustered with LINCS compounds")
                        source("lib/ColorMap.R", local = TRUE)
                        adds_SMI <<- added_smiles
                        color_dend()
                        }
                      
                }
      #5.2. Don't add compounds
      #5.2. Add compounds, but bypass clustering
      
        #########################################################################################################################################################
        ############################# Compute Tanimoto of User-Provided Candidates and LINCS Analogs ############################################################
        #########################################################################################################################################################
        else
                {
                ############################################################################################
                ##### User-Provided Compounds Not in SDF Format
                ############################################################################################
                  adds<<-input$AddedCompounds
                  if(!(grepl(".sdf", adds$name))){
                    adds_SMI <<- read.SMIset(adds$datapath)
                    #adds_SMI <- sdfset_add_smiles
                    #adds_SMI<-read.SMIset(paste(getwd(), adds$name, sep="/"))
                    #adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
                    #adds_csv <- adds_csv[c(2,1)]
                    #colnames(adds_csv) <- c("LSM_ID", "SMILES")
                    #adds_csv <<- adds_csv
                    #test2 <<- rbind(adds_csv, display)
                    #colnames(test2) <<- c("Compound_ID", "SMILES")
                    
                    print("In smiles before sdfset_add")
                    #print(adds_SMI)
                    cid(adds_SMI) <- makeUnique(cid(adds_SMI))
                    sdfset_add <- smiles2sdf(adds_SMI)
                    
                    #cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                    sdfset_add <<- sdfset_add
                    print("in smiles after sdfset_add")
                    #sdf_smiles <<- c(sdf_smiles, sdfset_add)
                    #bypass_clustering(sdf_smiles, sdfset_add)}
                    
                    
                    ############### FPsim and Define Target Gene #################################
                    
                    if ((input$Algorithm == "fpSim") & (input$Signature == "Define Target Gene")){
                      #adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      # 
                      # if (cid(sdfset_add) == ""){
                      #   cid(sdfset_add) <- "Compound"
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      # else{
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      
                      
                      ifelse( (cid(sdfset_add) == ""),{
                        cid(sdfset_add) <- "Compound"
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      },
                      {
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      })
                      
                      
                      #cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      added_labels <- cid(sdfset_add)
                      #added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      
                      #Edited 4-9-2021
                      #sums <- rowSums(sdfset_add)
                      #zeroes <- which(sums==0)
                      #sdfset_add <- sdfset_add[-zeroes,]
                      #added_labels <- added_labels[-zeroes]
                      
                      test_bypass_fpsim <-  bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      test_bypass_fpsim <- cbind(added_labels, test_bypass_fpsim)
                      colnames(test_bypass_fpsim) <- c("Compound", "LSM-ID", "Similarity")
                      
                      Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_fpsim))
                      
                      ######## Changed 11-18-2020
                      test_bypass_fpsim <- left_join(test_bypass_fpsim, max_cons_4, by = "LSM-ID")
                      test_bypass_fpsim <- test_bypass_fpsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      
                      
                      #for (i in 1:nrow(test_bypass_fpsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_fpsim$Concordance <- Concordance
                      #test_bypass_fpsim$Cell_Line <- Cell_Line
                      #test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      #test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      #test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance", "Similarity")]
                      #colnames(test_bypass_fpsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_fpsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                      
                      
                      for (i in 1:nrow(test_bypass_fpsim)){
                        test_bypass_fpsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_fpsim[i,2]), 1]
                      }
                      
                      test_bypass_fpsim <- test_bypass_fpsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_fpsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      for(i in 1:nrow(test_bypass_fpsim)){
                      #i <- 1
                       test_bypass_fpsim$sstar[i]  <- -(log10(1-(test_bypass_fpsim$Similarity[i] - 0.0001)) + log10(1-((test_bypass_fpsim$Concordance[i] - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_fpsim)[7] <- "S*"
                      test_bypass_fpsim$'S*'<- as.numeric(test_bypass_fpsim$'S*')
                      test_bypass_fpsim$'S*' <- round(test_bypass_fpsim$'S*', 3)


                      
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_fpsim$`LINCS Analog`<- makeLink5(test_bypass_fpsim$`LINCS Analog`)
                        #test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                        test_bypass_fpsim
                        return(test_bypass_fpsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      test_bypass_fpsim <<- test_bypass_fpsim
                      
                      
                      
                      
                      # Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      # for (i in 1:nrow(test_bypass_fpsim)){
                      #   Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      # }
                      # test_bypass_fpsim$Concordance <- Concordance
                      # test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      # test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      # 
                      # test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      # output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,4], test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      # test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                      #test_bypass_fpsim <<- test_bypass_fpsim
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
                      #stopApp()
                    }
                    
                    #else minsim
                    
                    ################### MinSim and Define Target Gene #######################################
                    
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Define Target Gene"))
                    { #load("./minSim_apfp_RObjects.RData")
                      #adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      print("In minsim and define target gene")
                      #added_labels <- sdfid(sdfset_add)
                      
                      # if (cid(sdfset_add) == ""){
                      #   cid(sdfset_add) <- "Compound"
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      # else{
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      
                      ifelse ((cid(sdfset_add) == ""),{
                        cid(sdfset_add) <- "Compound"
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      },{
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      })
                        #print("cid sfset_add", cid(sdfset_add))
                      
                      #cid(sdfset_add) <- "Test_Label"
                      added_labels <- cid(sdfset_add)
                      
                      #test_added_labels <<- added_labels
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      
                      apset_added <-sdf2ap(sdfset_add)
                      #test_apset_added_1 <<- apset_added
                      
                      #if (length(sdfset_add) > 1){
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      #}
                      #else{
                      #  print("Before fpset_added")
                      #  test_apset_added_2 <<- apset_added[1]
                      #  fpset_added <- desc2fp(apset_added)
                      #  print("After fpset_added")
                      #  
                      #  sdfset_add <- fpset_added
                        
                      #  print("after sdfset_add")
                      #}
                      print(dim(sdfset_add))
                      
            
                      test_bypass_minsim <- bypass_clustering_minsim(sdfset_add, lsm_rows,max_cons_4)
                      #test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim
                      
                      
                     
                      test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim)
                      
                      # Edited 03_25_2021
                      test_bypass_minsim <- test_bypass_minsim[which(!(is.na(test_bypass_minsim[,2]))),]
                      
                      
                      colnames(test_bypass_minsim) <- c("Compound", "LSM-ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_minsim))
                      
                      test_bypass_minsim <- left_join(test_bypass_minsim, max_cons_4, by = "LSM-ID")
                      test_bypass_minsim <- test_bypass_minsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance
                      #test_bypass_minsim$Cell_Line <- Cell_Line
                      #test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      #test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance", "Similarity")]
                      #colnames(test_bypass_minsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_minsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      #test_bypass_minsim$'LINCS Analog' <- makeLink(test_bypass_minsim$'LINCS Analog')
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),]
                      
                      
                      for (i in 1:nrow(test_bypass_minsim)){
                        test_bypass_minsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_minsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_minsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      
                      for(i in 1:nrow(test_bypass_minsim)){
                      #i <- 1
                       test_bypass_minsim$sstar[i]  <- -(log10(1-(test_bypass_minsim$Similarity[i] - 0.0001)) + log10(1-((test_bypass_minsim$Concordance[i] - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_minsim)[7] <- "S*"
                      test_bypass_minsim$'S*'<- as.numeric(test_bypass_minsim$'S*')
                      test_bypass_minsim$'S*' <- round(test_bypass_minsim$'S*', 3)

                      
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_minsim$`LINCS Analog`<- makeLink5(test_bypass_minsim$`LINCS Analog`)
                        #test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),]
                        test_bypass_minsim
                        return(test_bypass_minsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      
                      test_bypass_minsim <<- test_bypass_minsim
                      #stopApp()
                    }
                    ############################################################################################
                    ####### Added 11/27/2023 ###################################################################
                    if ((input$Algorithm == "fpSim") & (input$Signature == "Upload a Signature")){
                      #adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      # 
                      # if (cid(sdfset_add) == ""){
                      #   cid(sdfset_add) <- "Compound"
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      # else{
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      
                      
                      ifelse( (cid(sdfset_add) == ""),{
                        cid(sdfset_add) <- "Compound"
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      },
                      {
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      })
                      
                      
                      #cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      added_labels <- cid(sdfset_add)
                      #added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      
                      #Edited 4-9-2021
                      #sums <- rowSums(sdfset_add)
                      #zeroes <- which(sums==0)
                      #sdfset_add <- sdfset_add[-zeroes,]
                      #added_labels <- added_labels[-zeroes]
                      
                      test_bypass_fpsim <-  bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      test_bypass_fpsim <- cbind(added_labels, test_bypass_fpsim)
                      colnames(test_bypass_fpsim) <- c("Compound", "LSM-ID", "Similarity")
                      
                      Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_fpsim))
                      
                      ######## Changed 11-18-2020
                      test_bypass_fpsim <- left_join(test_bypass_fpsim, max_cons_4, by = "LSM-ID")
                      test_bypass_fpsim <- test_bypass_fpsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      
                      
                      #for (i in 1:nrow(test_bypass_fpsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_fpsim$Concordance <- Concordance
                      #test_bypass_fpsim$Cell_Line <- Cell_Line
                      #test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      #test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      #test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance", "Similarity")]
                      #colnames(test_bypass_fpsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_fpsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                      
                      
                      for (i in 1:nrow(test_bypass_fpsim)){
                        test_bypass_fpsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_fpsim[i,2]), 1]
                      }
                      
                      test_bypass_fpsim <- test_bypass_fpsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_fpsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      for(i in 1:nrow(test_bypass_fpsim)){
                        #i <- 1
                        test_bypass_fpsim$sstar[i]  <- -(log10(1-(test_bypass_fpsim$Similarity[i] - 0.0001)) + log10(1-((test_bypass_fpsim$Concordance[i] - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_fpsim)[7] <- "S*"
                      test_bypass_fpsim$'S*'<- as.numeric(test_bypass_fpsim$'S*')
                      test_bypass_fpsim$'S*' <- round(test_bypass_fpsim$'S*', 3)
                      
                      
                      
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_fpsim$`LINCS Analog`<- makeLink5(test_bypass_fpsim$`LINCS Analog`)
                        #test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                        test_bypass_fpsim
                        return(test_bypass_fpsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      test_bypass_fpsim <<- test_bypass_fpsim
                      
                      
                      
                      
                      # Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      # for (i in 1:nrow(test_bypass_fpsim)){
                      #   Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      # }
                      # test_bypass_fpsim$Concordance <- Concordance
                      # test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      # test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      # 
                      # test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      # output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,4], test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      # test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                      #test_bypass_fpsim <<- test_bypass_fpsim
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
                      #stopApp()
                    }
                    
                    #else minsim
                    
                    
                    
                    
                    
                    
                    ################## Minsim and Upload a Signature (SMILES)###################################
                    ################# MinSim and Define Target Gene #######################################
                    
                    
                    
                    
                    
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Upload a Signature"))
                    { #load("./minSim_apfp_RObjects.RData")
                      #adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      print("In minsim and Upload a Signature")
                      #added_labels <- sdfid(sdfset_add)
                      
                      # if (cid(sdfset_add) == ""){
                      #   cid(sdfset_add) <- "Compound"
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      # else{
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      
                      ifelse ((cid(sdfset_add) == ""),{
                        cid(sdfset_add) <- "Compound"
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      },{
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      })
                      #print("cid sfset_add", cid(sdfset_add))
                      
                      #cid(sdfset_add) <- "Test_Label"
                      added_labels <- cid(sdfset_add)
                      
                      #test_added_labels <<- added_labels
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      
                      apset_added <-sdf2ap(sdfset_add)
                      #test_apset_added_1 <<- apset_added
                      
                      #if (length(sdfset_add) > 1){
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      #}
                      #else{
                      #  print("Before fpset_added")
                      #  test_apset_added_2 <<- apset_added[1]
                      #  fpset_added <- desc2fp(apset_added)
                      #  print("After fpset_added")
                      #  
                      #  sdfset_add <- fpset_added
                      
                      #  print("after sdfset_add")
                      #}
                      print(dim(sdfset_add))
                      
                      
                      test_bypass_minsim <- bypass_clustering_minsim(sdfset_add, lsm_rows,max_cons_4)
                      #test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim
                      
                      
                      
                      test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim)
                      
                      # Edited 03_25_2021
                      test_bypass_minsim <- test_bypass_minsim[which(!(is.na(test_bypass_minsim[,2]))),]
                      
                      
                      colnames(test_bypass_minsim) <- c("Compound", "LSM-ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_minsim))
                      
                      test_bypass_minsim <- left_join(test_bypass_minsim, max_cons_4, by = "LSM-ID")
                      test_bypass_minsim <- test_bypass_minsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance
                      #test_bypass_minsim$Cell_Line <- Cell_Line
                      #test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      #test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance", "Similarity")]
                      #colnames(test_bypass_minsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_minsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      #test_bypass_minsim$'LINCS Analog' <- makeLink(test_bypass_minsim$'LINCS Analog')
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),]
                      
                      
                      for (i in 1:nrow(test_bypass_minsim)){
                        test_bypass_minsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_minsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_minsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      
                      for(i in 1:nrow(test_bypass_minsim)){
                        #i <- 1
                        test_bypass_minsim$sstar[i]  <- -(log10(1-(test_bypass_minsim$Similarity[i] - 0.0001)) + log10(1-((test_bypass_minsim$Concordance[i] - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_minsim)[7] <- "S*"
                      test_bypass_minsim$'S*'<- as.numeric(test_bypass_minsim$'S*')
                      test_bypass_minsim$'S*' <- round(test_bypass_minsim$'S*', 3)
                      
                      
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_minsim$`LINCS Analog`<- makeLink5(test_bypass_minsim$`LINCS Analog`)
                        #test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),]
                        test_bypass_minsim
                        return(test_bypass_minsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      
                      test_bypass_minsim <<- test_bypass_minsim
                      #stopApp()
                    }
                    
                    ###########################################################################################
                    ################
                    
                    
                    
                    
                    
                    
                    
                    ############################################################################################
                    ################### MinSim and Find Analogs in LINCS #######################################
                    
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Find Analogs in LINCS"))
                    { #load("./minSim_apfp_RObjects.RData")
                      #adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      #added_labels <- sdfid(sdfset_add)
                      index_analogs <<- 1
                      
                      ifelse ((cid(sdfset_add) == ""),{
                        cid(sdfset_add) <- "Compound"
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      },{cid(sdfset_add) <- makeUnique(cid(sdfset_add))})
                      
                      #cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      added_labels <- cid(sdfset_add)
                      
                      
                      #added_labels_2 <<- rep(added_labels, each=input$topN)
                      
                      lsm_rows <- 1:41572
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      
                      
                      #test_bypass_minsim <- bypass_clustering_minsim_sim_search(sdfset_add, lsm_rows)
                      #test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim)
                      added_labels_2 <<- rep(added_labels, each=input$topN)
                      test_bypass_minsim <- bypass_clustering_minsim_sim_search(sdfset_add, lsm_rows, input$topN)
                      test_bypass_minsim <- data.frame(added_labels_2, test_bypass_minsim)
                      colnames(test_bypass_minsim) <- c("User-added Compound", "LINCS Analog", "Similarity")
                      #test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      #test_bypass_minsim <- test_bypass_minsim[-which(test_bypass_minsim$Similarity < input$Similarity),]
                      
                      test_bypass_minsim <- test_bypass_minsim %>% filter(Similarity > input$Similarity)
                      
                      
                      #Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance 
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,1],-test_bypass_minsim[,3]),]
                      
                      for (i in 1:nrow(test_bypass_minsim)){
                        test_bypass_minsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_minsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,4,3)]
                      colnames(test_bypass_minsim) <- c("User-added Compound", "LINCS Analog", "LINCS Analog Name", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_minsim$`LINCS Analog`<- makeLink5(test_bypass_minsim$`LINCS Analog`)
                        #test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),]
                        return(test_bypass_minsim)}
                        , caption="Analogs in LINCS", rownames=FALSE, escape=FALSE)
                      #stopApp()
                      lsm_rows <<- which(rownames(lincs_fps_2) %in% test_bypass_minsim[,2])
                      test_bypass_minsim <<- test_bypass_minsim
                    }
                    
                    ################### FPSim and Find Analogs in LINCS #######################################
                    
                    else if ((input$Algorithm == "fpSim") & (input$Signature == "Find Analogs in LINCS"))
                    { #load("./minSim_apfp_RObjects.RData")
                      #adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      #added_labels <- sdfid(sdfset_add)
                      index_analogs <<- 1
                      lsm_rows <- 1:41572
                      
                      # if (cid(sdfset_add) == ""){
                      #   cid(sdfset_add) <- "Compound"
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      # else{
                      #   cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      # }
                      
                      ifelse ((cid(sdfset_add) == ""),{
                        cid(sdfset_add) <- "Compound"
                        cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      },{cid(sdfset_add) <- makeUnique(cid(sdfset_add))})
                      
                      #cid(sdfset_add) <- makeUnique(cid(sdfset_add))
                      added_labels <- cid(sdfset_add)
                      added_labels_2 <<- rep(added_labels, each=input$topN)
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      #test_bypass_fpsim <- bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      #test_bypass_fpsim <- data.frame(added_labels, test_bypass_fpsim)
                      
                      #Edited 4-9-2021
                      sums <- rowSums(sdfset_add)
                      zeroes <- which(sums==0)
                      sdfset_add <- sdfset_add[-zeroes,]
                      added_labels_2 <- added_labels_2[-zeroes]
                      
                      
                      print(dim(sdfset_add))
                      test_bypass_fpsim <- bypass_clustering_fpsim_sim_search(sdfset_add, lsm_rows, input$topN)
                      test_bypass_fpsim <- data.frame(added_labels_2, test_bypass_fpsim)
                      colnames(test_bypass_fpsim) <- c("User-added Compound", "LINCS Analog", "Similarity")
                      #test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      #test_bypass_fpsim <- test_bypass_fpsim[-which(test_bypass_fpsim$Similarity < input$Similarity),]
                      test_bypass_fpsim <- test_bypass_fpsim %>% filter(Similarity > input$Similarity)
                      #Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      #for (i in 1:nrow(test_bypass_fpsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      #}
                      #test_bypass_fpsim$Concordance <- Concordance 
                      #test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      #test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,3], decreasing=TRUE),]
                      #test_bypass_fpsim <- test_bypass_fpsim_sim_search[order(test_bypass_fpsim[,3], decreasing=TRUE),]
                      test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,1],-test_bypass_fpsim[,3]),]
                      
                      for (i in 1:nrow(test_bypass_fpsim)){
                        test_bypass_fpsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_fpsim[i,2]), 1]
                      }
                      
                      test_bypass_fpsim <- test_bypass_fpsim[,c(1,2,4,3)]
                      colnames(test_bypass_fpsim) <- c("User-added Compound", "LINCS Analog", "LINCS Analog Name", "Similarity")
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_fpsim$`LINCS Analog`<- makeLink5(test_bypass_fpsim$`LINCS Analog`)
                        #test_bypass_fpsim[order(test_bypass_fpsim[,3], decreasing=TRUE),]
                        return(test_bypass_fpsim)}
                        , caption="Analogs in LINCS", rownames=FALSE, escape=FALSE)
                      test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                      lsm_rows <<- which(rownames(lincs_fps_2) %in% test_bypass_fpsim[,2])
                      test_bypass_fpsim <<- test_bypass_fpsim
                    }
                  }
                  
                  ################################################################################
                  #### USer-Provided Candidates in SDF Format ####################################
                  ################################################################################
                  
                  else if(grepl(".sdf", adds$name)){
                    print("Added compounds in sdf format")
                    # inFile <- input$AddedCompounds
                    # old_name <- inFile$datapath
                    # dirstr <- dirname(inFile$datapath)
                    # new_name <- paste(dirstr, inFile$name, sep="/")
                    # file.rename(old_name,new_name)
                    #sdfset_add <<- read.SDFset(adds$name)
                    #print(dim(sdfset_add))
                    #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                    #added_labels <- sdfid(sdfset_add)
                    
                    ####################################
                    # ChemmineR Atom Pair Fingerprints
                    ##################################
                    #apset_added <-sdf2ap(sdfset_add)
                    #fpset_added <- as.matrix(desc2fp(apset_added))
                    
                    #sdfset_add <<- as.matrix(fpset_added)
                    #print(dim(sdfset_add))
                    #sdf_smiles <<- c(sdf_smiles, sdfset_add[1:10])
                    
                    
                    # If fpsim <- UI widget == fpsim
                    
                    ################### FPSim and Define Target Gene #######################################
                    
                    if ((input$Algorithm == "fpSim") & (input$Signature == "Define Target Gene")){
                      adds<<-input$AddedCompounds
                      
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      sdfset_add <- read.SDFset(adds$datapath)
                      added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      
                      test_fpset <<- sdfset_add
                      print(dim(sdfset_add))
                      
                      #Edited 4-9-2021
                      #sums <- rowSums(sdfset_add)
                      #zeroes <- which(sums==0)
                      #sdfset_add <- sdfset_add[-zeroes,]
                      #added_labels <- added_labels[-zeroes]
                      added_labels <<- added_labels
                      lsm_rows <<- lsm_rows
            
                      test_bypass_fpsim <-  bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      print("Made it this far.")
                      
                      test_bypass_fpsim <- cbind(added_labels, test_bypass_fpsim)
                      
                      #edited 4-9-2021
                      test_bypass_fpsim <- test_bypass_fpsim[which(!(is.na(test_bypass_fpsim[,2]))),]
                      
                      colnames(test_bypass_fpsim) <- c("Compound", "LSM-ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_fpsim))
                      
                      test_bypass_fpsim <- left_join(test_bypass_fpsim, max_cons_4, by = "LSM-ID")
                      test_bypass_fpsim <- test_bypass_fpsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      
                      #for (i in 1:nrow(test_bypass_fpsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_fpsim$Concordance <- Concordance
                      #test_bypass_fpsim$Cell_Line <- Cell_Line
                      #test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      #test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      #test_bypass_fpsim <- test_bypass_fpsim[,c("User-added Compound", "LINCS Analog", "Cell_Line", "Concordance", "Similarity")]
                      #test_bypass_fpsim <-test_bypass_fpsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance", "Similarity")]
                      #colnames(test_bypass_fpsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_fpsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                      
                      for (i in 1:nrow(test_bypass_fpsim)){
                        test_bypass_fpsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_fpsim[i,2]), 1]
                      }
                      
                      test_bypass_fpsim <- test_bypass_fpsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_fpsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      for(i in 1:nrow(test_bypass_fpsim)){
                      #i <- 1
                       test_bypass_fpsim$sstar[i]  <- -(log10(1-(test_bypass_fpsim$Similarity[i] - 0.0001)) + log10(1-((test_bypass_fpsim$Concordance[i] - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_fpsim)[7] <- "S*"
                      test_bypass_fpsim$'S*'<- as.numeric(test_bypass_fpsim$'S*')
                      test_bypass_fpsim$'S*' <- round(test_bypass_fpsim$'S*', 3)

                      
                      #return(test_bypass_fpsim)}
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_fpsim$`LINCS Analog`<- makeLink5(test_bypass_fpsim$`LINCS Analog`)
                        #test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                        test_bypass_fpsim
                        return(test_bypass_fpsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      test_bypass_fpsim <<- test_bypass_fpsim
                      
                      
                      
                      # for (i in 1:nrow(test_bypass_fpsim)){
                      #   Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      # }
                      # test_bypass_fpsim$Concordance <- Concordance
                      # test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      # test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      # 
                      # test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      # output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,4], test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      # test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                      #test_bypass_fpsim <<- test_bypass_fpsim
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
                      #stopApp()
                    }
                    
                    #else minsim
                    
                    ################### MinSim and Define Target Gene #######################################
                    
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Define Target Gene"))
                    { #load("./minSim_apfp_RObjects.RData")
                      
                      print("in minsim and define target gene")
                      adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      #Validate(need(File.Exists(adds$datapath), message = "File Not Found"))
                      
                      print("Before sdfset_add")
                      sdfset_add <- read.SDFset(adds$datapath)
                      
                      print("Made it past sdfset_add")
                      
                      added_labels <- sdfid(sdfset_add)
                      
                      print("Made it past added_labels")
                      
                      #sdfset_add <- read.SDFset("./sdf_test.sdf")
                      #sdfset_add_smiles <- sdf2smiles(sdfset_add[4])
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      
                      print(dim(sdfset_add))
                      test_bypass_minsim <- bypass_clustering_minsim(sdfset_add, lsm_rows,max_cons_4)
                      test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim)
                      
                      # Edited 03_23_2021
                      test_bypass_minsim <- test_bypass_minsim[which(!(is.na(test_bypass_minsim[,2]))),]
                      
                      colnames(test_bypass_minsim) <- c("Compound", "LSM-ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_minsim))
                      test_bypass_minsim <- left_join(test_bypass_minsim, max_cons_4, by = "LSM-ID")
                      test_bypass_minsim <- test_bypass_minsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                     
                      
                      
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance
                      #test_bypass_minsim$Cell_Line <- Cell_Line
                      #test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      #test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance","Similarity")]
                      #colnames(test_bypass_minsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_minsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),]
                      
                      for (i in 1:nrow(test_bypass_minsim)){
                        test_bypass_minsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_minsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_minsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      for(i in 1:nrow(test_bypass_minsim)){
                      #i <- 1
                       test_bypass_minsim$sstar[i]  <- -(log10(1-(test_bypass_minsim$Similarity[i] - 0.0001)) + log10(1-((test_bypass_minsim$Concordance[i] - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_minsim)[7] <- "S*"
                      test_bypass_minsim$'S*'<- as.numeric(test_bypass_minsim$'S*')
                      test_bypass_minsim$'S*' <- round(test_bypass_minsim$'S*', 3)

                      
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_minsim$`LINCS Analog`<- makeLink5(test_bypass_minsim$`LINCS Analog`)
                        #test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,6], test_bypass_minsim[,5], decreasing=TRUE),]
                        test_bypass_minsim
                        return(test_bypass_minsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      test_bypass_minsim <<- test_bypass_minsim
                      #stopApp()
                    }
                    ############################################################################################
                    ######### SDF && fpSim && Upload a Signature ###############################################
                    
                    if ((input$Algorithm == "fpSim") & (input$Signature == "Upload a Signature")){
                      adds<<-input$AddedCompounds
                      
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      sdfset_add <- read.SDFset(adds$datapath)
                      added_labels <- sdfid(sdfset_add)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      
                      test_fpset <<- sdfset_add
                      print(dim(sdfset_add))
                      
                      #Edited 4-9-2021
                      #sums <- rowSums(sdfset_add)
                      #zeroes <- which(sums==0)
                      #sdfset_add <- sdfset_add[-zeroes,]
                      #added_labels <- added_labels[-zeroes]
                      added_labels <<- added_labels
                      lsm_rows <<- lsm_rows
                      
                      test_bypass_fpsim <-  bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      print("Made it this far.")
                      
                      test_bypass_fpsim <- cbind(added_labels, test_bypass_fpsim)
                      
                      #edited 4-9-2021
                      test_bypass_fpsim <- test_bypass_fpsim[which(!(is.na(test_bypass_fpsim[,2]))),]
                      
                      colnames(test_bypass_fpsim) <- c("Compound", "LSM-ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_fpsim))
                      
                      test_bypass_fpsim <- left_join(test_bypass_fpsim, max_cons_4, by = "LSM-ID")
                      test_bypass_fpsim <- test_bypass_fpsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      
                      #for (i in 1:nrow(test_bypass_fpsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_fpsim$Concordance <- Concordance
                      #test_bypass_fpsim$Cell_Line <- Cell_Line
                      #test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      #test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      #test_bypass_fpsim <- test_bypass_fpsim[,c("User-added Compound", "LINCS Analog", "Cell_Line", "Concordance", "Similarity")]
                      #test_bypass_fpsim <-test_bypass_fpsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance", "Similarity")]
                      #colnames(test_bypass_fpsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_fpsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                      
                      for (i in 1:nrow(test_bypass_fpsim)){
                        test_bypass_fpsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_fpsim[i,2]), 1]
                      }
                      
                      test_bypass_fpsim <- test_bypass_fpsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_fpsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      
                      for(i in 1:nrow(test_bypass_fpsim)){
                        #i <- 1
                        test_bypass_fpsim$sstar[i]  <- -(log10(1-(test_bypass_fpsim$Similarity[i] - 0.0001)) + log10(1-((test_bypass_fpsim$Concordance[i] - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_fpsim)[7] <- "S*"
                      test_bypass_fpsim$'S*'<- as.numeric(test_bypass_fpsim$'S*')
                      test_bypass_fpsim$'S*' <- round(test_bypass_fpsim$'S*', 3)
                      
                      
                      #return(test_bypass_fpsim)}
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_fpsim$`LINCS Analog`<- makeLink5(test_bypass_fpsim$`LINCS Analog`)
                        #test_bypass_fpsim[order(test_bypass_fpsim[,5], test_bypass_fpsim[,4], decreasing=TRUE),]
                        test_bypass_fpsim
                        return(test_bypass_fpsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      test_bypass_fpsim <<- test_bypass_fpsim
                      
                      
                      
                      # for (i in 1:nrow(test_bypass_fpsim)){
                      #   Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      # }
                      # test_bypass_fpsim$Concordance <- Concordance
                      # test_bypass_fpsim$Concordance <- round(test_bypass_fpsim$Concordance, 3)
                      # test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      # 
                      # test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      # output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,4], test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      # test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                      #test_bypass_fpsim <<- test_bypass_fpsim
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
                      #stopApp()
                    }
                    
                    #else minsim
                    
                    ################### SDF && MinSim && Upload a Signature #######################################
                    
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Upload a Signature"))
                    { #load("./minSim_apfp_RObjects.RData")
                      adds<<-input$AddedCompounds
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      #Validate(need(File.Exists(adds$datapath), message = "File Not Found"))
                      sdfset_add <- read.SDFset(adds$datapath)
                      added_labels <- sdfid(sdfset_add)
                      
                      #sdfset_add <- read.SDFset("./sdf_test.sdf")
                      #sdfset_add_smiles <- sdf2smiles(sdfset_add[4])
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      
                      print(dim(sdfset_add))
                      test_bypass_minsim <- bypass_clustering_minsim(sdfset_add, lsm_rows,max_cons_4)
                      test_bypass_minsim <- data.frame(added_labels, test_bypass_minsim)
                      
                      # Edited 03_23_2021
                      test_bypass_minsim <- test_bypass_minsim[which(!(is.na(test_bypass_minsim[,2]))),]
                      
                      test <<- test_bypass_minsim
                      
                      colnames(test_bypass_minsim) <- c("Compound", "LSM-ID", "Similarity")
                      Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      Cell_Line <- vector("character", length=nrow(test_bypass_minsim))
                      test_bypass_minsim <- left_join(test_bypass_minsim, max_cons_4, by = "LSM-ID")
                      test_bypass_minsim <- test_bypass_minsim[, c(1, 2, 5, 6, 3)]
                      test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      
                      
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #  Cell_Line[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Cell Line"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance
                      #test_bypass_minsim$Cell_Line <- Cell_Line
                      #test_bypass_minsim$Concordance <- round(test_bypass_minsim$Concordance, 3)
                      #test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Cell_Line", "Concordance","Similarity")]
                      #colnames(test_bypass_minsim) <-c("User-added Compound", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      colnames(test_bypass_minsim) <-c("User-added Candidate", "LINCS Analog", "Cell Line", "Concordance", "Similarity")
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,5], test_bypass_minsim[,4], decreasing=TRUE),]
                      
                      for (i in 1:nrow(test_bypass_minsim)){
                        test_bypass_minsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_minsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,6,3,4,5)]
                      colnames(test_bypass_minsim) <- c("User-added Candidate", "LINCS Analog", "LINCS Analog Name", "Cell Line", "Concordance", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      for(i in 1:nrow(test_bypass_minsim)){
                        #i <- 1
                        test_bypass_minsim$sstar[i]  <- -(log10(1-(test_bypass_minsim$Similarity[i] - 0.0001)) + log10(1-((abs(test_bypass_minsim$Concordance[i]) - 0.0001)/lincs_max)))
                      }
                      colnames(test_bypass_minsim)[7] <- "S*"
                      test_bypass_minsim$'S*'<- as.numeric(test_bypass_minsim$'S*')
                      test_bypass_minsim$'S*' <- round(test_bypass_minsim$'S*', 3)
                      
                      ########### Remove Duplicates 04/17/2023 ####################
                      #test_bypass_minsim <- test_bypass_minsim[!duplicated(test_bypass_minsim), ]
                      ##############################################################################
                        
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_minsim$`LINCS Analog`<- makeLink5(test_bypass_minsim$`LINCS Analog`)
                        #test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,6], test_bypass_minsim[,5], decreasing=TRUE),]
                        test_bypass_minsim
                        return(test_bypass_minsim)}
                        , caption="My Candidates Ranked", rownames=FALSE, escape=FALSE)
                      test_bypass_minsim <<- test_bypass_minsim
                      #stopApp()
                    }
                    
                    ################### MinSim and Find Analogs in LINCS #######################################
                    
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Find Analogs in LINCS"))
                    { #load("./minSim_apfp_RObjects.RData")
                      index_analogs <<- 1
                      adds<<-input$AddedCompounds
                      sdfset_add <- read.SDFset(adds$datapath)
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      added_labels <- sdfid(sdfset_add)
                      added_labels <<- added_labels
                      
                      
                      added_labels_2 <<- rep(added_labels, each=input$topN)
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      #topN <<- input$topN
                      test_bypass_minsim <- bypass_clustering_minsim_sim_search(sdfset_add, lsm_rows, input$topN)
                      test_bypass_minsim <- data.frame(added_labels_2, test_bypass_minsim)
                      colnames(test_bypass_minsim) <- c("User-added Compound", "LINCS Analog", "Similarity")
                      #test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      #test_bypass_minsim <- test_bypass_minsim[-which(test_bypass_minsim$Similarity < input$Similarity),]
                      test_bypass_minsim <- test_bypass_minsim %>% filter(Similarity > input$Similarity)
                      #Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance 
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      #test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),]
                      test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,1], -test_bypass_minsim[,3]),]
                      for (i in 1:nrow(test_bypass_minsim)){
                        test_bypass_minsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_minsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,4,3)]
                      colnames(test_bypass_minsim) <- c("User-added Compound", "LINCS Analog", "LINCS Analog Name", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      #return(test_bypass_minsim)}
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      #callback <- c(
                      #  "('$table.dataTable.display tbody tr:odd').css('background-color', 'yellow');")
                      
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_minsim$`LINCS Analog`<- makeLink5(test_bypass_minsim$`LINCS Analog`)
                        #test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),]
                        test_bypass_minsim
                        return(test_bypass_minsim)}
                        , caption="Analogs in LINCS", rownames=FALSE, escape=FALSE)
                      #stopApp()
                      lsm_rows <<- which(rownames(lincs_fps_2) %in% test_bypass_minsim[,2])
                      test_bypass_minsim <<- test_bypass_minsim
                    }
                    
                    
                    ################### MinSim and Find Analogs in LINCS #######################################
                    
                    else if ((input$Algorithm == "minSim") & (input$Signature == "Find Analogs in LINCS"))
                    { #load("./minSim_apfp_RObjects.RData")
                      index_analogs <<- 1
                      adds<<-input$AddedCompounds
                      sdfset_add <- read.SDFset(adds$datapath)
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      added_labels <- sdfid(sdfset_add)
                      added_labels <<- added_labels
                      
                      
                      added_labels_2 <<- rep(added_labels, each=input$topN)
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      #topN <<- input$topN
                      test_bypass_minsim <- bypass_clustering_minsim_sim_search(sdfset_add, lsm_rows, input$topN)
                      test_bypass_minsim <- data.frame(added_labels_2, test_bypass_minsim)
                      colnames(test_bypass_minsim) <- c("User-added Compound", "LINCS Analog", "Similarity")
                      #test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      #test_bypass_minsim <- test_bypass_minsim[-which(test_bypass_minsim$Similarity < input$Similarity),]
                      test_bypass_minsim <- test_bypass_minsim %>% filter(Similarity > input$Similarity)
                      #Concordance <- vector("numeric", length=nrow(test_bypass_minsim))
                      #for (i in 1:nrow(test_bypass_minsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_minsim$LSM_ID[i]), "Concordance"]
                      #}
                      #test_bypass_minsim$Concordance <- Concordance 
                      #test_bypass_minsim <- test_bypass_minsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      #test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),]
                      test_bypass_minsim <- test_bypass_minsim[order(test_bypass_minsim[,1], -test_bypass_minsim[,3]),]
                      for (i in 1:nrow(test_bypass_minsim)){
                        test_bypass_minsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_minsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,4,3)]
                      colnames(test_bypass_minsim) <- c("User-added Compound", "LINCS Analog", "LINCS Analog Name", "Similarity")
                      test_bypass_minsim$Similarity <- round(test_bypass_minsim$Similarity, 3)
                      #return(test_bypass_minsim)}
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      #callback <- c(
                      #  "('$table.dataTable.display tbody tr:odd').css('background-color', 'yellow');")
                      
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_minsim$`LINCS Analog`<- makeLink5(test_bypass_minsim$`LINCS Analog`)
                        #test_bypass_minsim[order(test_bypass_minsim[,3], decreasing=TRUE),]
                        test_bypass_minsim
                        return(test_bypass_minsim)}
                        , caption="Analogs in LINCS", rownames=FALSE, escape=FALSE)
                      #stopApp()
                      lsm_rows <<- which(rownames(lincs_fps_2) %in% test_bypass_minsim[,2])
                      test_bypass_minsim <<- test_bypass_minsim
                    }
                    
                    ################### FPSim and Find Analogs in LINCS #######################################
                    
                    else if ((input$Algorithm == "fpSim") & (input$Signature == "Find Analogs in LINCS"))
                    { #load("./minSim_apfp_RObjects.RData")
                      index_analogs <<- 1
                      adds<<-input$AddedCompounds
                      sdfset_add <- read.SDFset(adds$datapath)
                      #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
                      added_labels <- sdfid(sdfset_add)
                      added_labels_2 <<- rep(added_labels, each=input$topN)
                      
                      ####################################
                      # ChemmineR Atom Pair Fingerprints
                      ##################################
                      apset_added <-sdf2ap(sdfset_add)
                      fpset_added <- as.matrix(desc2fp(apset_added))
                      
                      sdfset_add <- as.matrix(fpset_added)
                      print(dim(sdfset_add))
                      #test_bypass_fpsim <- bypass_clustering_fpsim(sdfset_add, lsm_rows)
                      #test_bypass_fpsim <- data.frame(added_labels, test_bypass_fpsim)
                      #colnames(test_bypass_fpsim) <- c("User-added Compound", "LINCS Analog", "Similarity")
                      
                      #Edited 4-9-2021
                      sums <- rowSums(sdfset_add)
                      zeroes <- which(sums==0)
                      sdfset_add <- sdfset_add[-zeroes,]
                      added_labels_2 <- added_labels_2[-zeroes]
                      
                      test_bypass_fpsim <- bypass_clustering_fpsim_sim_search(sdfset_add, lsm_rows, input$topN)
                      test_bypass_fpsim <- data.frame(added_labels_2, test_bypass_fpsim)
                      colnames(test_bypass_fpsim) <- c("User-added Compound", "LINCS Analog", "Similarity")
                      #test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      #test_bypass_fpsim <- test_bypass_fpsim[-which(test_bypass_fpsim$Similarity < input$Similarity),]
                      test_bypass_fpsim <- test_bypass_fpsim %>% filter(Similarity > input$Similarity)
                      #Concordance <- vector("numeric", length=nrow(test_bypass_fpsim))
                      #for (i in 1:nrow(test_bypass_fpsim)){
                      #  Concordance[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == test_bypass_fpsim$LSM_ID[i]), "Concordance"]
                      #}
                      #test_bypass_fpsim$Concordance <- Concordance 
                      #test_bypass_fpsim <- test_bypass_fpsim[,c("Compound", "LSM_ID", "Concordance","Similarity")]
                      #df_centroid$Representative[[i]] <- makeLink(df_centroid$Representative[[i]])
                      #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim[order(test_bypass_fpsim[,3], decreasing=TRUE),], caption="My Candidates Ranked", rownames=FALSE)
                      #test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,3], decreasing=TRUE),]
                      test_bypass_fpsim <- test_bypass_fpsim[order(test_bypass_fpsim[,1],-test_bypass_fpsim[,3]),]
                      for (i in 1:nrow(test_bypass_fpsim)){
                        test_bypass_fpsim$Analog_Name[i] <- lincs_compounds[which(lincs_compounds[,2] == test_bypass_fpsim[i,2]), 1]
                      }
                      
                      test_bypass_minsim <- test_bypass_minsim[,c(1,2,4,3)]
                      colnames(test_bypass_minsim) <- c("User-added Compound", "LINCS Analog", "LINCS Analog Name", "Similarity")
                      test_bypass_fpsim$Similarity <- round(test_bypass_fpsim$Similarity, 3)
                      shinyjs::hideElement("sim_search")
                      shinyjs::hideElement("logo")
                      shinyjs::showElement("CANDIDATES")
                      output$CANDIDATES <- DT::renderDataTable({
                        test_bypass_fpsim$`LINCS Analog`<- makeLink5(test_bypass_fpsim$`LINCS Analog`)
                        #test_bypass_fpsim[order(test_bypass_fpsim[,3], decreasing=TRUE),]
                        test_bypass_fpsim
                        return(test_bypass_fpsim)},
                        caption="Analogs in LINCS", rownames=FALSE, escape=FALSE)
                      
                      test_bypass_fpsim <<- test_bypass_fpsim
                      #stopApp()
                      lsm_rows <<- which(rownames(lincs_fps_2) %in% test_bypass_fpsim[,2])
                      test_bypass_fpsim <<- test_bypass_fpsim
                    }
                    # cluster_compounds()
                  # print("LINCS clustering complete!")
                  # output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
                }
    }#)
    }
#    else{print("Clustering bypassed")
    else{print("No added compounds, must cluster.")
     
    }
      
      
      
      
      
          
        
    
#        output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5))
       incProgress(7/7, message = "Complete!")
    })
      if (input$Signature != "Find Analogs in LINCS"){
      output$downloadCandidates <- renderUI({
         #if((input$AddedCompounds == 1))
        if(!(is.null(input$AddedCompounds))){
         downloadButton("CandidatesDownload", label = "My Candidates Ranked")}
       })
      }
      else if (input$Signature == "Find Analogs in LINCS"){
      output$downloadAnalogs <- renderUI({
        #if((input$AddedCompounds == 1))
        if(!(is.null(input$AddedCompounds))){
          downloadButton("AnalogsDownload", label = "Analogs in LINCS")}
      })
      }
      #output$downloadLINCS<- renderUI({
      #  if((input$Signature == "Define Target Gene") | (input$Signature == "Upload a Signature")){
      #  downloadButton("LINCSDownload", label = "LINCS Candidates Ranked")}
      #})
      if (input$Signature == "Define Target Gene"){
      output$downloadLINCS<- renderUI({
        if((input$Signature == "Define Target Gene") ){
          downloadButton("LINCSDownload", label = "LINCS Candidates Ranked")}
      })
      }
      # else if (input$Signature == "Upload a Signature"){
      # output$downloadUpload <- renderUI({
      #   if((input$Signature == "Upload a Signature")){
      #     downloadButton("downloadSigs", label = "LINCS Signatures Ranked")}
      # })
      # }
      # output$downloadSigs <- downloadHandler(
      #   filename = function() {
      #     'lincs_signatures_ranked.csv'
      #   },
      #   content = function(filename){
      #     write.csv(lincs_output, filename, row.names = FALSE)
      #   }
      # )
      
    }) # End input$Go

 
 
###################### Pasted in 02/06/2021 ###################################################################################################
###################### Chemical Similarity Analysis (Tab 2) ###################################################################################
###############################################################################################################################################
 
  observeEvent(input$Cluster,{
    withProgress(message = "Clustering.  May take a while...", value = 0, {
      #incProgress(1/1, message = "Clustering")
      shinyjs::hideElement("Representatives")
      shinyjs::hideElement("check_analogs")
      shinyjs::hideElement("check_clustersize")
      shinyjs::hideElement("check_clusternumber")
      if (index_analogs < 1){
      #if (index_analogs == 1){
        shinyjs::showElement("check_analogs")
        output$check_analogs <- renderUI({
          p("Clustering cannot be performed.  Please run Signature Connectivity Analysis (First Tab) first.")
        })
        return()
      }
      if (as.integer(input$ClusterSize) < 3){
        shinyjs::showElement("check_clustersize")
        output$check_clustersize<- renderUI({
          p("Minimum cluster size is 3.  Please enter another minimum cluster size")
        })
        return()
      }
      #if (!is.null(input$AddedCompounds)){  
      #  if(!(grepl(".sdf", adds$name))){
      pdf(NULL)
      if (!is.null(input$AddedCompounds)){
      #if (!(input$file_input == "")){
      #if (!(is.null(input$file_input))){
         
         #if(!(grepl(".sdf", adds$name))){
        #if(!(grepl(".sdf", input$file_input$name))){
          # 
         if(!(grepl(".sdf", adds$name))){
             adds_SMI <<- read.SMIset(adds$datapath)
          
          #incProgress(1/1, message = "In added, no sdf")
          incProgress(1/1, message = "Running SAR Analysis.  This may take a while...")
          #adds_SMI <<- read.SMIset(adds$datapath)
          #adds_SMI <<- read.SMIset(paste("../", adds$name, sep=""))
          #adds_SMI <- read.SMIset(paste("/srv/shiny-server/userfile/", input$file_input$name, sep=""))
          #adds_csv <- read.csv(adds$datapath, header = FALSE, sep = "\t")
         
          #adds_SMI <<- read.SMIset(paste(getwd(), adds$name, sep="/"))
          
            #adds_csv <- read.csv(paste(getwd(), adds$name, sep="/"), header = FALSE, sep = "\t")
          #adds_csv <- adds_csv[,c(2,1)]
          #colnames(adds_csv) <- c("LSM_ID", "SMILES")
          #adds_csv <<- adds_csv
          #test2 <<- rbind(adds_csv, display)
          #colnames(test2) <<- c("Compound_ID", "SMILES")
          
          sdfset_add <- smiles2sdf(adds_SMI)
          
          added_labels <- sdfid(sdfset_add)
          sdfset_add <<- sdfset_add
          ####################################
          # ChemmineR Atom Pair Fingerprints
          ##################################
          apset_added <-sdf2ap(sdfset_add)
          fpset_added <- as.matrix(desc2fp(apset_added))
          rownames(fpset_added) <- added_labels
          
          #sdfset_add <- as.matrix(fpset_added)
          #print(dim(sdfset_add))
          #lsm_rows <- 1:100
          #if(((nrow(fpset_added) + length(lsm_rows)) > 5000)){
          #  print("Too many compunds to cluster.  Reduce Number of Added Compounds and/or increase concordance threshold")
          #  shinyjs::showElement("cluster_check")
          #  output$cluster_check <- renderUI({
          #    p("Too many compunds to cluster.  Reduce Number of Added Compounds and/or increase concordance threshold")
          #  })
          #  return()
          #}
          
          fpsim_cluster(fpset_added, lsm_rows, input$NumberCmpds)
          #incProgress(1/2, detail = "Clustering")
          #test_cluster_fpsim <- cbind(added_labels, test_bypass_fpsim)
          #colnames(test_bypass_fpsim) <- c("Compound", "LSM_ID", "Similarity")
          #test_bypass_fpsim <<- test_bypass_fpsim
          #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
          
          if (nrow(fpset_added) >= input$NumberCmpds){
            color_added <- replicate(nrow(fpset_added),"blue")
            colorbar <- color_added[1:input$NumberCmpds]
            #color_lincs <- replicate((nrow(simMA) - nrow(fpset_added)), "green")
          }
          else {
            color_added <- replicate(nrow(fpset_added),"blue")
            color_lincs <- replicate((nrow(simMA) - nrow(fpset_added)), "green")
            colorbar <- c(color_added, color_lincs)
          }    
          #color_added <- replicate(nrow(fpset_added),"blue")
          #color_lincs <- replicate((nrow(simMA) - nrow(fpset_added)), "green")
          
          #colorbar <- c(color_added, color_lincs)
          
          #colorbar <<- colorbar
          shinyjs::showElement("distPlot")
          #pdf(NULL)
          output$distPlot <<- renderPlot({
            par(mar=c(5.1,0.5,4.1,0), xpd=FALSE)
            heatmap.2(1-simMA, Rowv=as.dendrogram(hc),
                      Colv=as.dendrogram(hc),
                      dendrogram="column",
                      col=colorpanel(40, "white","yellow","red"),
                      density.info="none",
                      trace="none",
                      labCol=FALSE,
                      labRow=FALSE,
                      main="Clustering of User-added and LINCS Candidates",
                      ColSideColors = colorbar
                      )
            #cexRow=0.5,
            #cexCol=0.5)
            legend("left",
                   legend=c("My Candidates", "LINCS"),
                   col=c("blue", "green"),
                   pch=19,
                   pt.cex=2)
          })
          
          
          #sdf_smiles <<- c(sdf_smiles, sdfset_add)
          #sdf_smiles <<- c(sdf_smiles, adds_SDF)
          #fpsim_cluster(fps)
          print("Your compounds have been clustered with LINCS compounds")
          #source("lib/ColorMap.R", local = TRUE)
          #
          #color_dend()
        }
        
        ###############Critical section
        else if(grepl(".sdf", adds$name)){
        #else if(grepl(".sdf", input$file_input$name)){
          incProgress(1/1, message = "Running SAR Analysis.  This may take a while...")
          #incProgress(1/3, message = "In added, sdf")
          print("Added compounds in sdf format")
          #sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
          #sdfset_add <<- read.SDFset(paste(getwd(), adds$name, sep="/"))
          #lsm_smiles_2 <- lsm_smiles
          #sdfset_add <<- read.SDFset(adds$datapath)
          #added_smiles <<- sdf2smiles(sdfset_add)
          #added_labels <- cid(added_smiles)
          #added_smiles_2 <- as.character(added_smiles[1:length(added_smiles)])
          #added_smiles_3 <- cbind(added_labels, added_smiles_2)
          
          #added_smiles_4 <- unname(added_smiles_3)
          #colnames(added_smiles_3) <- c("Compound_ID", "SMILES")
          
          #lsm_smiles_2_labels <- cid(lsm_smiles_2)
          #lsm_smiles_3 <- as.character(lsm_smiles_2[1:length(lsm_smiles_2)])
          #lsm_smiles_3b <- unname(lsm_smiles_3)
          #lsm_smiles_4 <- cbind(lsm_smiles_2_labels, lsm_smiles_3b)
          #colnames(lsm_smiles_4) <- c("Compound_ID", "SMILES")
          #colnames(display) <-  c("Compound_ID", "SMILES")
          
          #test2 <<- rbind(added_smiles_3, display)
          
          #colnames(test2) <<- c("Compound_ID", "SMILES")
          
          
          #sdf_smiles <<- c(sdf_smiles, sdfset_add)
          #cluster_compounds()
          #print("Your compounds have been clustered with LINCS compounds")
          
          ############# Hierarchical CLustering ###########################
          #incProgress(2/3, message = "In added, sdf, after ColorMap.R")
          
          #source("./lib/ColorMap.R", local = TRUE)
          
          #incProgress(2/3, message = "In added, sdf, after ColorMap.R")
          #adds_SMI <<- added_smiles
          #color_dend()
          
          #adds<<-input$AddedCompounds
          #sdfset_add <- read.SDFset(paste(getwd(), adds$name, sep="/"))
          sdfset_add <- read.SDFset(adds$datapath)
          #sdfset_add <- read.SDFset(paste("../", adds$name, sep=""))
          #sdfset_add <- read.SDFset(paste("/srv/shiny-server/userfile/", input$file_input$name, sep=""))
          #print(sdfid(sdfset_add))
          
          if (sdfid(sdfset_add)[1] == "")
            added_labels <- cid(sdfset_add)
          else{
            added_labels <- sdfid(sdfset_add)
          }
          added_labels <<- added_labels
          sdfset_add <<- sdfset_add
          ####################################
          # ChemmineR Atom Pair Fingerprints
          ##################################
          apset_added <-sdf2ap(sdfset_add)
          fpset_added <- as.matrix(desc2fp(apset_added))
          print(rownames(fpset_added))
          
          rownames(fpset_added)<- added_labels
          
          fpset_added <<- fpset_added
          #lsm_rows <- 1:100
          #if(((nrow(fpset_added) + length(lsm_rows)) > 5000)){
          #  print("Too many compunds to cluster.  Reduce Number of Added Compounds and/or increase concordance threshold")
          #  shinyjs::showElement("cluster_check")
          #  output$cluster_check <- renderUI({
          #    p("Too many compunds to cluster.  Reduce Number of Added Compounds and/or increase concordance threshold")
          #  })
          #  return()
          #}
          #sdfset_add <- as.matrix(fpset_added)
          #print(dim(sdfset_add))
          #lsm_rows <- 1:100
          #incProgress(2/3, message = "In added, sdf, before fpsim_cluster")
          print(paste("Number of Compounds", input$NumberCmpds))
          fpsim_cluster(fpset_added, lsm_rows, input$NumberCmpds)
          
          #fpsim_cluster(fpset_added, lsm_rows, TRUE)
          #incProgress(1/3, message = "In added, sdf, after fpsim_cluster")
          #incProgress(1/2, detail = "Clustering")
          #test_cluster_fpsim <- cbind(added_labels, test_bypass_fpsim)
          #colnames(test_bypass_fpsim) <- c("Compound", "LSM_ID", "Similarity")
          #test_bypass_fpsim <<- test_bypass_fpsim
          #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
          
          
          if (nrow(fpset_added) >= input$NumberCmpds){
          print("At 3033.")
          #color_added <- replicate(nrow(fpset_added),"blue")
          color_added <- replicate(input$NumberCmpds,"blue")
          colorbar <- color_added
          colorbar_test <<- colorbar
          print("At 3036.")
          #print("Length colorbar: ", length(col))
          #color_lincs <- replicate((nrow(simMA) - nrow(fpset_added)), "green")
          }
          else {
            color_added <- replicate(nrow(fpset_added),"blue")
            color_lincs <- replicate((nrow(simMA) - nrow(fpset_added)), "green")
            colorbar <- c(color_added, color_lincs)
          }          
          #colorbar <<- colorbar
          incProgress(3/3, message = "In added, sdf before distPlot")
          shinyjs::showElement("distPlot")
          output$distPlot <<- renderPlot({#pdf(NULL)
            par(mar=c(5.1,0.5,4.1,0), xpd=FALSE)
            heatmap.2(1-simMA, Rowv=as.dendrogram(hc),
                      Colv=as.dendrogram(hc),
                      col=colorpanel(40, "white","yellow","red"),
                      dendrogram="column",
                      density.info="none",
                      trace="none",
                      labCol=FALSE,
                      labRow=FALSE,
                      ColSideColors = colorbar,
                      main="Clustering of User-added and LINCS Candidates")
            #cexRow=0.5,
            #cexCol=0.5)
            legend("left",
                   legend=c("My Candidates", "LINCS"),
                   col=c("blue", "green"),
                   pch=19,
                   pt.cex=2)
          })
          
          #incProgress(1/1, message = "In added, no sdf after heatmap")
          incProgress(1/1, message = "Hierarchical Clustering Complete")
        }
        
        
        ########## MDS Plot #########################################################################
        #added_smiles <<- sdf2smiles(sdfset_add)
        #added_labels <- cid(added_smiles)
        
        ##################
        # Added 1/18/2021
        
        ################
        
        #added_smiles_2 <- as.character(added_smiles[1:length(added_smiles)])
        #added_smiles_3 <- cbind(added_labels, added_smiles_2)
        
        #added_smiles_4 <- unname(added_smiles_3)
        #colnames(added_smiles_3) <- c("Compound_ID", "SMILES")
        
        ##########Edit 1-17-2021
        #added_smiles_3 <<- added_smiles_3
        
        #################################################
        #lsm_smiles_2_labels <- cid(lsm_smiles_2)
        #lsm_smiles_3 <- as.character(lsm_smiles_2[1:length(lsm_smiles_2)])
        #lsm_smiles_3b <- unname(lsm_smiles_3)
        #lsm_smiles_4 <- cbind(lsm_smiles_2_labels, lsm_smiles_3b)
        #colnames(lsm_smiles_4) <- c("Compound_ID", "SMILES")
        
        #if (!(input$file_input == "")){
        #if (!(is.null(input$file_input))){
        if (!(is.null(adds))){    
          #Changed 4/10/22
          ###############################
          #added_ids <- sdfid(sdfset_add)
          
          #added_ids <- cid(sdfset_add)
          if (sdfid(sdfset_add)[1] == "")
            added_ids <- cid(sdfset_add)
          else{
            print("At 3109")
            added_ids <- sdfid(sdfset_add)
            print("at 3111")
          }
          
          #######################################################
          #colnames(display) <-  c("Compound_ID", "SMILES")
          
          
          # Changed 11/15/2023
          display <- display[,1]
          
          test2 <- c(added_ids, display)
          added_ids <<- added_ids
          
          #colnames(test2) <<- c("Compound_ID", "SMILES")
          test2 <- as.data.frame(test2)
          colnames(test2)[1] <- "Compound_ID"
          test2 <<- test2
          Added_Label <<- "My Candidates"
        }
        else{
          display <- display
          #test2 <- as.data.frame(d)
          #colnames(test2)[1] <- "Compound_ID"
          #test2 <<- test2
        } 
        withProgress(message = "Plotting Representatives", value=0, {
          #source("lib/Centroid.R", local=TRUE)
          source("./lib/Centroid2.R", local=TRUE)
          #source("/srv/shiny-server/lib/Centroid2.R")
          #Add a threshold option for cut_tree and Cluster size
          #6. Cut tree at user specified similarity score, identify representative based on minimum distance from all other members
          
          print("Right before centroid2. 3138")
          cut_tree((1-as.numeric(input$CutHeight)))
          
          #6.1. Include Added compounds
          print("Right after centroid2. 3143")
          if (!is.null(input$AddedCompounds)){
          #if (!(input$file_input == "")){  
            find_centroid_adds(ClusterMembers, (as.numeric(input$ClusterSize)-1))
          }
          #6.2. No compounds were added
          else{
            find_centroid(ClusterMembers, (as.numeric(input$ClusterSize)-1))
          }
          print("Right after find_centroid.  3152")
          df_centroid <- as.data.frame(unlist(centroid), stringsAsFactors = FALSE)
          colnames(df_centroid)="Representative"
          df_centroid<-df_centroid
          centroid_clusters(centroid, ClusterMembers)
          df_centroid <- cbind(centroid = df_centroid, Cluster = unlist(ClusterCentroid))
          print("Right after df_centroid. 3158" )
          
          makeLink <- function(val) {
            #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
            paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
          }
          
          makeLink2 <- function(val) {
            #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
            paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
          }
          
          
          
          
          
          if (!is.null(input$AddedCompounds)){
          #if (!(input$file_input == "")){
          #if (!(is.null(input$file_input))){
            #test5 <- data.frame("Compound_ID"=0, "SMILES"=0)
            test5 <- data.frame("Compound_ID"=0)
            print("In 3179.")
            df_centroid <<- df_centroid
            for (i in 1:length(df_centroid$Representative)){    
              if (length(which(test2$Compound_ID %in% df_centroid$Representative[[i]]))>0){
                print("before test 3.  3182")
                test3 <- which(test2$Compound_ID %in% df_centroid$Representative[[i]])
                print("after test3. 3184")
                for (j in 1:length(test3)){
                  #test4 <- test2[test3[[j]], c("Compound_ID", "SMILES")]
                  test4 <- test2[test3[[j]], "Compound_ID"]
                  test4 <- as.character(test4)
                  #if (test5$Compound_ID[[1]]!=0){
                  if (test5[1]!=0){  
                    test5 <- c(test5, test4)
                  }
                  else{
                    test5 <- test4
                  }
                }
              }
            }
          }
          else {
            print("Before test 5.  3202")
            test5 <- data.frame("LSM_ID"=0, "SMILES"=0)
            for (i in 1:length(df_centroid$Representative)){    
              if (length(which(display$LSM_ID %in% df_centroid$Representative[[i]]))>0){
                test3 <- which(display$LSM_ID %in% df_centroid$Representative[[i]])
                for (j in 1:length(test3)){
                  test4 <- display[test3[[j]], c("LSM_ID", "SMILES")]
                  if (test5$LSM_ID[[1]]!=0){
                    test5 <- rbind(test5, test4)
                  }
                  else{
                    test5 <- test4
                  }
                }
              }
            }
          }
          #colnames(test5) <- c("Representative", "SMILES")
          test5 <<- test5
          
          
          #df_centroid <- merge(test5, df_centroid, stringsAsFactors = FALSE)
          df_centroid <- df_centroid[which(df_centroid$Representative %in% test5),]
          df_centroid <- data.frame(lapply(df_centroid, as.character), stringsAsFactors = FALSE)
          df_centroid$Cluster <- as.integer(df_centroid$Cluster)
          df_centroid <- df_centroid[order(df_centroid$Cluster),]
          #df_centroid<<-df_centroid
          #df_centroid <- df_centroid[,-2]
          # for(i in 1:nrow(df_centroid)){
          #   
          #   df_centroid$Compound[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == df_centroid$Representative[i]), "Compound"]
          # } 
          
          for(i in 1:nrow(df_centroid)){
            if (grepl("LSM", df_centroid$Representative[i])){
              df_centroid$Compound[i] <- max_cons_4[which(max_cons_4$`LSM-ID` == df_centroid$Representative[i]), "Compound"]
            }
            else{
              df_centroid$Compound[i] <- df_centroid$Representative[i]
            }
          }
          
          
          df_centroid <- df_centroid[,c(1,3,2)]
          df_centroid <<- df_centroid
          ####### Changed ##
          # df_centroid <- df_centroid[-which(duplicated(df_centroid[,1])),]
          
          ##########################
          
          df_centroid <- df_centroid
          
          ########[-which(duplicated(df$V1)),]        
          shinyjs::showElement("Representatives")
          
          output$Representatives <<- DT::renderDataTable({
            for (i in 1:length(df_centroid$Representative)){
              if (!is.na((grep("LSM", df_centroid$Representative[[i]])) && (grep("LSM", df_centroid$Representative[[i]])==1))){
                df_centroid$Representative[[i]] <- makeLink(df_centroid$Representative[[i]])
              }
              else if(!is.null(grep("LSM", df_centroid$Representative[[i]]))){
                df_centroid$Representative[[i]] <- gsub("ZINC0", "ZINC", df_centroid$Representative[[i]])
                df_centroid$Representative[[i]] <- makeLink2(df_centroid$Representative[[i]])
              }
            }
            #colnames(df_centroid) <-c("Representative LSM ID", "Representative Name", "Cluster")
            colnames(df_centroid) <-c("Representative ID", "Representative Name", "Cluster")
            
            ######### Added
            
            return(df_centroid)
            #####################
            #df_centroid <- df_centroid[,-2]
            df_centroid <<- df_centroid
          }, escape = FALSE, rownames = FALSE)
          #print(df_centroid$Representative)
          df_centroid <<- df_centroid
          output$downloadReps <- renderUI({
            downloadButton("RepDownload", label = "Representatives")
          })
          
          
          print(paste("Got the centroids as ", centroid[[1]], sep=""))
          print("At 3275")
          incProgress(1/2, detail = "Plotting MDS")
          
          source("lib/MDS.R")
          #source("/srv/shiny-server/lib/MDS.R")
          #7. Generate MDS plot from representatives distance to one another, pie chart radius corresponds to cluster size
          if (length(centroid) < 3){
            shinyjs::showElement("check_clusternumber")
            output$check_clusternumber <- renderUI({
              shinyjs::hideElement("distPlot")
              #shinyjs::hideElement("MDSPlot")
              shinyjs::hideElement("Representatives")
              p("No clusters of that size.  Please enter a smaller minimum cluster size.")
              
              
            })
            return()
          }
          MDS_plot(centroid)
          #            centroid_clusters(centroid, ClusterMembers)
          #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(ClusterCentroid), hjust=1, vjust=2)) 
          #                             + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
          
          ############## Edited 1/18/2021
          shinyjs::showElement("MDSPlot")
          output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2))
                                       + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "My Candidates")) +
                                         scale_fill_manual(breaks=c("LINCS", "My Candidates"), values=c("green", "blue")) +
                                         ggtitle("Multidimentional Scaling Plot of User-Added and LINCS Candidates"))
          
          #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2))
          #                             + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("My Candidates", "LINCS")) +
          #                               scale_fill_manual(values=c("green", "blue")) +
          #                               ggtitle("Multidimentional Scaling Plot of User-Added and LINCS Candidates"))
          #output
          
          ######################################################
          #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid), hjust=1, vjust=2))
          #                          + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
          incProgress(2/2, detail = "Complete!")                                                         
        })
      } # if statement for added comounds for Clustering
      ##########################################################################################
      # Clustering - No added compounds
      else{
        incProgress(1/1, message = "In no added, no sdf")
        #fpset_added <- NULL
        #lsm_rows <- 1:100
        #if(length(lsm_rows) > 5000){
        #  print("Too many compunds to cluster. Increase concordance threshold")
        #  shinyjs::showElement("cluster_check")
        #  output$cluster_check <- renderUI({
        #    p("Too many compunds to cluster.  Increase concordance threshold")
        #  })
        #  return()
        #}
        fpsim_cluster_no_added(lsm_rows,input$NumberCmpds)
        #print("Got this far.")
        #incProgress(1/2, detail = "Clustering")
        #test_cluster_fpsim <- cbind(added_labels, test_bypass_fpsim)
        #colnames(test_bypass_fpsim) <- c("Compound", "LSM_ID", "Similarity")
        #test_bypass_fpsim <<- test_bypass_fpsim
        #output$CANDIDATES <- DT::renderDataTable(test_bypass_fpsim)
        
        
        #color_added <- replicate(nrow(fpset_added),"blue")
        color_lincs <- replicate(nrow(simMA), "green")
        #colorbar <- c(color_lincs)
        #colorbar <<- colorbar
        shinyjs::showElement("distPlot")
        output$distPlot <<- renderPlot(heatmap.2(1-simMA, Rowv=as.dendrogram(hc),
                                                 Colv=as.dendrogram(hc),
                                                 col=colorpanel(40, "white","yellow","red"),
                                                 density.info="none",
                                                 trace="none",
                                                 labCol=FALSE,
                                                 labRow=FALSE,
                                                 main="Clustering of LINCS Candidates"
                                                 #ColSideColors = colorbar
                                                 #cexRow=0.5,
                                                 #cexCol=0.5)
        )) # End heatmap
        
        colnames(display) <-  c("LSM_ID", "SMILES")
        
        test2 <<- display
        
        #colnames(test2) <<- c("Compound_ID", "SMILES")
        
        #Added_Label <<- input$AddedLabel
        withProgress(message = "Plotting Representatives", value=0, {
          #source("lib/Centroid.R", local=TRUE)
          source("lib/Centroid2.R", local=TRUE)
          #source("/srv/shiny-server/lib/Centroid2.R")
          #Add a threshold option for cut_tree and Cluster size
          #6. Cut tree at user specified similarity score, identify representative based on minimum distance from all other members
          cut_tree((1-as.numeric(input$CutHeight)))
          
          
          ########### Edited 1/18/2021
          #6.1. Include Added compounds
          if (!is.null(input$AddedCompounds)){
          #if (!(is.null(input$file_input))){  
            find_centroid_adds(ClusterMembers, (as.numeric(input$ClusterSize)-1))
          }
          #6.2. No compounds were added
          else{
            find_centroid(ClusterMembers, (as.numeric(input$ClusterSize)-1))
          }
          
          
          #####################################################################
          df_centroid <- as.data.frame(unlist(centroid), stringsAsFactors = FALSE)
          colnames(df_centroid)="Representative"
          df_centroid<-df_centroid
          centroid_clusters(centroid, ClusterMembers)
          df_centroid <- cbind(centroid = df_centroid, Cluster = unlist(ClusterCentroid))
          
          
          makeLink <- function(val) {
            #paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val, "</a>", sep="")
            paste("<a href='", "http://lincsportal.ccs.miami.edu/SmallMolecules/view/", val,"'>", val, "</a>", sep="")
          }
          
          makeLink2 <- function(val) {
            #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
            paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
          }
          
          
          
          test5 <- data.frame("LSM_ID"=0, "SMILES"=0)
          for (i in 1:length(df_centroid$Representative)){    
            if (length(which(display$LSM_ID %in% df_centroid$Representative[[i]]))>0){
              test3 <- which(display$LSM_ID %in% df_centroid$Representative[[i]])
              for (j in 1:length(test3)){
                test4 <- display[test3[[j]], c("LSM_ID", "SMILES")]
                if (test5$LSM_ID[[1]]!=0){
                  test5 <- rbind(test5, test4)
                }
                else{
                  test5 <- test4
                }
              } # end inner for loop
            }
          } # end outer for loop
          
          colnames(test5) <- c("Representative", "SMILES")
          test_test5_before_remove_dups <<- test5
          
          if(length(which(duplicated(test5$Representative))) > 0){
            test5 <- test5[-which(duplicated(test5$Representative)),]
          }
          test_test5_after_dups <<- test5
          test5$Representative <- as.character(test5$Representative)
          
          #df_centroid <<- df_centroid
          #print("Got this far.")
          test_df_centroid <<- df_centroid
          test_test5 <<- test5
          df_centroid <- merge(test5, df_centroid, by.x="Representative", by.y="Representative", stringsAsFactors = FALSE)
          #df_centroid
          #print("Got this far")
          
          
          #df_centroid <<- df_centroid
          
          
          
          #print("Got this far.") 
          
          #df_centroid <- data.frame(lapply(df_centroid, as.character), stringsAsFactors = FALSE)
          
          df_centroid$Cluster <- as.integer(df_centroid$Cluster)
          
          
          df_centroid <- df_centroid[order(df_centroid$Cluster),]
          
          # print("Got this far.")
          # #df_centroid<<-df_centroid
          # #
          
          if(index == 2){ 
            df_centroid <- df_centroid[,-2]
            
            #df_test <<- df_centroid
            max_cons_5 <- max_cons_4
            
            if(length(which(duplicated(max_cons_5$'Compound'))) > 0){
              max_cons_5 <- max_cons_5[-which(duplicated(max_cons_5[,1])),]
              #max_con_5_if <<- max_cons_5
            }
            for(i in 1:nrow(df_centroid)){
              df_centroid$Compound[i] <- max_cons_5[which(max_cons_5$`LSM-ID` == df_centroid$Representative[i]), "Compound"]
            }
            
            #df_test2 <- df_test
            #df_test2 <- df_test[,-2] 
            #Compound <- vector("character", length=nrow(df_test2))
            #for(i in 1:nrow(df_test2)){
            #  df_test2$Compound[i] <- max_cons_4[which(max_cons_4$'LSM-ID' == df_test2$Representative[i]), "Compound"]
            #}
            #
            #
            
            
            df_centroid <- df_centroid[,c(1,3,2)]
          }
          else if (index == 1) {  
            
            #df_centroid <- df_centroid[,-2]
            test_beginning_df_centroid <<- df_centroid
            #df_test <<- df_centroid
            max_cons_5 <- max_cons_4
            
            if(length(which(duplicated(max_cons_5$'LSM-ID'))) > 0){
              max_cons_5 <- max_cons_5[-which(duplicated(max_cons_5$'LSM-ID')),]
              #max_con_5_if <<- max_cons_5
            }
            test_max_cons_5 <<- max_cons_5
            for(i in 1:nrow(df_centroid)){
              df_centroid$Compound[i] <- max_cons_5[which(max_cons_5$`LSM-ID` == df_centroid$Representative[i]), "Compound"]
            }  
            
            df_centroid <- df_centroid[,c(1,4,3)]
            #test_df_centroid <<- df_centroid
          }
          
          #################### Testing
          #
          #for(i in 1:nrow(test_df_centroid)){
          #   test_df_centroid$Compound[i] <- test_max_cons_5[which(test_max_cons_5$`LSM-ID` == test_df_centroid$Representative[i]), "Compound"]
          # } 
          
          #df_centroid <- df_centroid[-which(duplicated(df_centroid[,2])),]
          #df_centroid <- df_centroid
          #print("Got this far.")
          
          #
          #test_df_centroid <<- df_centroid
          test_max_cons_5 <- max_cons_5
          ########[-which(duplicated(df$V1)),]              
          shinyjs::showElement("Representatives")
          output$Representatives <<- DT::renderDataTable({
            for (i in 1:length(df_centroid$Representative)){
              if (!is.na((grep("LSM", df_centroid$Representative[[i]])) && (grep("LSM", df_centroid$Representative[[i]])==1))){
                df_centroid$Representative[[i]] <- makeLink(df_centroid$Representative[[i]])
              }
              else if(!is.null(grep("LSM", df_centroid$Representative[[i]]))){
                df_centroid$Representative[[i]] <- gsub("ZINC0", "ZINC", df_centroid$Representative[[i]])
                df_centroid$Representative[[i]] <- makeLink2(df_centroid$Representative[[i]])
              }
            }
            colnames(df_centroid) <-c("Representative LSM ID", "Representative Name", "Cluster")
            
            ###################### Addded 11-13-2020
            #df_centroid <- df_centroid[-which(duplicated(df_centroid$Cluster)),]
            ######################
            
            return(df_centroid)
            #df_centroid <- df_centroid[,-2]
            
            df_centroid <<- df_centroid
          }, escape = FALSE, rownames = FALSE)
          #print(df_centroid$Representative)
          df_centroid <<- df_centroid
          output$downloadReps <- renderUI({
            downloadButton("RepDownload", label = "Representatives")
          })
          
          print(paste("Got the centroids as ", centroid[[1]], sep=""))
          incProgress(1/2, detail = "Plotting MDS")
          
          source("lib/MDS.R")
          #source("/srv/shiny-server/lib/MDS.R")
          #7. Generate MDS plot from representatives distance to one another, pie chart radius corresponds to cluster size
          if (length(centroid) < 3){
            shinyjs::showElement("check_clusternumber")
            output$check_clusternumber <- renderUI({
              shinyjs::hideElement("distPlot")
              #shinyjs::hideElement("MDSPlot")
              shinyjs::hideElement("Representatives")
              p("No clusters of that size.  Please enter a smaller minimum cluster size.")
              
              
            })
            return()
          }
          MDS_plot(centroid)
          #            centroid_clusters(centroid, ClusterMembers)
          #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(ClusterCentroid), hjust=1, vjust=2))
          #                             + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
          shinyjs::showElement("MDSPlot")
          output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=x, y=y, label=unlist(ClusterCentroid), hjust=1, vjust=2))
                                       + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "My Candidates"))  
                                       + scale_fill_manual(values=c("green", "blue")) +
                                         ggtitle("Multidimentional Scaling Plot of LINCS Candidates"))
          #output$MDSPlot<<- renderPlot(ggplot() + geom_text(mapping = aes(x=cluster_info$x1, y=cluster_info$y1, label=unlist(centroid), hjust=1, vjust=2))
          #                           + geom_scatterpie(mapping = aes(x=x1, y=y1, r=radius), data=cluster_info, cols=c("LINCS", "Added"))) #"P2Count", "P4Count", "LINCS"))
          incProgress(2/2, detail = "Complete!")
          
        }) # End MDS  
        
      } # End clustering - no added compounds
      output$downloadReps <- renderUI({
        downloadButton("RepDownload", label = "Representatives")
      })
      
      output$downloadClusters <- renderUI({
        downloadButton("ClustersDownload", label = "Clusters")
      })
    }) # Clustering Progress Bar
    
  }
  ) #Cluster ObserveEvent  
  
  # output$downloadCandidates <- renderUI({
  #   if((input$AddedCompounds == 1))
  #   downloadButton("CandidatesDownload", label = "My Candidates Ranked")
  # })
  
  output$CandidatesDownload <- downloadHandler(
    filename = function() {
      'my_candidates.csv'
    },
    content = function(filename){
      if(exists("test_bypass_minsim")){
        #write.csv(test_bypass_minsim[order(test_bypass_minsim[,4], test_bypass_minsim[,3], decreasing=TRUE),], filename, row.names=FALSE)
        write.csv(test_bypass_minsim[order(test_bypass_minsim[,6], test_bypass_minsim[,5], decreasing=TRUE),], filename, row.names=FALSE)
      }
      else{
        #write.csv(test_bypass_fpsim[order(test_bypass_fpsim[,4], test_bypass_fpsim[,3], decreasing=TRUE),], filename, row.names=FALSE)
        write.csv(test_bypass_fpsim[order(test_bypass_fpsim[,6], test_bypass_fpsim[,5], decreasing=TRUE),], filename, row.names=FALSE)
      }
    }
  )
  
  output$LINCSDownload <- downloadHandler(
    filename = function() {
      'lincs_candidates.csv'
    },
    content = function(filename){
      write.csv(lincs_output[,-5], filename, row.names = FALSE)
    }
  )

  output$AnalogsDownload <- downloadHandler(
    filename = function() {
      'table_download.csv'
    },
    content = function(filename){
      if(exists("test_bypass_minsim")){
        write.csv(test_bypass_minsim, filename, row.names=FALSE)
      }
      else{
        write.csv(test_bypass_fpsim, filename, row.names=FALSE)
      }
    })

  
  output$RepDownload <- downloadHandler(
    filename = "representatives.csv",
    content = function(filename){
      write.csv(df_centroid, filename, row.names=FALSE)
    }
  )
  
  output$ClustersDownload <- downloadHandler(
    filename = "clusters.csv",
    content = function(filename){
      write.csv(ClusterMembers, filename)
    }
  )
  
 
  #####Set it up so you can download Clusters the Representatives are part of
  
  output$RepDownload <- downloadHandler(
    filename = "representatives.csv",
    content = function(filename){
      write.csv(df_centroid, filename)
    }
    
  )
  output$ClusterDownload <- downloadHandler(
    filename = "clusters.csv",
    content = function(filename){
      write.csv(ClusterMembers, filename)
    }
  )
  
  output$Max_Scores<- downloadHandler(
    filename = "max_scores.csv",
    content = function(filename){
      write.csv(combined_score_output, filename)
    }
    
  )  
  # output Combined Scores
  #  source("./lib/Combined_Score.R")
  #  combined_score_output <<- combined_score(simMA, max_cons_2)
  #  output$Combined_Score <- DT::renderDataTable({
  #    
  #    return(combined_score_output)
  #df_concordances_o$`LSM-ID`<- makeLink(df_concordances_o$`LSM-ID`)
  #return(display)
  #return(df_concordances_o)
  #  },
  # escape = FALSE, rownames = FALSE)
  #  output$Similar_NCI <- eventReactive(input$GetNCI, {
  observeEvent(input$GetNCI, {
    print("Finding similar NCI compounds")
    source("./lib/SimilaritySearch.R")
    nci_list <- list()
    for(i in 1:length(centroid)){
      similarity_search(centroid[i])
      nci_list[[i]] <- output_nci}
    #nci_df <- as.data.frame(cbind(unlist(nci_list, recursive = FALSE), unlist(ClusterCentroid)))
    nci_df <- mapply(c,ClusterCentroid, nci_list)
    nci_df2 <- list()
    for(i in 1:length(nci_df)){
      if(is.character(nci_df[[i]])){ 
        nci_df2 <- list.append(nci_df2, nci_df[[i]])}
    }
    
    nci_df3 <- data.frame()
    for (j in 1:length(nci_df2)){
      for (k in 2:length(nci_df2[[j]])){
        cluster <- nci_df2[[j]][1]
        nsc <- nci_df2[[j]][k]
        cmpd_row <- cbind(nsc, cluster)
        nci_df3 <- rbind(nci_df3, cmpd_row)
      }
    } 
    #cluster_number <- unlist(ClusterCentroid)
    #nci_unlist <- unlist(nci_list)
    
    makeLink2 <- function(val) {
      #paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", "https://pubchem.ncbi.nlm.nih.gov/compound/", val, "</a>", sep="")
      paste("<a href='", "https://pubchem.ncbi.nlm.nih.gov/compound/", val,"'>", val, "</a>", sep="")
    }
    output$Similar_NCI <- DT::renderDataTable({
      
      #nci_df3$link <- makeLink2(nci_df3$nsc)
      nci_df3$nsc <- makeLink2(nci_df3$nsc)
      return(nci_df3)
      
    }, escape = FALSE)
    
    
    
    #    output$Similar_NCI <- renderDataTable(nci_df3)
    #print(nci_list)})
  })
  
  
  
  
  
  
  
  
  ###############################################################################################################################################
  ###################### Network Analysis (Tab 3) ###############################################################################################
  ###############################################################################################################################################
  
  
  
  ##############################################################################################################################  
  ################# Global STITCH ##############################################################################################
  ##############################################################################################################################  
  observeEvent(input$GlobalSTITCH, {
    #withProgress(message="Generating Global STITCH", value=0, {
    shinyjs::hideElement("STITCHPlot")
    shinyjs::hideElement("check_clustering")
    shinyjs::showElement("mSTITCHPlot")
    if (!(exists("simMA"))){
      shinyjs::showElement("check_clustering")
      output$check_clustering <- renderUI({
        p("Please perform chemical similarity analysis before network analysis.")
      })
      return()
    }
    print("Starting mSTITCH")
    source("./lib/mStitchFn.R")
    #source("/srv/shiny-server/lib/mStitchFn.R")
    withProgress(message = "Generating Global STITCH.  This may take a while...", value = 0, {
      Sys.sleep(0.25)
      
      S2L <- data.frame('lincspertid' = rownames(ClusterMembers), 
                        'clusterID' = ClusterMembers[ ,1], stringsAsFactors = FALSE)
      key <- df2
      lincsID <- key %>% select('compound', 'lincspertid') %>% unique()
      tmp <- left_join(S2L, lincsID, by = "lincspertid")
      ind <- which(is.na(tmp$compound), arr.ind = TRUE)
      tmp$compound[ind] <- tmp$lincspertid[ind]
      
      clusterList <- list()
      for(i in 1:length(unique(tmp$clusterID))){
        clusterList[i] <- tmp %>% filter(clusterID == i) %>% select(compound)
      }
      
      numbers <- as.numeric()
      for (i in 1:length(clusterList)){
        numbers[i] <- length(clusterList[[i]])
      }
      
      clus_size <- as.numeric()
      if(!(input$all_clusters)){
        clus_size <- 3
      } else {
        clus_size <- 1
      }
      
      ind <- as.data.frame(numbers) %>% dplyr::filter(numbers >= clus_size) %>% nrow()
      clusterList_tom <- clusterList[1:ind]
      
      chemX <- unlist(clusterList_tom)
      genes <- input$gene_knockdown
      
      #incProgress(0.6, message = "Gethering global gene-chemical interactions...Slow")
      
      stitch_all <- lapply(clusterList_tom, mStitch, gene = genes)
      doc <- bind_rows(stitch_all)  
      
      tmp <- strsplit(doc$X15, "|", fixed = TRUE) %>% sapply("[", 1) %>% 
        strsplit(":", fixed = TRUE) %>% sapply("[", 2) %>% 
        as.numeric()
      interScore <- data.frame("from" = doc$X3, "to" = doc$X4, "score" = tmp, "width" = (tmp*2)^2, 
                               "title" = paste0("<p>", "Score:", "<br>", tmp, "</p>"), stringsAsFactors = FALSE)
      #interScore <<- interScore
      
      mol <- unique(c(interScore$from, interScore$to)) 
      group <- toupper(mol) %>% replace(. %in% toupper(chemX), "Candidates in STITCH") %>%
        replace(. %in% toupper(genes), "Target Gene") %>%
        replace(!. %in% c("Candidates in STITCH", "Target Gene"), "STITCH gene/compound")
      
      link_modnodes <- character()
      for (i in seq_along(mol)) for (j in seq_along(chemX)) {
        if(grepl(toupper(genes), toupper(mol[i])) ){
          link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=', mol[i], '>', mol[i], '</a>')
        } else if(grepl("CHEMBL", mol[i]) ){
          link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://www.ebi.ac.uk/chembl/beta/compound_report_card/', mol[i], '>', mol[i], '</a>')
        } else if(grepl("ZINC", mol[i]) ){
          link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'http://zinc15.docking.org/substances/',  mol[i], '>', mol[i], '</a>')
        } else if(grepl("LSM", mol[i]) ){
          link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'http://lincsportal.ccs.miami.edu/SmallMolecules/view/', mol[i], '>', mol[i], '</a>')
        } else if(toupper(mol[i]) %in% toupper(chemX[j]) ){
          link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://pubchem.ncbi.nlm.nih.gov/search/#query=', mol[i], '>', mol[i], '</a>')
        } else if( is.na(link_modnodes[i]) ){
          link_modnodes[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://en.wikipedia.org/wiki/', mol[i], '>', mol[i], '</a>')
        }
      }
      
      molnodes <- data.frame(label = mol, id = mol, group = group, title = link_modnodes)
      
      # get stitch undetected compounds (caused by whatever reasons i.e., name/id confusing between chemInput and stitch) and their landing pages 
      ind <- toupper(chemX) %in% toupper(mol)
      stitch_undetect <- chemX[!ind]
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
      
      # when need unmatched inputs posted in the network, adding undetected compounds into molnodes
      molnodes <- rbind(molnodes, undetect_Node)
      molnodes <- unique(molnodes)
      # molnodes <<- molnodes
      
      #incProgress(0.2, message = "Generating global network")
      Sys.sleep(0.5)
      
      output$numberClust <- renderText(paste0("Global STITCH network plotted from ", length(clusterList_tom), 
                                              " clusters (>=", clus_size, " compounds/cluster) (",  ceiling((100*length(clusterList_tom)/length(clusterList))), 
                                              "%)"))
      
      
      output$mSTITCHPlot <- renderVisNetwork(visNetwork(molnodes, interScore) %>%
                                               visNodes(physics = TRUE) %>%
                                               visEdges(physics= TRUE, smooth = TRUE, color = list(color = "#008080", highlight = "blue")) %>% 
                                               visGroups(groupname = "Candidates in STITCH", shape = "ellipse",
                                                         #font = list(color = "#660000"), #"133f08"
                                                         font = list(color = "black"),
                                                         #color = list(background = "ffa7a7", border = "black", #"aefc6a" 
                                                         color = list(background = "magenta", border = "black", #"aefc6a" 
                                                                      highlight = list(background = "ffa7a7", 
                                                                                       border = "red"))) %>%
                                               visGroups(groupname = "Target Gene", shape = "circle",
                                                         #font =  list(color = "blue"),
                                                         font =  list(color = "black"),
                                                         #color = list(background = "fdff00", border = "black",
                                                         color = list(background = "yellow", border = "black",
                                                                      highlight = list(background = "feff8f", 
                                                                                       border = "red")) ) %>%
                                               visGroups(groupname = "STITCH gene/compound", shape = "dot", size = 10,
                                                         #font = list(color = "magenta", size = 10),
                                                         font = list(color = "black", size = 10),
                                                         #color = list(background = "ddd" , border = "black", #background ="#bb84ff" 
                                                         #             highlight = list(background = "caff44", 
                                                         color = list(background = "white" , border = "black", #background ="#bb84ff" 
                                                                      highlight = list(background = "caff44", 
                                                                                       border = "red"))) %>%
                                               visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                                                         #color = list(background = "82a9ff", border = "253085",  #4e91fc
                                                         font = list(color = "white", size = 10),
                                                         color = list(background = "blue", border = "253085",  #4e91fc
                                                                      highlight = list(background = "4e91fc", 
                                                                                       border = "red"))) %>%
                                               # visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
                                               #           font = list(color = "660000"), 
                                               #           color = list(background = "ffa7a7", border = "black", 
                                               #                        highlight = list(background = "ffa7a7", 
                                               #                                         border = "red"))) %>%
                                               # visGroups(groupname = "Gene KD", shape = "circle",
                                               #           font =  list(color = "blue"),
                                               #           color = list(background = "fdff00", border = "black",
                                               #                        highlight = list(background = "feff8f", 
                                               #                                         border = "red")) ) %>%
                                               # visGroups(groupname = "STITCH gene/compound", shape = "dot", size = 10,
                                             #           font = list(color = "magenta", size = 10),
                                             #           color = list(background = "ddd" , border = "black", 
                                             #                        highlight = list(background = "caff44", 
                                             #                                         border = "red"))) %>%
                                             # visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                                             #           color = list(background = "82a9ff", border = "253085",  
                                             #                        highlight = list(background = "4e91fc", 
                                             #                                         border = "red"))) %>%
                                             visIgraphLayout(randomSeed = 9, layout = "layout_components", physics = FALSE, smooth = FALSE) %>% 
                                               visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                                               visInteraction(navigationButtons = TRUE) %>%
                                               visLegend(width = 0.25) %>%
                                               visExport(type="png", name="export-network", float="left") # note: visExport only works under Shiny environment
      )
      incProgress(0.1, message = "Finish!")
      Sys.sleep(0.5)
      
    })
    #} 
    
  })    
  
  ######Not yet working, df2 is not a global variable yet.
  #     observeEvent(input$GetSTITCH, {
  #     print("Getting Script")
  #     #source("./lib/StitchNetFn3.R")
  #     source("./lib/StitchNetFn7.R")
  #     print("sourced")
  #     #stitchNet("sirolimus", "mtor")
  #     SelectedCluster <- input$ClusterNumber
  #     h <- 1
  #     CompoundLSM <- list()
  #     for (i in 1:length(ClusterMembers)){
  # #      if (ClusterMembers[[i]] %in% LINCSCompounds)
  #       if (ClusterMembers[[i]] == SelectedCluster) {
  #         CompoundLSM[[h]] <- labels(ClusterMembers)[[1]][[i]]
  #         h <- h+1
  #       }
  #     }
  #     CompoundLSM <<- CompoundLSM
  #     print("IDs obtained")
  #     h <- 1
  # 
  #     CompoundName <- list()
  #     for (j in 1:length(CompoundLSM)){
  #       if (grep("LSM-", CompoundLSM[[j]])==1){
  #         #for (f in 1:length(df2$lincspertid)){
  #           if (CompoundLSM[[j]] %in% df2$lincspertid){
  #             index <- which(df2$lincspertid==CompoundLSM[[j]])
  #             if (length(index)>1){
  #               index <- index[[1]]
  #             }
  #             #CompoundName[[h]] <- df2$stitchID[[index]]
  #             CompoundName[[h]] <- df2$compound[[index]]
  #             h <- h+1
  #         #  }
  #         }
  #       }
  #       else if (!("LSM" %in% CompoundLSM[[j]])){
  #         CompoundName[[h]] <- CompoundLSM[[j]]
  #         h <- h+1
  #       }
  #       
  #     }
  #     CompoundName <<- CompoundName
  #     
  #     chemInput <- unlist(CompoundName)
  #     stitchNet(chem = chemInput, gene = input$Gene, limit = input$Connections)
  #     print("function ran")
  #     output$STITCHPlot <- renderVisNetwork(visNetwork(molnodes, interScore) %>% 
  #                                             visNodes(physics = FALSE) %>%
  #                                             visEdges(physics= FALSE, smooth = FALSE, color = list(color = "gray", highlight = "blue")) %>%
  #                                             visGroups(groupname = "Input chemicals present in STITCH", shape = "ellipse",
  #                                                       color = list(background = "ffa7a7", border = "black", 
  #                                                                    highlight = list(background = "ffa7a7", 
  #                                                                                     border = "red"))) %>%
  #                                             visGroups(groupname = "Gene KD", shape = "circle",
  #                                                       color = list(background = "fdff00", border = "black",
  #                                                                    highlight = list(background = "feff8f", 
  #                                                                                     border = "red")) ) %>%
  #                                             visGroups(groupname = "STITCH connectors", shape = "dot", size = 10,
  #                                                       color = list(background = "ddd", border = "black", 
  #                                                                    highlight = list(background = "caff44", 
  #                                                                                     border = "red"))) %>%
  #                                             visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
  #                                                       color = list(background = "4e91fc", border = "253085", 
  #                                                                    highlight = list(background = "4e91fc", 
  #                                                                                     border = "red"))) %>%
  #                                             visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = FALSE) %>%
  #                                             visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
  #                                             visInteraction(navigationButtons = TRUE) %>%
  #                                             visLegend(width = 0.25) %>%
  #                                             visExport(type="png", name="export-network", float="left")
  #                                           )
  #     print("finished")
  # 
  #   }) #End of STITCH
  
  
  
  ############################################################################################################################
  # Revised STITCH ###########################################################################################################
  ############################################################################################################################
  observeEvent(input$GetSTITCH, {
    withProgress(message="Generating STITCH Network", value=0, {
    shinyjs::hideElement("check_clustering")
    shinyjs::hideElement("mSTITCHPlot")
    shinyjs::showElement("STITCHPlot")
    if (!(exists("simMA"))){
      shinyjs::showElement("check_clustering")
      output$check_clustering <- renderUI({
        p("Please perform chemical similarity analysis before network analysis.")
      })
      return()
    }
    print("Getting Script")
    #source("./lib/StitchNetFn3.R")
    #source("./lib/StitchNetFn9.R")
    source("./lib/StitchNetFn9_updated.R")
    #source("/srv/shiny-server/lib/StitchNetFn9.R")
    print("sourced")
    #stitchNet("sirolimus", "mtor")
    SelectedCluster <- input$ClusterNumber
    
    #######test#
    #SelectedCluster <- 1
    ##################
    
    h <- 1
    CompoundLSM <- list()
    for (i in 1:length(ClusterMembers)){
      #      if (ClusterMembers[[i]] %in% LINCSCompounds)
      if (ClusterMembers[[i]] == SelectedCluster) {
        CompoundLSM[[h]] <- labels(ClusterMembers)[[1]][[i]]
        h <- h+1
      }
    }
    CompoundLSM <<- CompoundLSM
    print("IDs obtained")
    h <- 1
    
    CompoundName <- list()
    for (j in 1:length(CompoundLSM)){
      # if (grep("LSM-", CompoundLSM[[j]])==1){
      
      ## add this ###
      if (length(grep("LSM-", CompoundLSM[[j]]))!=0 && grep("LSM-", CompoundLSM[[j]])==1){
        ###############
        
        #for (f in 1:length(df2$lincspertid)){
        if (CompoundLSM[[j]] %in% df2$lincspertid){
          index <- which(df2$lincspertid==CompoundLSM[[j]])
          if (length(index)>1){
            index <- index[[1]]
          }
          #CompoundName[[h]] <- df2$stitchID[[index]]
          CompoundName[[h]] <- df2$compound[[index]]
          h <- h+1
          #  }
        }
      }
      else if (!("LSM" %in% CompoundLSM[[j]])){
        CompoundName[[h]] <- CompoundLSM[[j]]
        h <- h+1
      }
      
    }
    
    ######## add this script to server ######
    if(length(CompoundName) == 0){
      CompoundName <- CompoundLSM
    }
    #########################################
    
    CompoundName <<- CompoundName
    
    chemInput <- unlist(CompoundName)
    #stitchNet(chem = chemInput, gene = input$Gene, limit = input$Connections)
    
    ## add this ########################        
    chemInput <- make.unique(chemInput)
    chemInput <<- chemInput
    stitchNet(chem = chemInput, gene = input$gene_knockdown, limit = input$Connections, confidence = input$confidence)
    #########test
    
    #chemInput_2 <- make.unique(chemInput)
    #stitchNet(chem = chemInput_2, gene = "EGFR", limit = 5, confidence = 400)
    
    #################################
    
    print("function ran")
    
    ##### add/replace this script to server #######################################    
    observeEvent(input$show_undetectNodes, {
      if(!input$show_undetectNodes){
        select_molnodes <- molnodes %>% filter(group != "Undetected by STITCH")
        assign("select_molnodes", select_molnodes, envir = .GlobalEnv)
      } else if(input$show_undetectNodes){
        select_molnodes <- molnodes 
        assign("select_molnodes", select_molnodes, envir = .GlobalEnv)
      }
    })
    
    # output$STITCHPlot <- renderVisNetwork(stitchNet(chem = chemInput, gene = input$Gene, limit = input$Connections))
    
    output$STITCHPlot <- renderVisNetwork(visNetwork(select_molnodes, interScore) %>%
                                            visNodes(physics = !(input$fixedNet)) %>%
                                            visEdges(physics = !(input$fixedNet), smooth = TRUE, color = list(color = "#008080", highlight = "blue")) %>% #"grey
                                            visGroups(groupname = "Candidates in STITCH", shape = "ellipse",
                                                      #font = list(color = "#660000"), #"133f08"
                                                      font = list(color = "black"),
                                                      #color = list(background = "ffa7a7", border = "black", #"aefc6a" 
                                                      color = list(background = "magenta", border = "black", #"aefc6a" 
                                                                   highlight = list(background = "ffa7a7", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Target Gene", shape = "circle",
                                                      #font =  list(color = "blue"),
                                                      font =  list(color = "black"),
                                                      #color = list(background = "fdff00", border = "black",
                                                      color = list(background = "yellow", border = "black",
                                                                   highlight = list(background = "feff8f", 
                                                                                    border = "red")) ) %>%
                                            visGroups(groupname = "STITCH gene/compound", shape = "dot", size = 10,
                                                      #font = list(color = "magenta", size = 10),
                                                      font = list(color = "black", size = 10),
                                                      #color = list(background = "ddd" , border = "black", #background ="#bb84ff" 
                                                      #             highlight = list(background = "caff44", 
                                                      color = list(background = "white" , border = "black", #background ="#bb84ff" 
                                                                   highlight = list(background = "caff44", 
                                                                                    border = "red"))) %>%
                                            visGroups(groupname = "Undetected by STITCH", shape = "ellipse", 
                                                      #color = list(background = "82a9ff", border = "253085",  #4e91fc
                                                      font = list(color = "white", size = 10),
                                                      color = list(background = "blue", border = "253085",  #4e91fc
                                                                   highlight = list(background = "4e91fc", 
                                                                                    border = "red"))) %>%
                                            visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = !(input$fixedNet), smooth = TRUE) %>% #"layout_nicely"
                                            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                                            visInteraction(navigationButtons = TRUE) %>%
                                            visLegend(width = 0.25) #%>%
                                          #visExport only works under Shiny environment, so enable a code below when running in Shiny environment
                                          #visExport(type="pdf", name="export-network", label = "Save network", float="left", background = "#fff")
    )
    ###########################################################################    
    
    
    
   
     print("Finished")
    })  
  })  
  
########
    
  }) # End of app
