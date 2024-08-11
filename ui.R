#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
# #########################################################
# # Install Bioconductor Packages
# ###############################
# 
# list.of.packages_bc <- c("ChemmineOB", "ChemmineR")
# new.packages_bc <- list.of.packages_bc[!(list.of.packages_bc %in% installed.packages()[,"Package"])]
# 
# if(length(new.packages_bc)){ 
#   
#   if (!requireNamespace("BiocManager", quietly=TRUE)){
#     BiocManager::install(version="3.11")
#     BiocManager::install(new.packages_bc, ask=FALSE)}
# }

#######################################################################
# Install all required packages
###############################

list.of.packages <- c(
  #"shiny",
  "shinyjs",
  "ggplot2",
  "plotly",
  "httr",
  "jsonlite",
  "DT",
  "visNetwork",
  #"ChemmineOB",
  #"ChemmineR",
  "plyr",
  "dendextend",
  "colorspace",
  "ggforce",
  "rlist",
  "scatterpie",
  "ggrepel",
  "bazar",
  "XML",
  "RCurl",
  "bitops",
  "scrapeR",
  "igraph",
  "circlize",
  "enrichR",
  "readr",
  "dplyr",
  "gplots"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=TRUE, INSTALL_opts = '--no-lock')

############## Bioconductor ################################################################
# 

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version="3.12")

list.of.packages_bm <- "BiocManager"
new.packages_bm <- list.of.packages_bm[!(list.of.packages_bm %in% installed.packages()[,"Package"])]
if(length(new.packages_bm)) install.packages(new.packages_bm, dependencies=TRUE, INSTALL_opts = '--no-lock', version="3.12")


list.of.packages_bc <- c("ChemmineOB", "ChemmineR")
new.packages_bc <- list.of.packages_bc[!(list.of.packages_bc %in% installed.packages()[,"Package"])]

if(length(new.packages_bc)){
  
  BiocManager::install(new.packages_bc, ask=FALSE)
}

##########################################################
# Load Libraries
#################
library(shiny)
library(shinyjs)
library(plotly)
library(visNetwork)

# Define UI for application that draws a histogram

shinyUI(
  fluidPage(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "customize.css")
    ),
  # Application title
  #titlePanel("Sig2Lead"),
  titlePanel(""),
      fluidPage(
        #useShinyjs(),
        #        img(src='s2l.png', align='right', width=100,height=100),
                
        tabsetPanel(id = "tabs",
          tabPanel("Signature Connectivity Analysis", value ="search", icon = NULL,
                   # fluidRow(
############# Add this ##############################################     
                 br(),
                  #tags$h1("This tab permits ranking of potential inhibitors of a target gene by employing transcriptional signature connectivity and chemical similarity analysis"),
                  tags$h1("Sig2Lead ranks potential inhibitors of a target gene by employing transcriptional signature connectivity and chemical similarity analysis"),
                  br(),
                  br(),
                     #fluidRow(
                     #column(3, 
                    #        wellPanel(
                     sidebarPanel(position="left",        
#####################################################################                   
                     #br(),
                            selectInput("Signature", NULL, list("Define Target Gene", "Upload a Signature", "Find Analogs in LINCS"), selected = "Input a Gene"),
                     br(),
                     
                     conditionalPanel(
                       condition = "input.Signature == 'Define Target Gene'",
                            textInput("gene_knockdown", "Input a Target Gene", value="egfr")
                     ),
                     conditionalPanel(
                        condition = "input.Signature == 'Upload a Signature'",
                            #useShinyjs(),
                            fileInput("UploadSignature", "Upload Signature", multiple = TRUE, accept = c(".txt", ".csv"), width = NULL, 
                                         buttonLabel = "Browse...", placeholder = "No file selected")
                            #actionButton('reset', 'Reset')
                     ),
                     conditionalPanel(
                        condition =  "input.Signature != 'NULL'",
                            fileInput("AddedCompounds", "Add candidate compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
                                         buttonLabel = "Browse...", placeholder = "No file selected")
                     ),
                            br(),
                     #conditionalPanel(
                    #   condition =  "input.Signature != 'Upload a Signature'",
                    #fileInput("AddedCompounds", "Add candidate compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
                    #           buttonLabel = "Browse...", placeholder = "No file selected")
                    # ),
                    conditionalPanel(
                        condition =  "input.Signature != 'Upload a Signature'",
                        downloadButton("ExampleFile", label = "Download Example SDF File")
                    ),
                    conditionalPanel(
                        condition =  "input.Signature != 'Upload a Signature'",
                        downloadButton("ExampleFileSmiles", label = "Download Example SMILES File")
                    ),
                    conditionalPanel(
                        condition =  "input.Signature == 'Upload a Signature'",
                        downloadButton("ExampleSig", label = "Download Example Signature File")
                    ),
                   uiOutput("mytest"),
                    # conditionalPanel(
                    #      condition = "input.AddedCompounds.name == null",
                    #         textInput("AddedLabel", "Label for added compounds", value="Added"),
                    #          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
                    # ),
                            br(),
                   #selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
                   
                  
                  br(),
                  actionButton(label="Show Advanced Options", "Advanced_button"),
                  br(),
                  br(),

                  
                   #checkboxInput("Advanced", "Advanced Options", selected),
                   conditionalPanel(condition="input.Advanced_button%2 ==  1 && input.Signature == 'Define Target Gene'",
                        #numericInput("topN", "Number of Analogs", value=1, min=1,max=5,step=1),
                        selectInput("ILINCS", "iLINCS Version", list("Current", "Legacy"), selected="Current"),
                        selectInput("ConOrDiscon","Activation or Inhibition", list("Inhibit", "Activate"), selected = "Inhibit"),            
                        conditionalPanel(condition="input.ConOrDiscon == 'Inhibit'",
                          selectInput("Algorithm","Chemical Similarity Search", list("minSim", "fpSim"), selected = "minSim"),
                          numericInput("Concordance", "Concordance Threshold", value=0.1, min=0.1, max=0.3, step=0.1)),
                        conditionalPanel(condition="input.ConOrDiscon == 'Activate'",
                                         selectInput("Algorithm","Chemical Similarity Search", list("minSim", "fpSim"), selected = "minSim"),
                                         numericInput("Concordance", "Discordance Threshold", value=-0.1, min=-0.5, max=-0.2, step=-0.1))
                        
                   ),  
                   conditionalPanel(condition="input.Advanced_button%2 ==  1 && input.Signature == 'Upload a Signature'",
                        #numericInput("topN", "Number of Analogs", value=1, min=1,max=5,step=1),
                        #selectInput("ILINCS", "iLINCS Version", list("Current", "Legacy"), selected="Current"),
                        selectInput("ConOrDiscon_2","Activation or Inhibition", list("Concordant", "Discordant"), selected="Concordant"),            
                        conditionalPanel(condition="input.ConOrDiscon_2 == 'Concordant'", numericInput("Concordance_2", "Concordance Threshold", value=0.1, min=0.1, max=0.3, step=0.1)),
                        conditionalPanel(condition="input.ConOrDiscon_2 == 'Discordant'", numericInput("Concordance_3", "Discordance Threshold", value=-0.1, min=-0.5, max=-0.1, step=0.1))
                  ),  
                  conditionalPanel(condition="input.Advanced_button%2 ==  1 && input.Signature == 'Find Analogs in LINCS'",
                    numericInput("topN", "Number of Analogs", value=1, min=1,max=5,step=1), 
                    selectInput("Algorithm","Chemical Similarity Search", list("minSim", "fpSim"), selected = "minSim"),
                    numericInput("Similarity", "Chemical Similarity Threshold", value=0.65, min=0.2, max=1.0, step=0.05)),
                    
                  br(),
                  # conditionalPanel(
                  #     condition = "if(!is.null(input.AddedCompounds))",
                  #     checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
                  #),     
                   #br(),
                            actionButton(label="Go!", "Go"),

                  br(),
                  br(),
                  br(),
                  # output$downloadCandidates <- renderUI({
                  #   req("input.Go==1")
                  #   downloadButton("CandidatesDownload", label = "My Candidates Ranked")
                  #   }),
                  uiOutput("downloadCandidates"),
                  #br(),
                  uiOutput("downloadLINCS"),
                  uiOutput("downloadUpload"),
                  uiOutput("downloadAnalogs"),
      
                  #downloadButton("ExampleFile", label = "Download Example Candidate File"),
                  #downloadButton("LINCSDownload", label="LINCS Candidates Ranked"),
                  
                  #downloadButton("SMILESDownload", label = "Download SMILES"),
                  #downloadButton("SDFDownload", label = "Download SDF"),
                  #downloadButton("ConTableDownload", label="Download Table"),
                  #downloadButton("NCIDownload", label = "Download Similar NCI")
                  
                  #tags$div(
                  #  tags$div(class='tlp',
                  #           icon("question-circle", lib = "font-awesome"),
                  #           tags$span("Learn More", id = "span1"),
                  #           tags$div(class="tooltip_text",
                  #                    tags$h3(id = "hover1", "For a set of candidate molecules defined in a user 
                  #                                            provided file and a target gene specified by the user, 
                  #                                            candidate drugs are ranked by similarity to their 'concordant' 
                  #                                            LINCS analogs that share transcriptional signature similarity 
                  #                                            with a knock-down of the target gene. Note that not all 
                  #                                            genes had their knock-down included in",
                  #                            tags$a(href="https://lincsproject.org/", "LINCS")
                  #                    )
                  #                    
                  #           )
                  #  )
                  #),

######################## Added 4/25/2023 ##############################################
conditionalPanel(
  condition = "input.Signature == 'Define Target Gene'",
  tags$div(
    tags$div(class='tlp',
             icon("question-circle", lib = "font-awesome"),
             tags$span("Learn More", id = "span1"),
             tags$div(class="tooltip_text",
                      tags$h3(id = "hover1", "For a set of candidate molecules defined in a user 
                                                              provided file and a target gene specified by the user, 
                                                              candidate drugs are ranked by similarity to their 'concordant' 
                                                              LINCS analogs that share transcriptional signature similarity 
                                                              with a knock-down of the target gene. Note that not all 
                                                              genes had their knock-down included in",
                              tags$a(href="https://lincsproject.org/", "LINCS.  "),  "  However a user-defined signature can be used as well."
                      )
                      
             )
    )
  ),
),
conditionalPanel(
  condition = "input.Signature == 'Upload a Signature'",
  tags$div(
    tags$div(class='tlp',
             icon("question-circle", lib = "font-awesome"),
             tags$span("Learn More", id = "span1"),
             tags$div(class="tooltip_text",
                      tags$h3(id = "hover1", "For a set of candidate molecules defined in a user 
                                              provided file and a user-provided loss-of-function (LOF) or gain-of-function (GOF) signature, 
                                              candidate drugs are ranked by similarity to their 'concordant' or 'discordant'
                                              LINCS analogs that share transcriptional signature similarity 
                                              with the provided signature. A concordant LINCS analog to a LOF signature 
                                              or a discordant LINCS analog to a GOF signature indicates an antagonist, 
                                              while a discordant LINCS analog to a LOF signature or a concordant LINCS 
                                              analog to a GOF signature indicates an agonist.",
                      )
                      
             )
    )
  ),
),
conditionalPanel(
  condition =  "input.Signature == 'Find Analogs in LINCS'",
  tags$div(
    tags$div(class='tlp',
             icon("question-circle", lib = "font-awesome"),
             tags$span("Learn More", id = "span1"),
             tags$div(class="tooltip_text",
                      tags$h3(id = "hover1", "For a set of candidate molecules defined in a user 
                                              provided file, candidate drugs are ranked by similarity 
                                              to all LINCS analogs." 
                                                              ,                                                              
                              tags$a(href="https://lincsproject.org/", "LINCS")
                      )
                      
             )
    )
  ),
),
                  
                  

############# Add this ######################################################################################                  
                ),

#fluidRow(
#column(9, 
#       htmlOutput("notify1"),
#       visNetworkOutput("KDnet_Plot", width = "900", height = "500px"),
#       htmlOutput("notify2")
       
#)),


# )

#############################################################################################################
      mainPanel(position="right",
                
          br(),
          br(),
          br(),
          br(),
          shinyjs::useShinyjs(),
          #imageOutput("logo"),
          imageOutput("logo2"),
          #uiOutput(outputId="logo"),
          uiOutput("check_connection"),
          uiOutput("sim_search"),
          uiOutput("check_kd"),
          uiOutput("check_sm"),
          
          uiOutput("check_format"),
          br(),
          br(),
          #tags$style(HTML('table.dataTable tr.selected td, table.dataTable td.selected [background-color:pink !important;}')),
          DT::dataTableOutput("CANDIDATES"),
          #downloadButton("CandidatesDownload", label = "Download Candidates"),
          #downloadButton("SMILESDownload", label = "Download SMILES"),
          #downloadButton("SDFDownload", label = "Download SDF"),
          #downloadButton("ExampleFile", label = "Download Example Candidate File"),
          br(),
          br(),
          DT::dataTableOutput("SMILES"),
          br(),
          br(),
          #downloadButton("ConTableDownload", label="Download Table"),
          br(),
          br(),
          DT::dataTableOutput("Similar_NCI"),
          #downloadButton("NCIDownload", label = "Download Similar NCI")
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          verbatimTextOutput("wait_message")
         
          
         
      )     # )
    ), #up to here
        
#        mainPanel(
#        textOutput("ItWorked"),
      #tabPanel("Candidate Compounds", value = "Candidate Compounds", icon =NULL,
               
       #  DT::dataTableOutput("CANDIDATES"),
         #downloadButton("SMILESDownload", label = "Download SMILES"),
         #downloadButton("SDFDownload", label = "Download SDF"),
         #downloadButton("ConTableDownload", label="Download Table")
      #),
      #tabPanel("LINCS Compounds", value = "LINCS Compounds", icon =NULL,
      #    DT::dataTableOutput("SMILES")
      #    downloadButton("SMILESDownload", label = "Download SMILES"),
      #    downloadButton("SDFDownload", label = "Download SDF"),
      #    downloadButton("ConTableDownload", label="Download Table")
      #  ),
#tabPanel("Consensus Ranking", value="Combined", icon=NULL,
         #          fluidRow(
         #            column(3, 
         #                   wellPanel(
         #                     
         #                     #####################################################################                   
         #                     br(),
         #                     selectInput("Signature", NULL, list("Input a Gene", "Upload a Signature", "Similarity Search"), selected = "Input a Gene"),
         #                     br(),
         #                     
         #                     conditionalPanel(
         #                       condition = "input.Signature == 'Input a Gene'",
         #                       textInput("gene_knockdown", "Input a Gene", value="bcl2a1")
         #                     ),
         #                     conditionalPanel(
         #                       condition = "input.Signature == 'Upload a Signature'",
         #                       fileInput("UploadSignature", "Upload Signature", multiple = FALSE, accept = c(".csv"), width = NULL, 
         #                                 buttonLabel = "Browse...", placeholder = "No file selected")
         #                     ),
         #                     br(),
         #                     fileInput("AddedCompounds", "Add compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
         #                               buttonLabel = "Browse...", placeholder = "No file selected"),
         #                     uiOutput("mytest"),
         #                     # conditionalPanel(
         #                     #      condition = "input.AddedCompounds.name == null",
         #                     #         textInput("AddedLabel", "Label for added compounds", value="Added"),
         #                     #          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
         #                     # ),
         #                     br(),
         #                     #selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
         #                     br(),
         #                     selectInput("Advanced", "See Advanced Options", list("Yes", "No"), selected = "No"),
         #                     conditionalPanel(condition="input.Advanced == 'Yes'",
         #                                      selectInput("Algorithm","Algorithm", list("minSim", "fpSim"), selected = "minSim"),
         #                                      selectInput("ConOrDiscon","Activation or Inhibition", list("Concordant", "Discordant"), selected = "Concordant"),
         #                                      numericInput("Concordance", "Concordance Threshold", value=0.2)
         #                     ),        
         #                     #conditionalPanel(
         #                     #     condition = "if(!is.null(input.AddedCompounds))",
         #                     #     checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
         #                     #),     
         #                     #br(),
         #                     actionButton(label="Go!", "Go")
         #                     ############# Add this ######################################################################################                  
         #                   )),
         #            #fluidRow(
         #            column(9, 
         #                   htmlOutput("notify1"),
         #                   visNetworkOutput("KDnet_Plot", width = "900", height = "500px"),
         #                   htmlOutput("notify2")
         #                   
         #            ))
         #          
         #          
         #          # )
         #          
         #          #############################################################################################################
         #          
         #          # )
#), #up to here
#DT::dataTableOutput("Combined_Score"),
#downloadButton("Max_Scores", label = "Download Scores")),
#        #downloadButton("NCIDownload", label = "Download Similar NCI"))

      tabPanel("Clustering and Visualization", value="heatmap", icon = NULL,
               #br(),
               #br(),
               br(),
               h1("	This tab permits clustering of concordant LINCS analogs and candidate compounds to determine classes of putative inhibitors"),
               br(),
               br(),
          sidebarPanel(position="left", width=5,
               
               #plotOutput("distPlot", width="850px", height="800px")
      #          fluidRow(
      #            column(3,
      #                   wellPanel(
      # 
      #                     #####################################################################
      #                     br(),
      #                     selectInput("Signature", NULL, list("Input a Gene", "Upload a Signature", "Similarity Search"), selected = "Input a Gene"),
      #                     br(),
      # 
      #                     conditionalPanel(
      #                       condition = "input.Signature == 'Input a Gene'",
      #                       textInput("gene_knockdown", "Input a Gene", value="bcl2a1")
      #                     ),
      #                     conditionalPanel(
      #                       condition = "input.Signature == 'Upload a Signature'",
      #                       fileInput("UploadSignature", "Upload Signature", multiple = FALSE, accept = c(".csv"), width = NULL,
      #                                 buttonLabel = "Browse...", placeholder = "No file selected")
      #                     ),
      #                     br(),
      #                     fileInput("AddedCompounds", "Add compounds in SMILES or SDF (Optional)", multiple = TRUE, accept = c(".txt", ".smi", ".csv", ".sdf"), width = NULL,
      #                               buttonLabel = "Browse...", placeholder = "No file selected"),
      #                     uiOutput("mytest"),
      #                     # conditionalPanel(
      #                     #      condition = "input.AddedCompounds.name == null",
      #                     #         textInput("AddedLabel", "Label for added compounds", value="Added"),
      #                     #          checkboxInput("Bypass", "Bypass Clustering", value=FALSE)
      #                     # ),
      #                     br(),
      #                     #selectInput("ConOrDiscon",NULL, list("Concordant", "Discordant"), selected = "Concordant"),
      #                     br(),
      #                     selectInput("Advanced", "See Advanced Options", list("Yes", "No"), selected = "No"),
      #                     conditionalPanel(condition="input.Advanced == 'Yes'",
      #                                      selectInput("Algorithm","Algorithm", list("fpSim","minSim"), selected = "fpSim"),
      #                                      selectInput("ConOrDiscon","Activation or Inhibition", list("Concordant", "Discordant"), selected = "Concordant"),
      #                                      numericInput("Concordance", "Concordance Threshold", value=0.2)
      #                     ),
      # 
      #                     actionButton(label="Go!", "Go")
      #                     ############# Add this ######################################################################################
      #                   )),
      #            #fluidRow(
      #            column(9,
      #                   htmlOutput("notify1"),
      #                   visNetworkOutput("KDnet_Plot", width = "900", height = "500px"),
      #                   htmlOutput("notify2")
      # 
      #            ))
      # 
      # 
      #          # )
      # 
      #          #############################################################################################################
      # 
      #          # )
      # ), #up to here
                actionButton("Cluster", label="Run SAR"),
                
                #br(),
                #br(),
#               actionButton(label="GenerateHeatmap", "GenerateHeatmap"),
               #textInput("CutHeight", "Tanimoto Similarity", value="0.75"),
               #textInput("ClusterSize", "Minimum Cluster Size", value="3"),
               #actionButton(label="MultiDimensional Scaling", "GetRepresentatives"),
               #br(),
               #br(),
               #br(),
                  actionButton(label="Show Advanced Options", "Advanced_button_mds"),
                  br(),
                  uiOutput("downloadReps"),
                  br(),
                 #br(),
                  uiOutput("downloadClusters"),
                  br(),
                  conditionalPanel(condition="input.Advanced_button_mds%2 ==  1",
                      numericInput("CutHeight", "Tanimoto Similarity", value="0.75"),
                      numericInput("ClusterSize", "Minimum Cluster Size", value="3"),
                      #checkboxInput("cluster_all_cmpds", value=FALSE, label="Cluster all Compounds"),
                      numericInput("NumberCmpds", "Number of Compounds to Cluster", value="500"),
                  br(),
                  #br(),
                  
                 
               #downloadButton("RepDownload", label="Download Representatives"),
               #downloadButton("ClusterDownload", label = "Download Clusters"),

          ),tags$div(
            tags$div(class='tlp',
                     icon("question-circle", lib = "font-awesome"),
                     tags$span("Learn More", id = "span1"),
                     tags$div(class="tooltip_text",
                              tags$h3(id = "hover1", "Top ranking candidate molecules defined 
                                      in a user provided file as well as their 'concordant' 
                                      LINCS analogs that share transcriptional signature similarity 
                                      with a knock-down of the target gene are clustered using 
                                      chemical similarity to identify structural classes of putative 
                                      drug candidates. Note that this step may take some time 
                                      if the number of top candidates to be included is increased." 
                              )
                              
                     )
            )
          ),),
          
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
  
          column(12,
               #plotOutput("distPlot", width = "850px", height = "800px"),
               #plotOutput("MDSPlot", width = "900px", height = "800px"),
               #uiOutput("cluster_check"),
               
               uiOutput("check_analogs"),
               uiOutput("check_clustersize"),
               uiOutput("check_clusternumber"),
               splitLayout(cellWidths = c("50%", "50%"),
               plotOutput("distPlot", width = "500px", height = "500px"),
               plotOutput("MDSPlot", width = "500px", height = "500px")),
               DT::dataTableOutput("Representatives")
               #downloadButton("RepDownload", label="Download Representatives"),
               #downloadButton("ClusterDownload", label = "Download Clusters"),
               #actionButton(label="Get Related NCI Compounds", "GetNCI")
               
              
          ),
          
     
      ),
      # tabPanel("MDS Plot", value="MDS", icon=NULL,
      #          plotOutput("MDSPlot", width = "900px", height = "800px"),
      #          DT::dataTableOutput("Representatives"),
      #          downloadButton("RepDownload", label="Download Representatives"),
      #          downloadButton("ClusterDownload", label = "Download Clusters"),
      #          actionButton(label="Get Related NCI Compounds", "GetNCI")
      #          ),


#      tabPanel("Similar Compounds", value="NCI", icon=NULL,
       #       DT::dataTableOutput("Similar_NCI"),
         #     downloadButton("NCIDownload", label = "Download Similar NCI")),
#### add this script to UI #### 

tabPanel("Network Analysis", value = "mSTITCH", icon = NULL,
         #fluidRow(
        #   column(3, wellPanel(
          br(),
          h1("This tab utilizes known relationships of gene targets and small molecules to identify networks of LINCS analogs, candidate compounds, and gene targets."),
          br(),
          br(),
         sidebarPanel( 
            selectInput("STITCH", NULL, list("Global STITCH", "Cluster Network STITCH"), multiple=FALSE),
            
            conditionalPanel(
              condition = "input.STITCH == 'Global STITCH'",
              checkboxInput("all_clusters", "Shows all clusters", value = FALSE, width = '400px'),
            actionButton(label="View Global STITCH Network", "GlobalSTITCH")
             ),
            br(),
            br(),
            br(),
            conditionalPanel(
              condition = "input.STITCH == 'Cluster Network STITCH'",
            textInput("ClusterNumber", "Cluster Number", value = 1),
            br(),
            selectInput("confidence", "Confidence Score", choices = c(150, 400, 700, 950), selected = 400),
            br(),
            sliderInput("Connections", "Numbers of Connections", min = 2, max = 25, value = 10, ticks = FALSE),
            br(),
            # textInput("Gene", "Gene of Interest", value = "bcl2a1"),
            checkboxInput("fixedNet", "Fixed node position", value = FALSE, width = '400px'),
            checkboxInput("show_undetectNodes", "Shows all", value = FALSE, width = '400px'),
            actionButton(label="View Selected Cluster Network", "GetSTITCH"),
            
            ), 
            tags$h1("Clusters or putative structural classes of top-ranking candidate molecules and their 'concordant' LINCS analogs are mapped into STiTCH drug-target network of known drug-target relationships."),
           #)),
           #column(9, 
            #        br(),      
            #      textOutput("numberClust"),    
           #)
         #),
         #fluidRow(
          #   column(12, 
           #       visNetworkOutput("mSTITCHPlot", width = "1000px", height = "750px"),
          # )
         ),
         mainPanel( 
           uiOutput("check_clustering"),  
           textOutput("numberClust"),
           visNetworkOutput("mSTITCHPlot", width = "1000px", height = "750px"),
           
           br(),
           br(),
           
           visNetworkOutput("STITCHPlot", width = "900px", height = "600px")
           
           #br(),
           #br()
           
          )
          ),
tabPanel("Help",
         mainPanel(
           #includeHTML("./www/Sig2Lead_v2_user_manual_for_help.htm")
           includeHTML("./www/Sig2Lead_Version_2.3_User_Manual.htm")
         )
)

# tabPanel("Selected Cluster Network", value = "STITCH", icon = NULL,
#          fluidRow(
#            column(2, wellPanel(
#              textInput("ClusterNumber", "Cluster Number", value = 1),
#              br(),
#              selectInput("confidence", "Confidence Score", choices = c(150, 400, 700, 950), selected = 400),
#              br(),
#              sliderInput("Connections", "Numbers of Connections", min = 2, max = 25, value = 10, ticks = FALSE),
#              br(),
#              # textInput("Gene", "Gene of Interest", value = "bcl2a1"),
#              checkboxInput("fixedNet", "Fixed node position", value = FALSE, width = '400px'),
#              checkboxInput("show_undetectNodes", "Shows all", value = FALSE, width = '400px'),
#              actionButton(label="View with STITCH", "GetSTITCH")
#            )),
#            column(9, 
#                   br(),
#                   visNetworkOutput("STITCHPlot", width = "900px", height = "600px")
#            )
#          ))


############################
        )
      )
  )
)


