# Load the necessary libraries
library(shiny)
library(shinydashboard)
library(colourpicker)
library(markdown)
library(shinymanager)
library(BiocManager)
library(dplyr)
options(repos = BiocManager::repositories())
# library(conflicted)
# setwd("E:/B-Cell-Lymphoma/ShinyAppPublication/OnkoDrugScreening")
source("./dataLoading.r")
source("./drug_comb_analysis.r")
source("./plot_functions.r")
## deal with the conflict function issue.
# conflict_prefer("select", "dplyr")
# conflict_prefer("rename", "dplyr")
# conflict_prefer("intersect", "base")
convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  mi
}
# Define UI components separately
sidebar <- dashboardSidebar(
  width = 300,
  tags$head(
    tags$style(HTML("
      .file-input { margin-bottom: 0px; } #tab {margin-top:-30px;}
      .select-input { margin-bottom: 2px; }
      .slider-input { margin-bottom: 2px; }
      .action-button { margin-bottom: 2px; }
      .numeric-input { margin-bottom: 2px; }
      .form-group { margin-bottom: 0px; margin-top: 0px;}
      
    "))
  ),

  fileInput("target_primary_data", "Target cancer data (CSV/MAF):",
            multiple = FALSE,
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv",
              ".maf",
              ".maf.gz",
              ".vcf"
            )),
  fileInput("comparison_primary_data", "Comparison cancer data (CSV/MAF):",
            multiple = FALSE,
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv",
              ".maf",
              ".maf.gz",
              ".vcf"
            )),
  fileInput("cellline_data", "Cellline data (CSV/MAF):",multiple = FALSE,
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv",
              ".maf",
              ".maf.gz",
              ".vcf"
            )),
  selectInput("cancer_selector", 
              "Select Target Disease:", 
              choices = unique(DrugClassified_DB$Cancer)),
  sliderInput("Percentage_slider", "Percentage Threshold", min = 0, 
              max = 100, value = 20),
  # selectInput("Testtype_selecter", "Select Test Type",
  #             choices = c("greater", "two.sided", "less"),
  #             selected = "greater"), # Default selection set to 'greater'
  fluidRow(
    column(width = 6,
           selectInput("Testtype_selecter", "Select Test Type",
                       choices = c("greater", "two.sided", "less"),
                       selected = "greater"), # Default selection set to 'greater'
           
    ),
    column(width = 6, 
           selectInput("Response_type_selector", "Select Response",
                       choices = c("Positive","Negative","All"),
                       selected = "All",
                       multiple = FALSE
    )
  )
  ),
  
  tags$label("Image Output Resolution Parameters"),
  fluidRow(
    column(width = 6,
           numericInput("img_w", "width:", 1980, min = 100, max = 4096)
           
    ),
    column(width = 6, 
           numericInput("img_h", "height:", 1080, min = 100, max = 4096)
    )
  ),
  fluidRow(
    column(width = 6,
           actionButton("test_data", "Test_Data", icon = icon("arrow-circle-up"))
    ),
    column(width = 6,
           actionButton("submit", "Submit", icon = icon("arrow-circle-right"))
    )
  ),

  
  sidebarMenu(
    ### Volcano plot Item ----------------
    convertMenuItem(menuItem("VolcanoPlot",tabName = "volcanoplot",
                             #icon=icon("volcanoplot"),
                             # select column for pval
                             selectInput("pvalue_col",
                                         "Input column for significance (y axis)",
                                         choices = c("P.value",
                                                     "AdjP.value"
                                         ),
                                         #pval_cols,
                                         multiple = FALSE
                             ),
                             # set pvalue threshold 
                             sliderInput("pvalue_threshold",
                                         "Set significance threshold",
                                         min = 0,
                                         max = 1,
                                         value = .05),
                             
                             # set log2fc threshold 
                             sliderInput("logfc_threshold",
                                         "Select logfc threshold",
                                         min = 0,
                                         max = 2,
                                         value = 0.5,
                                         step = .1),
                             # drug selector menu
                             uiOutput("drugcombs_selector"),
                             
                             # show/hide logfc and pval line
                             checkboxInput("show_pvalue_threshold",
                                           "Show significance threshold line",
                                           value = TRUE),
                             
                             # output ui for axis label inputs
                             uiOutput("y_axis_labeler"),
                             uiOutput("x_axis_labeler"),
                             
                             # show/hide logfc lines
                             checkboxInput("show_logfc_threshold",
                                           "Show logfc threshold line",
                                           value = TRUE)
    ),tabName = "volcanoplot"),
    ### Heatmap Item ----------------
    convertMenuItem(menuItem("Heatmap", tabName = "heatmap",
                             #icon=icon("heatmap"),
                             sliderInput("cellline_threshold",
                                         "Set percentage threshold based on Cellline",
                                         min = 20, 
                                         max = 100, 
                                         post  = " %", 
                                         value = 30),
                             
                             colourInput("col_high", "Select high proportion colour", "red"),
                             colourInput("col_mid", "Select medium proportion colour", "orange"),
                             colourInput("col_low", "Select low proportion colour", "yellow")
    ),tabName = "heatmap"),
    ### CirclePlot Item ----------------
    convertMenuItem(menuItem("CirclePlot",tabName = "circleplot"
    ),tabName = "circleplot"),
    ### AlluvialDiagram Item ----------------
    convertMenuItem(menuItem("AlluvialDiagram",tabName = "alluvialdiagram",
                             #icon=icon("alluvialdiagram"),
                             uiOutput("drug_selector")),tabName = "alluvialdiagram"),
    convertMenuItem(menuItem("UpsetPlot",tabName = "upsetplot"
    ),tabName = "upsetplot"),
    convertMenuItem(menuItem("BarPlot",tabName = "barplot"                  
    ),tabName = "barplot"),
    
  
    convertMenuItem(menuItem("Target_Table",tabName = "Target_Table",uiOutput("Target_Response_selector")
                             # selectInput("Target_Response_selector",
                             #             "Select the predict response type",
                             #             choices = c("Positive","Negative","All"),
                             #             selected = "All",
                             #             multiple = FALSE
                             # )
    ),tabName = "Target_Table"),
    convertMenuItem(menuItem("Comparison_Table",tabName = "Comparison_Table",uiOutput("Comparison_Response_selector")
                             # selectInput("Comparison_Response_selector",
                             #             "Select the predict response type",
                             #             choices = c("Positive","Negative","All"),
                             #             selected = "All",
                             #             multiple = FALSE
                             # )
    ),tabName = "Comparison_Table"),
    convertMenuItem(menuItem("Cellline_Table",tabName = "Cellline_Table",uiOutput("Cellline_Response_selector")
                             # selectInput("Cellline_Response_selector",
                             #             "Select the predict response type",
                             #             choices = c("Positive","Negative","All"),
                             #             selected = "All",
                             #             multiple = FALSE
                             # )
    ),tabName = "Cellline_Table"),
    convertMenuItem(menuItem("DrugComb_Analysis_Table",tabName = "DrugComb_Analysis_Table",
                             checkboxInput("show_de",
                                           "Show only significantly different features",
                                           FALSE)
                             
                             #icon=icon("alluvialdiagram"),
    ),tabName = "DrugComb_Analysis_Table"),
    convertMenuItem(menuItem("About",tabName = "about"
                             #icon=icon("alluvialdiagram"),
    ),tabName = "about")

  )
)

body <- dashboardBody(
  # tabItems(
  #   tabItem(tabName = "table",
  #           h2("Uploaded File Contents"),
  #           tableOutput("table")
  #   )
  # )
  tabItems(
    tabItem(tabName = "volcanoplot",
            verbatimTextOutput("click_info",
                               placeholder = TRUE),
            
            # output ggplot volcano
            plotOutput("volcano_plot",
                       width = "100%",
                       height = "600px",
                       hover = "volcano_hover",
                       click = "volcano_click",
                       dblclick = "volcano_dbl_click",
                       brush = brushOpts(
                         id = "volcano_brush",
                         resetOnNew = TRUE)),
            
            # Download button for plot
            downloadButton('download_volcano', 'Download volcano plot as PDF'),
            
            br(),
            br(),
            
            # HIGHLIGHTED DRUGS TABLE -----
            dataTableOutput("drug_highlight_tbl")
    ),
    tabItem(tabName = "heatmap",
            # output ggplot heatmap
            plotOutput("heatmap_plot",width = "100%",height = "600px"),
            # Download button for plot
            downloadButton('download_Heatmap', 'Download Heatmap plot as Image'),
            
            br(),
            br(),
            
    ),
    tabItem(tabName = "circleplot",
            # output circle
            plotOutput("circle_plot",width = "100%",height = "800px"),
            # Download button for plot
            downloadButton('download_Circle', 'Download Circle plot as Image'),
            br(),
            br(),
    ),
    tabItem(tabName = "alluvialdiagram",
            # output ggplot Alluvial Diagram
            plotOutput("alluvial_plot",width = "100%",height = "800px"),
            # Download button for plot
            downloadButton('download_alluvial', 'Download alluvial plot as PDF'),
            
            br(),
            br(),
    ),
    tabItem(tabName = "upsetplot",
            # output circle
            plotOutput("upset_plot",
                       width = "100%",
                       height = "800px"),
            
            # Download button for plot
            downloadButton('download_Upset', 'Download Upset plot as Image'),
            
            br(),
            br(),
    ),
    tabItem(tabName = "barplot",
            # output circle
            plotOutput("bar_plot",
                       width = "100%",
                       height = "800px"),
            
            # Download button for plot
            downloadButton('download_Bar', 'Download Upset plot as Image'),
            
            br(),
            br(),
    ),
    tabItem(tabName = "Target_Table",
            uiOutput("Target_Table_Title"),
            # h2("Target Drug Table"),
            dataTableOutput("Target_Table"),
            downloadButton("downloadTarget", "Download Target Data Table"),
    ),
    tabItem(tabName = "Comparison_Table",
            uiOutput("Comparison_Table_Title"),
            dataTableOutput("Comparison_Table"),
            downloadButton("downloadComparison", "Download Comparison Data Table"),
    ),
    tabItem(tabName = "Cellline_Table",
            uiOutput("Cellline_Table_Title"),
            dataTableOutput("Cellline_Table"),
            downloadButton("downloadCellline", "Download Cellline Data Table"),
    ),
    tabItem(tabName = "DrugComb_Analysis_Table",
            h2("Drug-Comb analysis Table"),
            dataTableOutput("DrugComb_Analysis_Table"),
            downloadButton("downloadDrugComb", "Download DrugComb Analysis Table"),
    ),
    tabItem(tabName = "about",
            includeMarkdown("about.md")
    )
  )
)

# Combine the UI components in dashboardPage
ui <- dashboardPage(
  dashboardHeader(
    title = span("Onko_DrugCombScreen"),
    titleWidth =300
  ),
  sidebar,
  body
)

# Define server logic
options(shiny.maxRequestSize = 1024*1024^2)   ### upload file size limitation
server <- function(input, output,session) {
  # This reactive function will capture the data once the files are uploaded
  
  # A helper function to read the data based on file extension
  read_data <- function(inFile) {
    if (is.null(inFile)) {
      return(NULL)
    }
    
    ext <- tolower(tools::file_ext(inFile$name))
    if (ext == "csv") {
      df <- read.csv(inFile$datapath, header = TRUE, sep = ",")
    } else if (ext == "maf" || ext == "gz") {
      df <- MAF_to_SNVTable(inFile$datapath,txDBPath)
      # Placeholder for code to handle MAF or VCF files
      # This is a simplified example; actual reading will depend on the file structure
      # df <- read.table(inFile$datapath, header = TRUE, sep = "\t")
    } else {
      stop("File type not supported.")
    }
    
    return(df)
  }
  
  # Define reactiveValues to store data for each input
  Data_Input <- reactiveValues(    ### reactiveValues() but not the reactive()
    target_primary_data = NULL,
    comparison_primary_data = NULL,
    cellline_data = NULL,
    Target_Table = NULL,
    Comparison_Table = NULL,
    Cellline_Table = NULL,
    DrugComb_Analysis_Table = NULL,
    data_PlotDF = NULL,
    CirclePlot_DFList = NULL
    
    
  )
  selected_PatCancerType <- reactive({
    input$cancer_selector
  })
  Percentage_Cutoff <- reactive({
    input$Percentage_slider
  })
  Testtype <- reactive({
    input$Testtype_selecter
  })
  Response_type_selection <- reactive({
    input$Response_type_selector
  })
  #### Get target_primary_data,comparison_primary_data,cellline_data and their names ------------
  # Observe changes in target_primary_data
  observeEvent(input$target_primary_data, {
    Data_Input$target_primary_data <- read_data(input$target_primary_data)
    Data_Input$Target_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$target_primary_data)
  })
  
  # Observe changes in comparison_primary_data
  observeEvent(input$comparison_primary_data, {
    Data_Input$comparison_primary_data <- read_data(input$comparison_primary_data)
    Data_Input$Comparison_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$comparison_primary_data)
  })
  
  # Observe changes in cellline_data
  observeEvent(input$cellline_data, {
    Data_Input$cellline_data <- read_data(input$cellline_data)
    Data_Input$Cellline_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$cellline_data)
  })
  
  # Define a reactive value to hold the filename
  target_filename <- reactiveVal(NULL)
  comparison_filename <- reactiveVal(NULL)
  cellline_filename <- reactiveVal(NULL)

  # Observe and update target_filename based on file input
  observeEvent(input$target_primary_data, {
    req(input$target_primary_data)  # Ensure the file is uploaded
    string <- tools::file_path_sans_ext(basename(input$target_primary_data$name))
    match <- regexpr("[[:alpha:]]+[[:digit:]]*", string)
    target_filename(regmatches(string, match))
  })
  
  # Observe and update comparison_filename based on file input
  observeEvent(input$comparison_primary_data, {
    req(input$comparison_primary_data)
    string <- tools::file_path_sans_ext(basename(input$comparison_primary_data$name))
    match <- regexpr("[[:alpha:]]+[[:digit:]]*", string)
    comparison_filename(regmatches(string, match))
  })
  
  # Observe and update cellline_filename based on file input
  observeEvent(input$cellline_data, {
    req(input$cellline_data)
    string <- tools::file_path_sans_ext(basename(input$cellline_data$name))
    match <- regexpr("[[:alpha:]]+[[:digit:]]*", string)
    cellline_filename(regmatches(string, match))
  })
  
  ####---------------------
  # observeEvent(input$Response_type_selector, {
  #   if(!is.null(Data_Input$Target_Table) & !is.null(Data_Input$Comparison_Table) & !is.null(Data_Input$Comparison_Table)){
  #     PredictType <- input$Response_type_selector
  #     Data_Input$Target_Table <- Get_DrugPredictsType_DF(Data_Input$Target_Table, PredictType)
  #     Data_Input$Comparison_Table <- Get_DrugPredictsType_DF(Data_Input$Comparison_Table, PredictType)
  #     Data_Input$Cellline_Table <- Get_DrugPredictsType_DF(Data_Input$Cellline_Table, PredictType)}
  # })
  # Observe changes in the cancer_selector dropdown
  observeEvent(c(input$cancer_selector,input$Response_type_selector), {
    
    # Update tables (or perform other operations) when cancer_selector value changes
    if (!is.null(Data_Input$target_primary_data)) {
      Data_Input$Target_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$target_primary_data)
      if(!is.null(Data_Input$Target_Table))
        Data_Input$Target_Table <- Get_DrugPredictsType_DF(Data_Input$Target_Table, Response_type_selection())
    }
    
    if (!is.null(Data_Input$comparison_primary_data)) {
      Data_Input$Comparison_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$comparison_primary_data)
      if(!is.null(Data_Input$Comparison_Table))
        Data_Input$Comparison_Table <- Get_DrugPredictsType_DF(Data_Input$Comparison_Table, Response_type_selection())
    }
    
    if (!is.null(Data_Input$cellline_data)) {
      Data_Input$Cellline_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$cellline_data)
      if(!is.null(Data_Input$Cellline_Table))
        Data_Input$Cellline_Table <- Get_DrugPredictsType_DF(Data_Input$Cellline_Table, Response_type_selection())
    }
  })
  
  
  output$Target_Response_selector <- renderUI({
    selectInput("Target_Response_selector",
                "Select the predict response type",
                choices = c("Positive","Negative","All"),
                selected = Response_type_selection(),
                multiple = FALSE
    )
  })
  output$Comparison_Response_selector <- renderUI({
    selectInput("Comparison_Response_selector",
                "Select the predict response type",
                choices = c("Positive","Negative","All"),
                selected = Response_type_selection(),
                multiple = FALSE
    )
  })
  output$Cellline_Response_selector <- renderUI({
    selectInput("Cellline_Response_selector",
                "Select the predict response type",
                choices = c("Positive","Negative","All"),
                selected = Response_type_selection(),
                multiple = FALSE
    )
  })
  
  observeEvent(input$Target_Response_selector, {if (!is.null(Data_Input$target_primary_data)) {
    Data_Input$Target_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$target_primary_data)
    Data_Input$Target_Table <- Get_DrugPredictsType_DF(Data_Input$Target_Table, input$Target_Response_selector)}
  })
  observeEvent(input$Comparison_Response_selector, {if (!is.null(Data_Input$comparison_primary_data)) {
    Data_Input$Comparison_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$comparison_primary_data)
    Data_Input$Comparison_Table <- Get_DrugPredictsType_DF(Data_Input$Comparison_Table, input$Comparison_Response_selector)
  }
    
  })
  observeEvent(input$Cellline_Response_selector, {if (!is.null(Data_Input$cellline_data)) {
    Data_Input$Cellline_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$cellline_data)
    Data_Input$Cellline_Table <- Get_DrugPredictsType_DF(Data_Input$Cellline_Table, input$Cellline_Response_selector)}
  })
  
  # Render the table for target_primary_data
  output$Target_Table <- renderDataTable({
    if (is.null(Data_Input$Target_Table)) return(NULL)
    Data_Input$Target_Table
    # Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$target_primary_data)
  },options = list(scrollX = TRUE))
  output$Target_Table_Title <- renderUI({
    h2(paste("Target Data -", target_filename()), "Table")
  })
  output$downloadTarget <- downloadHandler(
    filename = function() {
      paste0("primary_data-",target_filename(),"_",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      # data <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$target_primary_data)
      write.csv(Data_Input$Target_Table, file, row.names = FALSE,na="")
    }
  )

  
  output$Comparison_Table <- renderDataTable({
    if (is.null(Data_Input$Comparison_Table)) return(NULL)
    Data_Input$Comparison_Table
    # Additional processing for comparison_primary_data using selected_PatCancerType()
    # Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$comparison_primary_data)
    
  },options = list(scrollX = TRUE))
  output$Comparison_Table_Title <- renderUI({
    h2(paste("Comparison Data -", comparison_filename()), "Table")
  })
  output$downloadComparison <- downloadHandler(
    filename = function() {
      paste0("comparison_data-",comparison_filename(),"_",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(Data_Input$Comparison_Table, file, row.names = FALSE,na="")
    }
  )
  
  output$Cellline_Table <- renderDataTable({
    if (is.null(Data_Input$Cellline_Table)) return(NULL)
    # Additional processing for cellline_data using selected_PatCancerType()
    # Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), Data_Input$cellline_data)
    Data_Input$Cellline_Table
  },options = list(scrollX = TRUE))
  output$Cellline_Table_Title <- renderUI({
    h2(paste("Cellline Data -", cellline_filename()), "Table")
  })
  output$downloadCellline <- downloadHandler(
    filename = function() {
      paste0("cellline_data-",cellline_filename(),"_",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(Data_Input$Cellline_Table, file, row.names = FALSE,na="")
    }
  )
  
  observeEvent(input$submit, {
    # Use withProgress to show a progress indicator
    withProgress(message = 'Processing... Please wait', value = 0, {
      # Optionally, set an initial progress amount
      incProgress(0.1)  # Shows initial progress
      
      # Check if any of the tables are NULL
      if (is.null(Data_Input$Target_Table) || is.null(Data_Input$Comparison_Table)) {
        Data_Input$DrugComb_Analysis_Table <- NULL  # Ensures the table output is cleared
        incProgress(0.9, detail = "Clearing results...")  # Completes the progress
      } else if(is.null(Data_Input$Cellline_Table)) {
        Data_Input$DrugComb_Analysis_Table <- DrugComb_analysis(Data_Input$Target_Table, Data_Input$Comparison_Table, Freq_cutoff = Percentage_Cutoff(), testType = Testtype())
        incProgress(0.5, detail = "Processing without cell line data...")  # Adjusts the progress to complete
        Data_Input$CirclePlot_DFList <- CirclePlot_DF_Create(Data_Input$DrugComb_Analysis_Table)
        incProgress(0.9, detail = "Processing CirclePlot_DFlist ...")  # Adjusts the progress to complete
      } else {
        Data_Input$DrugComb_Analysis_Table <- DrugComb_analysis(Data_Input$Target_Table, Data_Input$Comparison_Table, Data_Input$Cellline_Table, Freq_cutoff = Percentage_Cutoff(), testType = Testtype())
        incProgress(0.4, detail = "Processing with cell line data...")  # Adjusts the progress to complete
        Data_Input$data_PlotDF <- PlotDF_Create(Data_Input$DrugComb_Analysis_Table)
        incProgress(0.7, detail = "Processing PlotDF ...")  # Adjusts the progress to complete
        Data_Input$CirclePlot_DFList <- CirclePlot_DF_Create(Data_Input$DrugComb_Analysis_Table)
        incProgress(0.9, detail = "Processing CirclePlot_DFList...")  # Adjusts the progress to complete
      }
    })
  })

  observeEvent(input$test_data, {
   
    Data_Input$target_primary_data <- test_Target_data
    Data_Input$comparison_primary_data <- test_Comparison_data
    Data_Input$cellline_data <- test_Cellline_data
    
    Data_Input$Target_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), test_Target_data)
    Data_Input$Comparison_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), test_Comparison_data)
    Data_Input$Cellline_Table <- Get_MTBreporter_DF(DrugClassified_DB, selected_PatCancerType(), test_Cellline_data)
    # print(Data_Input$Target_Table)

    
    withProgress(message = 'Processing... Please wait', value = 0, {
      # Optionally, set an initial progress amount
    incProgress(0.1)  # Shows initial progress

    Data_Input$DrugComb_Analysis_Table <- DrugComb_analysis(Data_Input$Target_Table, Data_Input$Comparison_Table, Data_Input$Cellline_Table, Freq_cutoff = Percentage_Cutoff(), testType = Testtype())
    incProgress(0.4, detail = "Processing with cell line data...")  # Adjusts the progress to complete
    Data_Input$data_PlotDF <- PlotDF_Create(Data_Input$DrugComb_Analysis_Table)
    incProgress(0.7, detail = "Processing PlotDF ...")  # Adjusts the progress to complete
    Data_Input$CirclePlot_DFList <- CirclePlot_DF_Create(Data_Input$DrugComb_Analysis_Table)
    incProgress(0.9, detail = "Processing CirclePlot_DFList...")  # Adjusts the progress to complete

      
    })

    target_filename("test_Target")
    comparison_filename("test_Comparison")
    cellline_filename("test_Cellline")
    

  })
  
  
  output$DrugComb_Analysis_Table <- renderDataTable({
    Data_Input$DrugComb_Analysis_Table
  },options = list(scrollX = TRUE,pageLength = 3,autoWidth = TRUE))
  
  is_de <- reactive({
    data <-  Data_Input$data_PlotDF
    abs(data[[c("log2oddsRatio")]]) >= input$logfc_threshold & data[[input$pvalue_col]] <= input$pvalue_threshold
  })
  
  de_DrugComb_Analysis_Table <- reactive({
    data <- Data_Input$DrugComb_Analysis_Table 
    if (input$show_de) {
      filter(data, is_de())
    } else {
      data
    }
  })
  # Render the DrugComb_Analysis_Table with significant difference
  output$DrugComb_Analysis_Table <- renderDataTable({
    de_DrugComb_Analysis_Table()
  },options = list(scrollX = TRUE,pageLength = 3,autoWidth = TRUE))
  
  
  output$downloadDrugComb <- downloadHandler(
    filename = function() {
      paste0("DrugComb_analysis","_",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(de_DrugComb_Analysis_Table(), file, row.names = FALSE,na="")
    }
  )



  #' Region Start: Volcano Plot and Interaction Handlers
  # ------------------------------------------------------------------
  # X AND Y AXES LABELER -----
  
  # capture pvalue column selected and default value with it
  reactive_pvalue_value <- reactive({
    paste0("-log10(", input$pvalue_col, ")")
  })
  
  # enter custom x (logfc) axis label
  output$x_axis_labeler <- renderUI({
    textInput("x_axis_lab",
              "Specify X axis label",
              value = "log2oddsRatio",
              placeholder = "ex: log2oddsRatio")
  })
  
  # enter custom x (logfc) axis label
  output$y_axis_labeler <- renderUI({
    textInput("y_axis_lab",
              "Specify Y axis label",
              value = reactive_pvalue_value(),
              placeholder = "ex: log_pval")
  })
  
  # HIGHLIGHTED DRUG TABLE -----
  
  # select drugs to highlight
  output$drugcombs_selector <- renderUI({
    data <- Data_Input$data_PlotDF
    selectInput("selected_drugs",
                "Select feature(s) to label",
                sort(data[["Drug_comb"]]),
                multiple = TRUE,
                selectize= TRUE)
  })
  
  # initialize drug_list$clicked_drug_list as NULL
  # This will reactively update
  drug_list <- reactiveValues(clicked_drugs_list = NULL)
  
  # store clicked drug info
  clicked_drug <- reactive({
    if (is.null(data_w_log_pval())) {
      return(NULL) # Early exit if data_w_log_pval is NULL
    } else {
      nearPoints(data_w_log_pval(),
                 input$volcano_click,
                 xvar = "log2oddsRatio",
                 yvar = "log_pval",
                 maxpoints = 1) %>%
        dplyr::select("Drug_comb")
    }
  })
  
  # when a point is clicked on the volcano plot
  # add drug to clicked drug list
  # if the point has been clicked twice, remove from list
  observeEvent(input$volcano_click, {
    # Check if data for the plot is available
    if (!is.null(Data_Input$data_PlotDF)) {
      # Proceed with the existing logic only if data is available
      if (is.null(input$selected_drugs)) {
        drug_list$clicked_drugs_list <- NULL
      }
      if (is.null(drug_list$clicked_drugs_list)) {
        drug_list$clicked_drugs_list <- clicked_drug()
      } else {
        drug_present <- clicked_drug() %in% input$selected_drugs
        if (drug_present) {
          present_idx <- !grepl(clicked_drug(), input$selected_drugs)
          drug_list$clicked_drugs_list <- input$selected_drugs[present_idx]
        } else {
          drug_list$clicked_drugs_list <- c(clicked_drug(), input$selected_drugs)
        }
      }
    }
  })
  
  
  observe({
    if(is.null(input$shinymanager_where) || (!is.null(input$shinymanager_where) && input$shinymanager_where %in% "application")){
      data <- Data_Input$data_PlotDF
      updateSelectInput(session, 
                        "selected_drugs",
                        label = "Select feature(s) to label",
                        choices = sort(data[["Drug_comb"]]),
                        selected = drug_list$clicked_drugs_list)
    }
  })
  
  # reactive function that subsets data by highlighted_drug vector
  highlight_drugs_data <- reactive({
    data <- Data_Input$data_PlotDF
    if (length(input$selected_drugs) > 0) {
      highlight_drugs_data <- data[data[["Drug_comb"]] %in% input$selected_drugs, c("Drug_comb", "log2oddsRatio", input$pvalue_col)]
    } else {
      highlight_drugs_data <- data.frame(NA, NA, NA)
      names(highlight_drugs_data) <- c("Drug_comb", "log2oddsRatio", input$pvalue_col)
    }
  })
  
  # render a data table of highlighted drug info
  output$drug_highlight_tbl <- renderDataTable({
    highlight_drugs_data()
  })
  
  
  
  #this is the value that will be input into volcanoPlot()
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # when there is a double click on the plot
  # if brush is null, nothing happens,
  # if brush is not null, assign values to ranges
  observeEvent(input$volcano_dbl_click, {
    # Check if data for the plot is available
    if (!is.null(Data_Input$data_PlotDF)) {
      brush <- input$volcano_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    }
  })
  
  # volcano plot in reactive function (is this necessary?? can't be sure.)
  reactive_volcano <- reactive({
    # if(is.null(Data_Input$target_primary_data)) return(NULL)
    plotVolcano(data = Data_Input$data_PlotDF, 
                pvalue_col = input$pvalue_col, 
                pvalue_thresh = input$pvalue_threshold, 
                logfc_thresh = input$logfc_threshold,
                de_vec = is_de(),
                show_logfc_thresh = input$show_logfc_threshold,
                show_pvalue_thresh= input$show_pvalue_threshold,
                highlight_drugcombs =  input$selected_drugs,
                x_label = input$x_axis_lab,
                y_label = input$y_axis_lab,
                xlim = ranges$x,
                ylim = ranges$y,
                specific_category = paste0(target_filename(),"_vs_",comparison_filename())
                
    )
  })
  
  # output volcano plot
  output$volcano_plot <- renderPlot({
    if (is.null(Data_Input$data_PlotDF)) return(NULL)
    reactive_volcano()
  })
  # DISPLAY DRUG INFO ON HOVER OVER -----
  
  # Create -log10 pvalue column
  data_w_log_pval <- reactive({
    # Check if the data_PlotDF is not NULL
    if (!is.null(Data_Input$data_PlotDF)) {
      # Proceed with data transformation since data is available
      data <- Data_Input$data_PlotDF
      reduced_data <- data %>%
        mutate(log_pval = -log10(as.numeric(data[[input$pvalue_col]])))
      return(reduced_data)
    } else {
      # Return NULL or an empty data frame if data_PlotDF is NULL
      return(NULL)
    }
  })
  
  # Collect nearpoint info and reduce to only Drug_comb, log2oddsRatio and pvalue_col
  point_info <- reactive({
    # Ensure data_w_log_pval() is not NULL before proceeding
    if (!is.null(data_w_log_pval())) {
      nearpoint_out <- nearPoints(data_w_log_pval(), input$volcano_hover, xvar = "log2oddsRatio", yvar = "log_pval", maxpoints = 1)
      # Make sure to use correct column name without .data pronoun and directly select columns
      nearpoint_out %>%
        dplyr::select("Drug_comb", "log2oddsRatio", input$pvalue_col)
    } else {
      # Return NULL or an informative placeholder if the data is not available
      return(data.frame(Drug_comb = NA, log2oddsRatio = NA, Pvalue = NA))
    }
  })
  # render printed text
  output$click_info <- renderPrint({
    # Fetch the information from point_info() reactive expression
    info <- point_info()
    
    # Check if info has any rows to print
    if (!is.null(info) && nrow(info) > 0) {
      print(info)
    } else {
      # Optionally, print a message indicating that no data is available
      # This is more user-friendly than showing a blank or an error
      print("Hover over the plot to see candidate drug combination details.")
    }
  })

  
  # DOWNLOAD HANDLER ----- VOLCANO
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0(target_filename(),"_vs_",comparison_filename(),"_volcano_plot_", Sys.Date(), ".pdf")
    },
    
    content = function(file) {
      ggsave(file, reactive_volcano(), device = "pdf", width = input$img_w/64, height = input$img_h/64, units = "in")
    })
  # ------------------------------------------------------------------
  #' Region End
  
  ### Heatmap plot in reactive function ----------------
  reactive_heatmap <- reactive({
    data <- Data_Input$DrugComb_Analysis_Table
    plotHeatmap(data = data,
                Cellline_cutoff = input$cellline_threshold,
                col_high = input$col_high,
                col_mid = input$col_mid,
                col_low = input$col_low
    )
  })

  # output heatmap plot
  output$heatmap_plot <- renderPlot({
    if (is.null(Data_Input$DrugComb_Analysis_Table)) return(NULL)
    reactive_heatmap()
  })
  # DOWNLOAD HANDLER ----- Heatmap

  output$download_Heatmap <- downloadHandler(
    filename = function() {
      paste0(target_filename(),"_Heatmap_", Sys.Date(), ".png")
    },

    # content = function(file) {
    #   ggsave(file, reactive_heatmap(), device = "pdf", width = 10, height = 5, units = "in")
    #   dev.off()
    # }
    content = function(file) {
      png(filename = file,width = input$img_w, height = input$img_h)

      draw(reactive_heatmap())
      dev.off()
    }
  )
  ### Circle Plot in reactive function ----------------
  reactive_circle <- reactive({
    plotCircle(CirclePlot_DF = Data_Input$CirclePlot_DFList[["CirclePlot_DF"]],
               DrugNameDF = Data_Input$CirclePlot_DFList[["DrugNameDF"]],
               title_name = target_filename()
    )
  })
  # output circle plot
  output$circle_plot <- renderPlot({
    if (is.null(Data_Input$CirclePlot_DFList)) return(NULL)
    reactive_circle()
  })
  # DOWNLOAD HANDLER ----- Circle
  output$download_Circle <- downloadHandler(
    filename = function() {
      paste0("Circle_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(filename = file,width = input$img_w, height = input$img_h)
      # reactive_circle()
      plotCircle(CirclePlot_DF = Data_Input$CirclePlot_DFList[["CirclePlot_DF"]],
                 DrugNameDF = Data_Input$CirclePlot_DFList[["DrugNameDF"]],
                 title_name = target_filename()
      )
      dev.off()
    }
  )
  
  
  ### Alluvial diagram  in reactive function ----------------
  # select drugs to highlight
  output$drug_selector <- renderUI({
    
    selectInput("selected_alluvialdrug",
                "Select a drug for alluvial diagram",
                sort(Data_Input$Target_Table[["Classified_Drug_Name"]]),
                # options = list(`server` = TRUE),
                multiple = FALSE,
                #selected = "PARP Inhibitor",
                selectize= TRUE
    )
  })
    reactive_alluvial <- reactive({
      plotAlluvial(drugDF = Data_Input$Target_Table,
                   classified_drugName = input$selected_alluvialdrug,
                   Alluvial_title = paste(input$selected_alluvialdrug,"in",input$subtype_drug,target_filename())
                   # Alluvial_title = target_filename()
      )
    })
    
    output$alluvial_plot <- renderPlot({
      if (is.null(Data_Input$Target_Table)) return(NULL)
      reactive_alluvial()
    })
    
    output$download_alluvial <- downloadHandler(
      filename = function() {
        paste0(target_filename(),"_alluvial_plot_", Sys.Date(),".png")
      },
      content = function(file) {
        png(filename = file,width = input$img_w, height = input$img_h)
        
        print(reactive_alluvial())
        dev.off()
      }
    )
  
  ### Upset Plot in reactive function ----------------
  reactive_upset <- reactive({
    
    plotUpset(drugDF = Data_Input$Target_Table,
              title_name = target_filename()
    )
  })
  
  # output bar plot
  output$upset_plot <- renderPlot({
    if (is.null(Data_Input$Target_Table)) return(NULL)
    reactive_upset()
  })
  
  output$download_Upset <- downloadHandler(
    filename = function() {
      paste0(target_filename(), "_Upset_plot_",Sys.Date(), ".png")
    },

    content = function(file) {
      png(filename = file,width = input$img_w, height = input$img_h)
      
      print(reactive_upset())
      dev.off()
    }
  )
  ### Bar plot in reactive function ---------------- 
  reactive_bar <- reactive({
    
    plotfrequencyBarplot(drugDF = Data_Input$Target_Table,
                         Barplot_title = target_filename()
    )
  })
  
  # output circle plot
  output$bar_plot <- renderPlot({
    if (is.null(Data_Input$Target_Table)) return(NULL)
    reactive_bar()
  })
  
  output$download_Bar <- downloadHandler(
    filename = function() {
      paste0(target_filename(), "_SingleDurgFreq_Barplot_",Sys.Date(), ".png")
    },
    content = function(file) {
      png(filename = file,width = input$img_w, height = input$img_h)
      
      print(reactive_bar())
      dev.off()
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)