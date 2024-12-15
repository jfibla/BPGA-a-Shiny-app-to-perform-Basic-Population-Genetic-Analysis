# BPGA-a-Shiny-app-to-perform-Basic-Population-Genetic-Analysis
#~/bpga_ngroc/app.R


library(shiny)
library(qqman)
library(readr)
library(data.table)
library(ggplot2)
library(shinyjs)
library(dplyr)
library(tidyr)
library(rmarkdown)
library(knitr)
library(ggforce)
library(mapplots)
library(maps)
library(ggh4x)
library(conflicted)

conflicts_prefer(
  maps::map(),
  dplyr::count())

# Increase max upload size to 500MB
options(shiny.maxRequestSize = 500 * 1024^2)


ui <- fluidPage(
  #suppress warnning message at console                                         <<<<<<<<<<<<<<<<<<<<<<<
  tags$head(
    tags$style(
      HTML("
        body {
          background-image: url('fons.jpg');
          background-size: cover;
          background-repeat: no-repeat;
          background-attachment: fixed;
        }
      ")
    )
  ),
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  useShinyjs(),  # Initialize shinyjs
  titlePanel("Basic Population Genetic Analysis - (BPGA)"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel("output.panel_visible === false",  #>>>>> cond1
                       ########################## Text explanation project
                       absolutePanel(id = "Text explanation", class = "panel panel-default", fixed = F,
                                     draggable = T, top = 100, left = 40, right =40, bottom = "auto",
                                     width = "auto", height ="auto",
                                     style = "background-color: #E8E8E8;
                         opacity: 0.85;
                         z-index: 1000;
                         padding: 20px 20px 20px 20px;
                         margin: auto;
                         border-radius: 4px,
                         box-shadow: 0 0 10px rgba(0,0,0,0.2);
                         padding-bottom: 2mm;
                         padding-top: 1mm;",       
                                     h3("Project description"),
                                     fluidRow(column(11, p(style="text-align: justify;","This Shiny app provides an interactive platform for visualizing population genetic structure. Users can upload their own datasets or 
                                                           use example files to explore population relationships, genetic diversity, and ancestral components. The app dynamically reads PLINK binary files and integrates them with worldwide 
                                                           reference populations from the 1000 Genomes Project (1000G) and the Human Genome Diversity Project (HGDP), enabling basic population analyses such as Principal Component Analysis (PCA), 
                                                           ADMIXTURE, and FST analyses.")),
                                              column(11, p(style="text-align: justify;","Designed for both researchers and students in population genetics, the app offers a user-friendly interface for exploring genetic 
                                                           structure through multiple methods. It facilitates interactive data exploration and produces publication-ready plots.")),
                                              column(11,strong("Click on "), 
                                                     actionButton("start_analysis", "Start analysis"),
                                                     strong("button and enjoi"))
                                                     #actionButton("reload_button", "Reload session"))
))
), #>>>>> /cond1

#tags$a(href="https://github.com/jfibla/BPGA-a-Shiny-app-to-perform-Basic-Population-Genetic-Analysis.git", 
#       h5( "Code available")),
                                     
      # to generate user folder
   #   actionButton("start_analysis", "Start analysis"),
      actionButton("reload_button", "Reload session"),
    #  h4("Main folders:"),
      verbatimTextOutput("folder_path"),
conditionalPanel(condition = "output.panel_visible === true",     #>>>>> cond2
    selectInput("files", "Select input file:",
                  c("Example single pop" = "exp1","Example multi pop" = "exp2",
                    "User file" = "usr")),
conditionalPanel(condition = "input.files == 'usr'",             #>>>>> cond3
      fileInput("upload_file", "Load user population data (compressed .zip or .gz)", 
                accept = c(".zip", ".gz")),
      textInput("popID", "Assign population identifier (use a three-letter code like 'USR')", "USR"),
      textInput("lat", "Assign latitude"),
      textInput("lon", "Assign longitude")
      ),                                                         #>>>>> /cond3
      h3("Process input File"),
tags$hr(style="border-color: grey;"),

      radioButtons("selector", "Select option:", 
              choices = list("Only loaded population" = "usr","Merge with Woldwide populations" = "wwp")),
conditionalPanel(condition = "input.selector == 'usr'",    #>>>>> cond3.1  
      actionButton("process_file", "Decompress file"),
),                                                         #>>>>> /cond3.1
conditionalPanel(condition = "input.selector == 'wwp'",           #>>>>> cond4
      actionButton("run_merge", "Decompress and merge files")
      ),                                                          #>>>>> /cond4
      h3("Perform population analysis"),
      uiOutput("pop1_select"),  # Dynamic SelectInput for POP1
      shinyjs::hidden(actionButton("select_all_pop1", "Select all")),
      tags$hr(style="border-color: grey;"),
 
      h4("PCA analysis"),
      sliderInput("sliderM1", "Set MAF for prunning:", min = 0.05, max = 0.50, value = 0.05, step = 0.05),
      actionButton("run_pca", "Run PCA analysis"),
verbatimTextOutput("log_output"),
tags$hr(style="border-color: grey;"),

#conditionalPanel(condition = "output.run_pca === true",           #>>>>> cond5
      h4("Admixture analysis"),
      sliderInput("sliderM2", "Set MAF:", min = 0.05, max = 0.50, value = 0.05, step = 0.05),
verbatimTextOutput("maf_output"),  
sliderInput("sliderK", "Set K values:", min = 1, max = 10, value = 1),

      actionButton("run_admix", "Run ADMIXTURE"),
verbatimTextOutput("cvErrors"),
tags$hr(style="border-color: grey;"),
      h4("FST analysis: compare two populations"),
 #     h4("FST analysis"),
      uiOutput("pop2_checkbox"),  # Dynamic Checkbox for POP2
      shinyjs::hidden(actionButton("select_all_pop2", "Remove selected")),
      actionButton("run_fst", "Run FST analysis"),
verbatimTextOutput("max_fst_output"),  # Output for displaying the max FST rowtags$hr(style="border-color: grey;"),
h4("Decoration plots"),
radioButtons("popcolor", "Colored:", 
             choices = list("Superpopulation" = "POP1", "Subpopulation" = "POP2")),
sliderInput("sliderR", "Set Pie radius:", min = 0, max = 100, value = 4, step = 1),
sliderInput("sliderE", "Expand map:", min = 0, max = 2, value = 0.2, step = 0.05),
#       ),                                                               #>>>>> /cond5
uiOutput("sorter"),
#selectInput("sorter", "Sort admixture boxplot:",
#                 c("By first component" = "V1",
#                   "By second component" = "V2")),
tags$hr(style="border-color: grey;"),
      h4("Report"),
      downloadButton("downloadReport", "Download Report"),
    )                                                             #>>>>> /cond2
    ),

###################################################################### 
mainPanel(
      
      textOutput("status"),
      uiOutput("progress_bar_ui"),
      plotOutput("PCA_plot"),
      plotOutput("ADMIXTURE_plot"),
      plotOutput("box_plot"), 
      plotOutput("map_plot", height = "700px"),
      plotOutput("manhattan_plot")
      
    )
  )
)
################################################################################  
########################### SERVER #############################################  
################################################################################  

server <- function(input, output, session) {
  # value of slider M as min slider M1
  observeEvent(input$sliderM1, {
    updateSliderInput(session, "sliderM2", min = input$sliderM1)
  })
  
  # Reactive value to track PLINK completion
 
  
  # Reactive value to track button state
  panel_visible <- reactiveVal(FALSE)
  
  observeEvent(input$start_analysis, {
    # Toggle the value of the panel visibility
    panel_visible(!panel_visible())
  })
  
  # Send the value of panel_visible to the client
  output$panel_visible <- reactive({
    panel_visible()
  })
  outputOptions(output, "panel_visible", suspendWhenHidden = FALSE)
  
  
  observeEvent(input$show_content, {
    req(input$folder) # Ensure the folder path is provided
    
    # Render the folder content
    output$folder_content <- renderText({
      folder_path <- input$folder
      
      # Check if the folder exists
      if (!dir.exists(folder_path)) {
        return(paste("Folder does not exist:", folder_path))
      }
      
      # Get the list of files and folders inside the directory
      folder_content <- list.files(folder_path, full.names = TRUE)
      
      if (length(folder_content) == 0) {
        return(paste("The folder is empty:", folder_path))
      }
      
      # Format and display the content
      paste("Content of", folder_path, ":\n", paste(folder_content, collapse = "\n"))
    })
  })
  processed_file_path <- reactiveVal(NULL)  # Store processed file path
  # Disable the PCA, ADMIXTURE, and FST buttons by default
  shinyjs::disable("run_pca")
  shinyjs::disable("run_admix")
  shinyjs::disable("run_fst")
  
  ################################################################################  
 # create a random folder and subfolders
    random_folder <- reactiveVal()
    
    # Crear una carpeta aleatòria quan l'usuari comenci l'anàlisi
    observeEvent(input$start_analysis, {
      # Generar un nom únic per la carpeta principal
      folder_path <- tempfile(pattern = "user_output_", tmpdir = tempdir())
      dir.create(folder_path)  # Crear la carpeta principal
      
      # Crear subcarpetes dins la carpeta principal
      subfolders <- c("decompressed_files","pca", "merged", "temporal", "admx", "fst","plots")
      lapply(subfolders, function(subfolder) {
        dir.create(file.path(folder_path, subfolder))  # Crear cada subcarpeta
      })
      
      # Assignar la carpeta principal al reactiu
      random_folder(list(
        main_folder = folder_path,
        subfolders = subfolders
      ))
      
      # Mostrar el camí de la carpeta principal i les subcarpetes
      output$folder_path <- renderText({
        paste("Temporal user folder:", folder_path)
      })
      
      output$subfolder_paths <- renderText({
        paste(
          "Subfolders created:\n",
          paste(file.path(folder_path, subfolders), collapse = "\n")
        )
      })
  #  }) 
    

      ################################################################################  
      # Function to decompress and label the uploaded file
      observeEvent(input$process_file, {
        # Ensure only one progress object is created at the start
        progress <- shiny::Progress$new()
        progress$set(message = "Processing file...", value = 0)
        on.exit({
          progress$close()  # Ensure progress bar is closed at the end of processing
        }, add = TRUE)
        
        if (input$selector == "usr") {
          shinyjs::enable("run_pca")
        }
        
        if (input$files == "exp1" || input$files == "exp2"){                   ############### if the example file is selected
          # go to folder with the example file
          if (input$files == "exp1"){
            ex_file_path <- "example/example_single.zip"  # load example_single files from "example" folder
            file_ext <- tools::file_ext(ex_file_path) # extract file extension
            base_name <- "example" # extract base name
          } else {
            ex_file_path <- "example/example_multi.zip"  # load example_multi files from "example" folder
            file_ext <- tools::file_ext(ex_file_path)
            base_name <- "example"
          }
          if(input$selector == "wwp") {                 ####### if example file is merged with WWW
            unzip(ex_file_path, exdir = file.path(folder_path, "decompressed_files")) #>>>>>> send to "decompressed_files" folder
            decompressed_file_path <- file.path(folder_path, "decompressed_files", base_name)
            # Store the processed file path
            processed_file_path(decompressed_file_path)
            
            # rename files as INPUT mantenint l'extensió original
            file_extension <- tools::file_ext(decompressed_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "decompressed_files", paste0("INPUT", file_extension))
            
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(decompressed_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(decompressed_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(decompressed_file_path, ".bim"), paste0(new_name, ".bim")) 
            # Retorna el nou camí del fitxer processat
            processed_file_path <- new_name
            
            # recover codes from first column
            df2_reactive <- reactive({
              fam <- fread(paste0(file.path(folder_path),"/decompressed_files/INPUT.fam"), header = F)
              fam.1 <- fam[, 1:2]
              colnames(fam.1) = c("IDI", "IDII")
              fam.2 <- fam.1[, 1]
              colnames(fam.2) <- "V1"
              fam2 <- fam.2 %>%
                separate(V1, into = c("POP1", "POP2"), sep = "_")
              df2 <- as.data.frame(cbind(fam2,fam.2))
              return(df2)
            })
            
            # UI dinàmic per a la selecció de POP1
            output$pop1_select <- renderUI({
              df2 <- df2_reactive()
              selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
            })
            
            # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
            output$pop2_checkbox <- renderUI({
              req(input$pop1)  # Assegura que POP1 està seleccionat
              df2 <- df2_reactive()
              pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
              selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
            })
            
            # Update status to indicate file processed successfully
            output$status <- renderText("File processed successfully.")
            
          } else {                                  ####### if example file is processed alone
            
            unzip(ex_file_path, exdir = file.path(folder_path, "merged"))  #>>>>>>> send to "temporal" folder
            temporal_file_path <- file.path(folder_path, "merged", base_name)
            processed_file_path(temporal_file_path)
            
            # rename file as "MERGED" mantenint l'extensió original
            file_extension <- tools::file_ext(temporal_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "merged", paste0("MERGED", file_extension))
            
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(temporal_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(temporal_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(temporal_file_path, ".bim"), paste0(new_name, ".bim"))
            
            # Retorna el nou camí del fitxer processat
            processed_file_path <- new_name
            
            # recover codes from first column
            df2_reactive <- reactive({
              fam <- fread(paste0(file.path(folder_path),"/merged/MERGED.fam"), header = F)
              fam.1 <- fam[, 1:2]
              colnames(fam.1) = c("IDI", "IDII")
              fam.2 <- fam.1[, 1]
              colnames(fam.2) <- "V1"
              fam2 <- fam.2 %>%
                separate(V1, into = c("POP1", "POP2"), sep = "_")
              df2 <- as.data.frame(cbind(fam2,fam.2))
              return(df2)
            })
            
            # UI dinàmic per a la selecció de POP1
            output$pop1_select <- renderUI({
              df2 <- df2_reactive()
              selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
            })
            
            # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
            output$pop2_checkbox <- renderUI({
              req(input$pop1)  # Assegura que POP1 està seleccionat
              df2 <- df2_reactive()
              pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
              selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
            })
          }
          # File processed successfully
          output$status <- renderText("File processed successfully.")
          #    progress$inc(1)  # Increment progress
          
        } else {                                 ################ selection of file loaded by the user
          req(input$upload_file)
          
          usr_file_path <- input$upload_file$datapath
          file_ext <- tools::file_ext(usr_file_path)
          
          base_name <- NULL
          
          if(input$selector == "wwp") {             ########### if loaded file must be merged with WWW 
            
            if (file_ext == "zip") {
              unzip(usr_file_path, exdir = file.path(folder_path, "decompressed_files")) #>>>>>> send to "decompressed_files" folder
              base_name <- tools::file_path_sans_ext(input$upload_file$name)
            } else if (file_ext == "gz") {
              base_name <- tools::file_path_sans_ext(basename(usr_file_path))
              system(paste("gunzip", usr_file_path))
            } else {
              stop("Unsupported file type.")
            }
            decompressed_file_path <- paste0(file.path(folder_path, "decompressed_files", base_name))
            processed_file_path(decompressed_file_path)
            print(base_name)
            # to assign new code to .fam file
            user.fam <- fread(paste0(decompressed_file_path,".fam"), header = F)
            colnames(user.fam) <- c("V1","V2","V3","V4","V5","V6")
            user.id <- paste0(input$popID,"_",input$popID)
            user.fam$V1 <- user.id
            write.table(user.fam, paste0(decompressed_file_path,".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
            
            # Renomena amb el nom "INPUT" mantenint l'extensió original
            file_extension <- tools::file_ext(decompressed_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "decompressed_files", paste0("INPUT", file_extension))
            print(new_name)
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(decompressed_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(decompressed_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(decompressed_file_path, ".bim"), paste0(new_name, ".bim"))
            
          } else {                                ############  if loaded file is processed alone
            if (file_ext == "zip") {
              unzip(usr_file_path, exdir = file.path(folder_path, "merged"))  #>>>>>> send to "merged" folder
              base_name <- tools::file_path_sans_ext(input$upload_file$name)
            } else if (file_ext == "gz") {
              base_name <- tools::file_path_sans_ext(basename(usr_file_path))
              system(paste("gunzip", usr_file_path))
            } else {
              stop("Unsupported file type.")
            }
            merged_file_path <- paste0(file.path(folder_path, "merged", base_name))
            processed_file_path(merged_file_path)
            
            # to assign new code to .fam file
            user.fam <- fread(paste0(merged_file_path,".fam"), header = F)
            colnames(user.fam) <- c("V1","V2","V3","V4","V5","V6")
            user.id <- paste0(input$popID,"_",input$popID)
            user.fam$V1 <- user.id
            write.table(user.fam, paste0(merged_file_path,".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
            
            # Defineix el nou camí amb el nom "MERGED" mantenint l'extensió original
            file_extension <- tools::file_ext(merged_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "merged", paste0("MERGED", file_extension))
            
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(merged_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(merged_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(merged_file_path, ".bim"), paste0(new_name, ".bim"))
            
            # recover codes from first column
            df2_reactive <- reactive({
              fam <- fread(paste0(file.path(folder_path, "merged/MERGED.fam")), header = F)
              fam.1 <- fam[, 1:2]
              colnames(fam.1) = c("IDI", "IDII")
              fam.2 <- fam.1[, 1]
              colnames(fam.2) <- "V1"
              fam2 <- fam.2 %>%
                separate(V1, into = c("POP1", "POP2"), sep = "_")
              df2 <- as.data.frame(cbind(fam2,fam.2))
              return(df2)
            })
            
            # UI dinàmic per a la selecció de POP1
            output$pop1_select <- renderUI({
              df2 <- df2_reactive()
              selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
            })
            
            # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
            output$pop2_checkbox <- renderUI({
              req(input$pop1)  # Assegura que POP1 està seleccionat
              df2 <- df2_reactive()
              pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
              selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
            })
          }
          # Final step - Update status and progress bar
          output$status <- renderText("File processed successfully.")
          progress$inc(1)  # Ensure progress is fully completed 
        }
      })
#####################################################################################  
      # Function to merge files
      observeEvent(input$run_merge, {
        # Ensure only one progress object is created at the start
        progress <- shiny::Progress$new()
        progress$set(message = "Processing file...", value = 0)
        on.exit({
          progress$close()  # Ensure progress bar is closed at the end of processing
        }, add = TRUE)
        
        if (input$selector == "usr") {
          shinyjs::enable("run_pca")
        }
        
        if (input$files == "exp1" || input$files == "exp2"){                   ############### if the example file is selected
          # go to folder with the example file
          if (input$files == "exp1"){
            ex_file_path <- "example/example_single.zip"  # load example_single files from "example" folder
            file_ext <- tools::file_ext(ex_file_path) # extract file extension
            base_name <- "example" # extract base name
          } else {
            ex_file_path <- "example/example_multi.zip"  # load example_multi files from "example" folder
            file_ext <- tools::file_ext(ex_file_path)
            base_name <- "example"
          }
          if(input$selector == "wwp") {                 ####### if example file is merged with WWW
            unzip(ex_file_path, exdir = file.path(folder_path, "decompressed_files")) #>>>>>> send to "decompressed_files" folder
            decompressed_file_path <- file.path(folder_path, "decompressed_files", base_name)
            # Store the processed file path
            processed_file_path(decompressed_file_path)
            
            # rename files as INPUT mantenint l'extensió original
            file_extension <- tools::file_ext(decompressed_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "decompressed_files", paste0("INPUT", file_extension))
            
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(decompressed_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(decompressed_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(decompressed_file_path, ".bim"), paste0(new_name, ".bim")) 
            # Retorna el nou camí del fitxer processat
            processed_file_path <- new_name
            
            # recover codes from first column
            df2_reactive <- reactive({
              fam <- fread(paste0(file.path(folder_path),"/decompressed_files/INPUT.fam"), header = F)
              fam.1 <- fam[, 1:2]
              colnames(fam.1) = c("IDI", "IDII")
              fam.2 <- fam.1[, 1]
              colnames(fam.2) <- "V1"
              fam2 <- fam.2 %>%
                separate(V1, into = c("POP1", "POP2"), sep = "_")
              df2 <- as.data.frame(cbind(fam2,fam.2))
              return(df2)
            })
            
            # UI dinàmic per a la selecció de POP1
            output$pop1_select <- renderUI({
              df2 <- df2_reactive()
              selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
            })
            
            # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
            output$pop2_checkbox <- renderUI({
              req(input$pop1)  # Assegura que POP1 està seleccionat
              df2 <- df2_reactive()
              pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
              selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
            })
            
            # Update status to indicate file processed successfully
            output$status <- renderText("File processed successfully.")
            
          } else {                                  ####### if example file is processed alone
            
            unzip(ex_file_path, exdir = file.path(folder_path, "merged"))  #>>>>>>> send to "temporal" folder
            temporal_file_path <- file.path(folder_path, "merged", base_name)
            processed_file_path(temporal_file_path)
            
            # rename file as "MERGED" mantenint l'extensió original
            file_extension <- tools::file_ext(temporal_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "merged", paste0("MERGED", file_extension))
            
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(temporal_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(temporal_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(temporal_file_path, ".bim"), paste0(new_name, ".bim"))
            
            # Retorna el nou camí del fitxer processat
            processed_file_path <- new_name
            
            # recover codes from first column
            df2_reactive <- reactive({
              fam <- fread(paste0(file.path(folder_path),"/merged/MERGED.fam"), header = F)
              fam.1 <- fam[, 1:2]
              colnames(fam.1) = c("IDI", "IDII")
              fam.2 <- fam.1[, 1]
              colnames(fam.2) <- "V1"
              fam2 <- fam.2 %>%
                separate(V1, into = c("POP1", "POP2"), sep = "_")
              df2 <- as.data.frame(cbind(fam2,fam.2))
              return(df2)
            })
            
            # UI dinàmic per a la selecció de POP1
            output$pop1_select <- renderUI({
              df2 <- df2_reactive()
              selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
            })
            
            # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
            output$pop2_checkbox <- renderUI({
              req(input$pop1)  # Assegura que POP1 està seleccionat
              df2 <- df2_reactive()
              pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
              selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
            })
          }
          # File processed successfully
          output$status <- renderText("File processed successfully.")
          #    progress$inc(1)  # Increment progress
          
        } else {                                 ################ selection of file loaded by the user
          req(input$upload_file)
          
          usr_file_path <- input$upload_file$datapath
          file_ext <- tools::file_ext(usr_file_path)
          
          base_name <- NULL
          
          if(input$selector == "wwp") {             ########### if loaded file must be merged with WWW 
            
            if (file_ext == "zip") {
              unzip(usr_file_path, exdir = file.path(folder_path, "decompressed_files")) #>>>>>> send to "decompressed_files" folder
              base_name <- tools::file_path_sans_ext(input$upload_file$name)
            } else if (file_ext == "gz") {
              base_name <- tools::file_path_sans_ext(basename(usr_file_path))
              system(paste("gunzip", usr_file_path))
            } else {
              stop("Unsupported file type.")
            }
            decompressed_file_path <- paste0(file.path(folder_path, "decompressed_files", base_name))
            processed_file_path(decompressed_file_path)
            print(base_name)
            # to assign new code to .fam file
            user.fam <- fread(paste0(decompressed_file_path,".fam"), header = F)
            colnames(user.fam) <- c("V1","V2","V3","V4","V5","V6")
            user.id <- paste0(input$popID,"_",input$popID)
            user.fam$V1 <- user.id
            write.table(user.fam, paste0(decompressed_file_path,".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
            
            # Renomena amb el nom "INPUT" mantenint l'extensió original
            file_extension <- tools::file_ext(decompressed_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "decompressed_files", paste0("INPUT", file_extension))
            print(new_name)
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(decompressed_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(decompressed_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(decompressed_file_path, ".bim"), paste0(new_name, ".bim"))
            
          } else {                                ############  if loaded file is processed alone
            if (file_ext == "zip") {
              unzip(usr_file_path, exdir = file.path(folder_path, "merged"))  #>>>>>> send to "merged" folder
              base_name <- tools::file_path_sans_ext(input$upload_file$name)
            } else if (file_ext == "gz") {
              base_name <- tools::file_path_sans_ext(basename(usr_file_path))
              system(paste("gunzip", usr_file_path))
            } else {
              stop("Unsupported file type.")
            }
            merged_file_path <- paste0(file.path(folder_path, "merged", base_name))
            processed_file_path(merged_file_path)
            
            # to assign new code to .fam file
            user.fam <- fread(paste0(merged_file_path,".fam"), header = F)
            colnames(user.fam) <- c("V1","V2","V3","V4","V5","V6")
            user.id <- paste0(input$popID,"_",input$popID)
            user.fam$V1 <- user.id
            write.table(user.fam, paste0(merged_file_path,".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
            
            # Defineix el nou camí amb el nom "MERGED" mantenint l'extensió original
            file_extension <- tools::file_ext(merged_file_path)  # Obtenim l'extensió
            new_name <- file.path(folder_path, "merged", paste0("MERGED", file_extension))
            
            # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
            file.rename(paste0(merged_file_path, ".bed"), paste0(new_name, ".bed"))
            file.rename(paste0(merged_file_path, ".fam"), paste0(new_name, ".fam"))
            file.rename(paste0(merged_file_path, ".bim"), paste0(new_name, ".bim"))
            
            # recover codes from first column
            df2_reactive <- reactive({
              fam <- fread(paste0(file.path(folder_path, "merged/MERGED.fam")), header = F)
              fam.1 <- fam[, 1:2]
              colnames(fam.1) = c("IDI", "IDII")
              fam.2 <- fam.1[, 1]
              colnames(fam.2) <- "V1"
              fam2 <- fam.2 %>%
                separate(V1, into = c("POP1", "POP2"), sep = "_")
              df2 <- as.data.frame(cbind(fam2,fam.2))
              return(df2)
            })
            
            # UI dinàmic per a la selecció de POP1
            output$pop1_select <- renderUI({
              df2 <- df2_reactive()
              selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
            })
            
            # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
            output$pop2_checkbox <- renderUI({
              req(input$pop1)  # Assegura que POP1 està seleccionat
              df2 <- df2_reactive()
              pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
              selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
            })
          }
          # Final step - Update status and progress bar
          output$status <- renderText("File processed successfully.")
          progress$inc(1)  # Ensure progress is fully completed 
        }
##################################### Merge step ##################################  
        # Obtain procesed file
        decompressed_file_path <- paste0(file.path(folder_path, "decompressed_files/INPUT"))
        processed_file_path(decompressed_file_path)
        bfile_path <- processed_file_path()
       # print(bfile_path)
        
        merge_file_path <- paste0(file.path(folder_path, "merged"))

        # Initialize progress bar
        progress <- shiny::Progress$new()
        progress$set(message = "Merging files...", value = 0)
        on.exit({
          progress$close()  # Ensure progress bar is closed at the end of processing
        }, add = TRUE)
       # on.exit(progress$close())  # Ensure progress bar closes on exit
        
        # Function to update progress
        update_progress <- function(amount) {
          progress$inc(amount)
        }
        
        # Set status text to running
        output$status <- renderText("Files merged")
        
        user.bim <- fread(paste0(bfile_path,".bim"))
        # extract from reference those comon to user
        system(paste0("www/plink19 --bfile REFERENCE/Reference --extract ",bfile_path,".bim --aec --make-bed --out ",merge_file_path,"/to_merge1"))
        # extract from user those comon to reference-extracted 
        system(paste0("www/plink19 --bfile ",bfile_path," --extract ",merge_file_path,"/to_merge1.bim --aec --make-bed --out ",merge_file_path,"/to_merge2"))
        
        # merge step 1
        system(paste0("www/plink19 --bfile ",merge_file_path,"/to_merge1 --bmerge ",merge_file_path,"/to_merge2 --aec --make-bed --out ",merge_file_path,"/MERGE1"))
        # Check if ".missnp" file exist
        # Set the file path
        miss_file_path <- paste0(merge_file_path,"/MERGE1-merge.missnp")
        
        # Check if the file exists
        if (file.exists(miss_file_path)) { # if exist procede as:
          # merge step 1
          system(paste0("www/plink19 --bfile  ",merge_file_path,"/to_merge1 --bmerge temporal/to_merge2 --aec --make-bed --out ",merge_file_path,"/MERGE1"))
          # exclude step
          system(paste0("www/plink19 --bfile  ",merge_file_path,"/to_merge1 --exclude  ",merge_file_path,"/MERGE1-merge.missnp --aec --make-bed --out ",merge_file_path,"/preMERGE_a"))
          system(paste0("www/plink19 --bfile  ",merge_file_path,"/to_merge2 --exclude  ",merge_file_path,"/MERGE1-merge.missnp --aec --make-bed --out ",merge_file_path,"/preMERGE_b"))
          system(paste0("www/plink19 --bfile  ",merge_file_path,"/preMERGE_a --bmerge  ",merge_file_path,"/preMERGE_b --aec --make-bed --out ",merge_file_path,"/MERGED"))
          # merge step 2 
        } else {                      # if not exist procede as:
          # merge step 2    
          system(paste0("www/plink19 --bfile  ",merge_file_path,"/to_merge1 --bmerge  ",merge_file_path,"/to_merge2 --aec --make-bed --out ",merge_file_path,"/MERGED"))
        }
        
        # generate dataframe to selector populations
        df2_reactive <- reactive({
          fam <- fread(paste0(merge_file_path,"/MERGED.fam"), header = F)
          fam.1 <- fam[, 1:2]
          colnames(fam.1) = c("IDI", "IDII")
          fam.2 <- fam.1[, 1]
          colnames(fam.2) <- "V1"
          fam2 <- fam.2 %>%
            separate(V1, into = c("POP1", "POP2"), sep = "_")
          df2 <- as.data.frame(cbind(fam2,fam.2))
          return(df2)
        })
        
        # UI dinàmic per a la selecció de POP1
        output$pop1_select <- renderUI({
          df2 <- df2_reactive()
          selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
        })
        
        # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
        output$pop2_checkbox <- renderUI({
          req(input$pop1)  # Assegura que POP1 està seleccionat
          df2 <- df2_reactive()
          pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
          selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
          # checkboxGroupInput("pop2", "Select POP2 Groups:", choices = unique(pop2_choices))
        }) 
        shinyjs::enable("run_pca")
        showNotification("Files merged successfully! You can now perform PCA analyses by selecting target populations.", type = "message", duration = 10)
      })

      #####################################################################################
      ######################################PCA ###########################################
      #####################################################################################
      observeEvent(input$run_pca, {
        # Reset PLINK status to FALSE when starting
      #  run_pca(FALSE)
        
        req(processed_file_path())  # Assegurar que el fitxer es va processar
        req(input$pop1)  # Assegurar que hi ha alguna superpoblació seleccionada
        
        # Inicialitzar el bar de progrés
        progress <- shiny::Progress$new()
        progress$set(message = "Running PCA by PLINK...", value = 0)
       # on.exit(progress$close())  # Assegurar que es tanca el bar
        on.exit({
          progress$close()  # Ensure progress bar is closed at the end of processing
        }, add = TRUE)
        
        # Funció per actualitzar el bar de progrés
        update_progress <- function(amount) {
          progress$inc(amount)
        }
        
        # Obtenir el fitxer processat
        merge_file_path <- paste0(file.path(folder_path, "merged"))
        bfile_path <- paste0(merge_file_path,"/MERGED")
        pca_files <- paste0(file.path(folder_path, "pca"))
        
        # Comprovar que els fitxers necessaris existeixen
        required_files <- paste0(bfile_path, c(".bed", ".bim", ".fam"))
        if (!all(file.exists(required_files))) {
          showNotification("Missing required PLINK files (.bed, .bim, .fam).", type = "error", duration = 10)
          return(NULL)
        }
        
        # Executar PLINK, assegurant-se que es mostra el missatge d'error si falla
        pop.pre1 <- fread(paste0(bfile_path,".fam"), header = F)
        pop.pre1.1 <- pop.pre1[, 1:2]
        colnames(pop.pre1.1) = c("IDI", "IDII")
        pop.pre2 <- pop.pre1[, 1]
        colnames(pop.pre2) <- "V1"
        pop2 <- pop.pre2 %>%
          separate(V1, into = c("POP1", "POP2"), sep = "_")
        pop2 <- cbind(pop2,pop.pre1.1)
        write.table(pop2,paste0(pca_files,"/pop_data.txt"), quote = FALSE, row.names = FALSE, col.names = T)
        
        to_keep <- subset(pop2, POP1 %in% input$pop1)
        write.table(to_keep[, 3:4], paste0(pca_files,"/to_keep.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
        write.table(to_keep[, 1:2], paste0(pca_files,"/to_color.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        plink_result <-  tryCatch({
          # Step 0: Filter the data by selected superpopulation(s) unless "Select all" is checked
          if (input$select_all_pop1=="Select all") {
            # If "Select All" is checked, skip the --keep step
            system(paste0("www/plink19 --bfile", bfile_path, " --make-bed --out ",pca_files,"/step1"))
            update_progress(0.20)  # Update progress to 20%
          } else {
            pop.pre1 <- fread(paste0(pca_files,"/to_keep.txt"), header = F)
            system(paste0("www/plink19 --bfile ", bfile_path, " --keep ",pca_files,"/to_keep.txt --make-bed --out  ",pca_files,"/step1"))
            update_progress(0.20)  # Update progress to 20%
          } 
          # Step 1: Run PLINK with MAF filter
          n = input$sliderM1
          print(n)
          system(paste0("www/plink19 --bfile  ",pca_files,"/step1 --maf ", n ," --make-bed --out  ",pca_files,"/step2"))
          update_progress(0.20)  # Update progress to 40%
          
          # Step 2: Perform linkage disequilibrium (LD) pruning
          system(paste0("www/plink19 --bfile  ",pca_files,"/step2 --indep-pairwise 200 25 0.3 --out  ",pca_files,"/prune1"))
          update_progress(0.20)  # Update progress to 60%
          
          # Step 3: Extract pruned SNPs and create a new binary file
          system(paste0("www/plink19 --bfile  ",pca_files,"/step2 --extract  ",pca_files,"/prune1.prune.in --make-bed --out  ",pca_files,"/step3"))
          update_progress(0.20)  # Update progress to 80%
          
          # Step 4: Perform PCA
          system(paste0("www/plink19 --bfile  ",pca_files,"/step3 --pca 10 header tabs var-wts --out  ",pca_files,"/input_PCA"))
          update_progress(0.20)  # Update progress to 100%
          
          output$status <- renderText("PCA analysis completed successfully.")
          
        }, error = function(e) {
          showNotification("PLINK analysis failed. Check the logs.", type = "error", duration = 10)
          output$status <- renderText(paste("Error:", e$message))
        })
############################################ PCA log file
        # Directory for temporary files
        log_file <- file.path(pca_files, "step3.log")
        
        # Display PLINK output in the UI
        output$plink_output <- renderText({
          paste(plink_result, collapse = "\n")
        })
        
        # Wait until the log file is generated
        if (file.exists(log_file)) {
          # Allow user to search log file dynamically
          keyword <- "variants and"
          lines_before <- 0
          lines_after <- 0
          
          # Construct grep command
          grep_command <- paste0(
            "grep ",
            ifelse(lines_before > 0, paste0("-B ", lines_before, " "), ""),
            ifelse(lines_after > 0, paste0("-A ", lines_after, " "), ""),
            shQuote(keyword), " ", shQuote(log_file)
          )
       
          # Extract matching lines from the log file
          log_lines <- tryCatch(
            system(grep_command, intern = TRUE),
            error = function(e) "No matching lines found."
          )
          
          # Display the extracted lines in the UI
          output$log_output <- renderText({
            paste(log_lines, collapse = "\n")
          })
        } else {
          output$log_output <- renderText({
            "Log file not found. Ensure PLINK ran successfully."
          })
        }
   
        ############################################
        
        update_progress(1)  # Finalitzar el progrés
        
        # Read the PCA results from the PLINK output file
        PCA_results <- fread(paste0(pca_files,"/input_PCA.eigenvec"), header = TRUE)
        
        
        # Render the PCA plot
        output$PCA_plot <- renderPlot({
          to_color <- fread(paste0(pca_files,"/to_color.txt"), header = F)
          colnames(to_color) <- c("POP1", "POP2")
          
          # Get unique items in the order they first appear
          unique_names <- unique(unique(to_color$POP1))
          # Reorder the Name column based on unique values
          to_color$POP1 <- factor(to_color$POP1, levels = unique_names)
          # Get unique items in the order they first appear
          unique_names2 <- unique(unique(to_color$POP2))
          # Reorder the Name column based on unique values
          to_color$POP2 <- factor(to_color$POP2, levels = unique_names2)
          
          ploted_files <- paste0(file.path(folder_path, "plots"))
          
          if (input$popcolor == "POP1") {
            PCA <- ggplot(PCA_results, aes(x = PC1, y = PC2, color = factor(to_color$POP1))) +
              geom_point(size = 3, shape=1, stroke = 1.5) +
              labs(title = "PCA Plot", x = "PC1", y = "PC2") +
              guides(color = guide_legend(title = "Superpopulation")) 
          } else if (input$popcolor == "POP2") {
            PCA <-  ggplot(PCA_results, aes(x = PC1, y = PC2, color = factor(to_color$POP2), shape = factor(to_color$POP1))) +
              geom_point(size = 3, stroke = 1.5) +
              labs(title = "PCA Plot", x = "PC1", y = "PC2") +
              guides(color = guide_legend(title = "Population"),
                     shape = guide_legend(title = "Superpulation")) 
          }
          plot(PCA)
          ggsave(paste0(ploted_files,"/PCA.png"), plot = PCA)
        })
        
        shinyjs::enable("run_admix")
        shinyjs::enable("run_fst")
        showNotification("PCA analysis finished. You can now run ADMIXTURE or FST analysis.", type = "message", duration = 10)
        # Update the status to TRUE when done
     #   run_pca(TRUE)
      })  
      
      #####################################################################################
      ############################### ADMIXTURE ###########################################
      #####################################################################################
      observeEvent(input$run_admix, {
        req(input$run_pca)  # Ensure PLINK is run first
        
        # Initialize progress bar
        progress <- shiny::Progress$new()
        progress$set(message = "Running ADMIXTURE analysis. This may take a while......", value = 0)
       # on.exit(progress$close())  # Ensure progress bar closes on exit
        on.exit({
          progress$close()  # Ensure progress bar is closed at the end of processing
        }, add = TRUE)
        
        # Function to update progress
        update_progress <- function(amount) {
          progress$inc(amount)
        }
        
        # Set status text to running
        output$status <- renderText("Running ADMIXTURE analysis...")
        
        admx_files <- paste0(file.path(folder_path, "admx"))
        pca_files <- paste0(file.path(folder_path, "pca"))
        
        # adjust MAF
        m = input$sliderM2
        MAF_filter <-  tryCatch({
        system(paste0("www/plink19 --bfile ",pca_files,"/step3 --maf ", m, " --make-bed --out ",admx_files,"/input_admx"))  
        update_progress(0.05 / (input$sliderM2 - 0.05))  # Update progress incrementally based on M
        })
        ############################################ to show maf log file
        # Directory for temporary files
        log_file <- file.path(admx_files, "input_admx.log")
        
        # Display PLINK output in the UI
        output$maf_output <- renderText({
          paste(MAF_filter, collapse = "\n")
        })
        
        # Wait until the log file is generated
        if (file.exists(log_file)) {
          # Allow user to search log file dynamically
          keyword <- "variants and"
          lines_before <- 0
          lines_after <- 0
          
          # Construct grep command
          grep_command <- paste0(
            "grep ",
            ifelse(lines_before > 0, paste0("-B ", lines_before, " "), ""),
            ifelse(lines_after > 0, paste0("-A ", lines_after, " "), ""),
            shQuote(keyword), " ", shQuote(log_file)
          )
          
          # Extract matching lines from the log file
          log_lines <- tryCatch(
            system(grep_command, intern = TRUE),
            error = function(e) "No matching lines found."
          )
          
          # Display the extracted lines in the UI
          output$maf_output <- renderText({
            paste(log_lines, collapse = "\n")
          })
        } else {
          output$maf_output <- renderText({
            "Log file not found. Ensure PLINK ran successfully."
          })
        }
        
        ############################################
        
        # detect number of processor Np. Use Np-1 
        Np <- parallel::detectCores() - 1
        print(Np)
        
        # Run ADMIXTURE for the selected K value
        for (i in 1:input$sliderK) {
          system(paste("www/admixture --cv", paste0("-j",Np," ",admx_files,"/input_admx.bed"), i, "| tee", paste0("log", i, ".out"))) 
          update_progress(1.0 / (input$sliderK - 1))  # Update progress incrementally based on K
        }
        
        # Identify the .Q file and move to temporary admx folder
        main_dir <- getwd()  # Get the current working directory
        
        for (j in 1:10) {  # Loop through the specified range
          q_file_name <- paste0("input_admx.", j, ".Q")  # Construct the file name
          q_file_path <- file.path(main_dir, q_file_name)  # Create the full path for the source file
          destination_path <- file.path(admx_files, q_file_name)  # Create the destination path

          # Move the .Q file if it exists
          if (file.exists(q_file_path)) {  
            success <- file.rename(q_file_path, destination_path)  # Attempt to rename (move) the file
            if (success) {
              cat(paste("Moved .Q file to:", destination_path, "\n"))  # Success message
            } else {
              cat(paste("Failed to move .Q file from", q_file_path, "to", destination_path, "\n"))  # Failure message
            }
          } else {
            cat(paste(".Q file does not exist at:", q_file_path, "\n"))  # File not found message
          }
        }
        
        # Identify the log file and move to temporary admx folder
        for (t in 1:10) {  # Loop through the specified range
          log_file_name <- paste0("log", t, ".out")  # Construct the file name
          log_file_path <- file.path(main_dir, log_file_name)  # Create the full path for the source file
          destination_path <- file.path(admx_files, log_file_name)  # Create the destination path
          
          # Move the .Q file if it exists
          if (file.exists(log_file_path)) {  
            success <- file.rename(log_file_path, destination_path)  # Attempt to rename (move) the file
            if (success) {
              cat(paste("Moved log file to:", destination_path, "\n"))  # Success message
            } else {
              cat(paste("Failed to move log file from", log_file_path, "to", destination_path, "\n"))  # Failure message
            }
          } else {
            cat(paste("log file does not exist at:", log_file_path, "\n"))  # File not found message
          }
        }
        
        # Initialize a list to store the log contents
        all_logs <- list()  # Initialize a list to store logs
        
        # Loop through the range of K values
        for (r in 1:input$sliderK) {
          log_file <- paste0(admx_files, "/log", r, ".out")  # Construct the log file path
          
          # Check if the file exists before reading it
          if (file.exists(log_file)) {
            adm_log <- readLines(log_file)  # Read the log file lines
            all_logs[[paste0("K", r)]] <- adm_log  # Store log content
          } else {
            warning(paste("Log file not found:", log_file))  # Warning if file not found
          }
        }
        
        # Extract lines with "CV error" from all logs
        cv_error_lines <- unlist(lapply(all_logs, function(log) grep("CV error", log, value = TRUE)))
       # print(cv_error_lines)
        # Output the CV error lines to the verbatimTextOutput
        output$cvErrors <- renderPrint({
          if (length(cv_error_lines) > 0) {
            cat(cv_error_lines, sep = "\n")  # Print all CV error lines
          } else {
            cat("No CV error lines found.")  # Message if no lines found
          }
        })

       # remove undesired files
        unlink("1")
        unlink("2")
        unlink("3")
        unlink("4")
        unlink("5")
        unlink("6")
        unlink("7")
        unlink("8")
        unlink("9")
        unlink("10")
        unlink("log")
        unlink("*.P")
        unlink("log*")
        unlink("*.out")
        
        # UI dinàmic per a la selecció del admixture component
        output$sorter <- renderUI({
          selectInput("sorter", "Sort admixture boxplot by 'Ancestry component': ", choices = 1:input$sliderK, multiple = F,selected=1)
        }) 
        
        # Render the ADMIXTURE plot
        output$ADMIXTURE_plot <- renderPlot({
          # Load data
          runs <- fread(paste0(admx_files,"/input_admx.", input$sliderK, ".Q"))
          to_color <- fread(paste0(pca_files,"/to_color.txt"), header = F)
          colnames(to_color) <- c("POP1", "POP2")
          
          # Combine population and ancestry data
          admixture_data <- cbind(to_color, runs)
          ##############################################################################
          data <- admixture_data
          # Ordenar el dataframe per la primera columna després de les columnes de població
          if (input$popcolor == "POP1"){
            data_sorted <- data %>%
              group_by(POP1) %>%
              arrange(V1, .by_group = TRUE) %>%
              ungroup()
            
            # Afegir una columna d'índex per a cada grup
            data_sorted <- data_sorted %>%
              group_by(POP1) %>%
              mutate(Index = row_number()) %>%
              ungroup() 
            
            # Transformar el dataframe a format llarg
            data_long <- data_sorted %>%
              pivot_longer(cols = starts_with("V"), names_to = "Variable", values_to = "Value")
            
            # Crear el gràfic
            ADM <- ggplot(data_long, aes(x = Index, y = Value, fill = Variable)) +
              geom_bar(stat = "identity", position = "stack") +
              facet_nested(~POP1, scales = "free_x") +
              labs(x = "Individuals", y = "Ancestry component", title = "Admixture Plot",fill = "Components") +
              #  theme_minimal() +
              theme(legend.position = "bottom")
            
          } else {
            data_sorted <- data %>%
              group_by(POP2) %>%
              arrange(V1, .by_group = TRUE) %>%
              ungroup()
            
            # Afegir una columna d'índex per a cada grup
            data_sorted <- data_sorted %>%
              group_by(POP2) %>%
              mutate(Index = row_number()) %>%
              ungroup() 
            
            # Transformar el dataframe a format llarg
            data_long <- data_sorted %>%
              pivot_longer(cols = starts_with("V"), names_to = "Variable", values_to = "Value")
            
            ADM <- ggplot(data_long, aes(x = Index, y = Value, fill = Variable)) +
              geom_bar(stat = "identity", position = "stack") +
              facet_nested(~POP1 + POP2, scales = "free_x", space = "free") +
              labs(
                x = "Individuals",
                y = "Ancestry component",
                title = "Admixture Plot",
                fill = "Components"
              ) +
              # theme_minimal() +
              theme(
                legend.position = "bottom",
                plot.margin = margin(10, 10, 10, 10), # Adjust as needed (top, right, bottom, left)
                panel.spacing = unit(0, "lines"),
                strip.text.y = element_text(angle = 45) # Adjust text rotation for row labels
              )
            
          }  
          print(ADM)
          ploted_files <- paste0(file.path(folder_path, "plots"))
          
          # Save the plot as PNG
          ggsave(paste0(ploted_files,"/ADM.png"), plot = ADM)
##############################################################################
        })
        output$box_plot <- renderPlot({
          # Load data
          runs <- fread(paste0(admx_files,"/input_admx.", input$sliderK, ".Q"))
          to_color <- fread(paste0(pca_files,"/to_color.txt"), header = F)
          colnames(to_color) <- c("POP1", "POP2")
          
          # Combine population and ancestry data
          admixture_data <- cbind(to_color, runs)
          # Calculate mean values of first (V1) component by POP2
          for (i in 1:input$sorter) { 
            ii <- paste0("V",i)
            MEAN_POP2 <- aggregate(admixture_data[[ii]], by = list(admixture_data$POP2), FUN = mean)
          }
          #order by mean value
          MEAN_POP2 = MEAN_POP2[order(MEAN_POP2$x),]
          MEAN_POP2$Group.1
          level_orderX <- c(MEAN_POP2$Group.1)
          admixture_data <- merge(admixture_data, MEAN_POP2, by.x = "POP2", by.y = "Group.1", all.x = TRUE)
          
          # Sort admixture data based on the selected population column
          admixture_data <- admixture_data[order(admixture_data[["x"]]), ]
          admixture_data <- admixture_data %>%
            relocate(x, .before = V1)
          # Reshape the data for ggplot (long format)
          admixture_long <- pivot_longer(admixture_data, cols = 4:ncol(admixture_data))
          
          if (input$popcolor == "POP1") {
            p <- ggplot(admixture_long, aes(x = name, y = value, fill = POP1)) + 
              geom_boxplot() +
              labs(x = "Ancestry component", y = "Ancestry Proportion", fill = "Population") 
            
          } else {
            
            p <- ggplot(admixture_long, aes(x = name, y = value, fill = factor(POP2, level = level_orderX))) + 
              geom_boxplot() +
              facet_grid(cols = vars(POP1),scales = "free", space = "free") +
              scale_x_discrete(expand = c(0, 0.5)) +
              labs(x = "Ancestry component", y = "Ancestry Proportion", fill = "Population") +
              ggforce::facet_row(vars(POP1), scales = 'free', space = 'free')
            
          }    
          
          print(p)
          ploted_files <- paste0(file.path(folder_path, "plots"))
          
          # Save the plot as PNG
          ggsave(paste0(ploted_files,"/box_ADM.png"), plot = p)
        })
        ################################################################################     
        output$map_plot <- renderPlot({
          # Load data
          runs <- fread(paste0(admx_files,"/input_admx.", input$sliderK, ".Q"))
          to_color <- fread(paste0(pca_files,"/to_color.txt"), header = F)
          colnames(to_color) <- c("POP1", "POP2")
          coord <- fread("www/coord.txt", header = TRUE)
          
          # Add user-defined data if applicable
          if (input$files == "usr") {
            user_data <- data.frame(
              POP1 = input$popID,
              POP2 = input$popID,
              Pop_name = "user population",
              LAT = as.numeric(input$lat),
              LON = as.numeric(input$lon)
            )
            coord <- rbind(coord, user_data)
          }
          
          # Count population sizes
          pop_counts <- to_color %>%
            count(POP2, name = "pop_size")
          
          # Merge data
          admixture_data <- merge(to_color, pop_counts, by = "POP2") %>%
            cbind(runs) %>%
            merge(coord, by = c("POP1", "POP2"))
          
          Npop <- nrow(pop_counts)
          qpop <- matrix(NA, ncol = input$sliderK, nrow = Npop)
          
          # Calculate mean ancestry proportions and coordinates for each population
          coord_pop <- matrix(NA, ncol = 2, nrow = Npop)
          for (i in seq_len(Npop)) {
            pop_subset <- admixture_data$POP2 == pop_counts$POP2[i]
            qpop[i, ] <- colMeans(runs[pop_subset, ])
            coord_pop[i, ] <- colMeans(admixture_data[pop_subset, c("LON", "LAT")])
          }
          
          # Define plot limits
          xlim <- range(coord_pop[, 1]) + c(-1, 1) * input$sliderE * diff(range(coord_pop[, 1]))
          ylim <- range(coord_pop[, 2]) + c(-1, 1) * input$sliderE * diff(range(coord_pop[, 2]))
          
          # Create base map
          basemap(xlim, ylim)
          plot(coord_pop, xlab = "LON", ylab = "LAT", type = "n", xlim = xlim, ylim = ylim)
          map(add = TRUE, col = "grey90", fill = TRUE)
          
          # Scale pie chart radius based on population size
          max_pop_size <- max(pop_counts$pop_size)
          radius_scaling_factor <- input$sliderR / 20
          
          for (i in seq_len(Npop)) {
            scaled_radius <- sqrt(pop_counts$pop_size[i] / max_pop_size) * radius_scaling_factor
            add.pie(
              z = qpop[i, ],
              x = coord_pop[i, 1],
              y = coord_pop[i, 2],
              labels = "",
              radius = scaled_radius,
              col = rainbow(input$sliderK)
            )
            # Add POP2 labels
            text(
              x = coord_pop[i, 1],
              y = coord_pop[i, 2] + scaled_radius * 1.1,  # Offset vertically
              labels = pop_counts$POP2[i],
              cex = 0.7,  # Adjust size as needed
              col = "black",
              font = 2  # Make the font bold
            )
          }
          
          width <- 2400  # Width in pixels
          height <- 1400 # height in pixels

          ploted_files <- paste0(file.path(folder_path, "plots"))
          
          png(paste0(ploted_files,"/map_ADM.png"), width = width, height = height, res = 150)
          
          # Replot within the PNG device
          basemap(xlim, ylim)
          plot(coord_pop, xlab = "LON", ylab = "LAT", type = "n", xlim = xlim, ylim = ylim)
          map(add = TRUE, col = "grey90", fill = TRUE)
          
          for (i in seq_len(Npop)) {
            scaled_radius <- sqrt(pop_counts$pop_size[i] / max_pop_size) * radius_scaling_factor
            add.pie(
              z = qpop[i, ],
              x = coord_pop[i, 1],
              y = coord_pop[i, 2],
              labels = "",
              radius = scaled_radius,
              col = rainbow(input$sliderK)
            )
            text(
              x = coord_pop[i, 1],
              y = coord_pop[i, 2] + scaled_radius * 1.5,  # Offset vertically
              labels = pop_counts$POP2[i],
              cex = 0.9,  # Adjust size for output
              col = "black",
              font = 2
            )
          }
          
          dev.off()
          # }
        })
      })
#######################################################################################
####################################### FST ###########################################
#######################################################################################
      # Run GCTA Fst test 
      observeEvent(input$run_fst, {
        req(input$run_pca)  # Ensure PLINK is run first
        
        # Initialize progress bar
        progress <- shiny::Progress$new()
        progress$set(message = "Running FST analysis", value = 0)
        # on.exit(progress$close())  # Ensure progress bar closes on exit
        on.exit({
          progress$close()  # Ensure progress bar is closed at the end of processing
        }, add = TRUE)
        
        pca_files <- paste0(file.path(folder_path, "pca"))
        b_files <- paste0(file.path(folder_path, "pca/step3"))
        fst_files <- paste0(file.path(folder_path, "fst"))
        fst_in_files <- paste0(file.path(folder_path, "fst/to_fst.txt"))
        fst_out_files <- paste0(file.path(folder_path, "fst/FST_output"))
        
        # Step 0: Filter the data by selected superpopulation(s) unless "Select All" is checked
        if (input$select_all_pop2=="Select all") {
          pop2 <- fread(paste0(pca_files,"/pop_data.txt"), header = TRUE)
          to_fst <- to_fst[, c(3,4,2)]
          write.table(to_fst, paste0(fst_files,"/to_fst.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
          system(paste("www/gcta64 --bfile ",b_files," --fst --sub-popu ",fst_in_files," --out ",fst_out_files))
        } else {
          pop2 <- fread(paste0(pca_files,"/pop_data.txt"), header = TRUE)
          to_fst <- subset(pop2, POP2 %in% input$pop2)
          to_fst <- to_fst[, c(3,4,2)]
          write.table(to_fst, paste0(fst_files,"/to_fst.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
          system(paste("www/gcta64 --bfile ",b_files," --fst --sub-popu ",fst_in_files," --out ",fst_out_files))
        } 
        
        # Generate and render the Manhattan plot
        output$manhattan_plot <- renderPlot({
          pop_names <- fread(fst_in_files, header = F)
          names_pop <- unique(pop_names$V3)
          FST_output <- fread(paste0(fst_out_files,".fst"), header = TRUE)
          FST_output <- FST_output[complete.cases(FST_output[, c("Chr", "bp", "Fst")]), ]
          FST_output <- FST_output[is.finite(FST_output$Fst), ]
          
          # Convert necessary columns to numeric if they aren't already
          FST_output$Chr <- as.numeric(as.character(FST_output$Chr))
          FST_output$bp <- as.numeric(as.character(FST_output$bp))
          FST_output$Fst <- as.numeric(as.character(FST_output$Fst))
          FST_output <- as.data.frame(FST_output)
          # Find the row with the maximum Fst value
          max_fst_row <- FST_output[which.max(FST_output$Fst), ]
          print(max_fst_row)
          
          
            output$max_fst_output <- renderPrint({
              if (nrow(max_fst_row) > 0) {
                print(max_fst_row)  # Print the max FST row
              } else {
                cat("No FST data available.")  # Message if no data available
              }
            })

          # Print the SNP value of the row with the maximum Fst
          max_fst_snp <- max_fst_row$SNP
          max_fst_snp
          max1 <- max(FST_output$Fst)
          max <- max(FST_output$Fst) + (max(FST_output$Fst)*0.1)
          qqman::manhattan(
            FST_output,
            chr = "Chr",
            bp = "bp",
            p = "Fst",
            snp = "SNP",
            highlight = max_fst_snp,
            col = c("gray10", "red"),
            chrlabs = c(1:23),
            annotateTop = TRUE,
            ylim = c(0,max),
            annotatePval = max1,
            ylab = "Fst",
            main = paste0("Manhattan plot of Fst-scores of population: ",names_pop),
            logp = F)
          # to save as png
          ploted_files <- paste0(file.path(folder_path, "plots"))
          png(paste0(ploted_files,"/FST.png"), width = 2400, height = 1200, res = 150)
          qqman::manhattan(
            FST_output,
            chr = "Chr",
            bp = "bp",
            p = "Fst",
            snp = "SNP",
            highlight = max_fst_snp,
            col = c("gray10", "red"),
            chrlabs = c(1:23),
            annotateTop = TRUE,
            ylim = c(0,max),
            annotatePval = max1,
            ylab = "Fst",
            main = paste0("Manhattan plot of Fst-scores of population: ",names_pop),
            logp = F)
          dev.off()
        })
      })
      
###########################################
  
      # Download report
      output$downloadReport <- downloadHandler(
        filename = function() {
          paste("Analysis_Report_", Sys.Date(), ".html", sep = "")
        },
        content = function(file) {
          # Path to plots
          ploted_files <- paste0(file.path(folder_path, "plots"))
          # Generate a report using Rmarkdown
          rmarkdown::render("report_template.Rmd", output_file = file,
                            params = list(
                              pca_plot = paste0(ploted_files,"/PCA.png"),
                              admixture_plot = paste0(ploted_files,"/ADM.png"),
                              box_ADM = paste0(ploted_files,"/box_ADM.png"),
                              map_ADM = paste0(ploted_files,"/map_ADM.png"),
                              manhattan_plot = paste0(ploted_files,"/FST.png")
                            ))
        }
      )  
  
      ###########################################
      observeEvent(input$clear_folders, {
        base_folder <- file.path(folder_path) # Ensure folder_path is valid
        
        # Check if base_folder is valid
        if (!dir.exists(base_folder)) {
          output$clear_status <- renderText({
            paste("Base folder does not exist:", base_folder)
          })
          return() # Stop execution if base folder does not exist
        }
        
        subfolders <- c("admx","temporal", "merged", "pca", "plots", "fst", "decompressed_files","temporal")
        status <- c()
        
        for (subfolder in subfolders) {
          subfolder_path <- file.path(base_folder, subfolder)
          
          if (dir.exists(subfolder_path)) {
            # Remove all files in the subfolder
            unlink(file.path(subfolder_path, "*"), recursive = TRUE, force = TRUE)
            status <- c(status, paste("Cleared files in:", subfolder_path))
          } else {
            status <- c(status, paste("Subfolder does not exist:", subfolder_path))
          }
        }
        
        # Render the status messages
        output$clear_status <- renderText({
          paste(status, collapse = "\n")
        })
      })
      
      
      observeEvent(input$reload_button, {
        session$reload()
      })
  }) # end

  ###########################################
}
# Reset All button functionality

shinyApp(ui = ui, server = server)

