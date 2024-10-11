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

# Increase max upload size to 500MB
options(shiny.maxRequestSize = 500 * 1024^2)

ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  titlePanel("Basic Population Genetic Analysis - (BPGA)"),
  sidebarLayout(
    sidebarPanel(
    #  actionButton("remove", "Remove previous files"),
      actionButton("reset_all", "Reset All"),
      actionButton("reset_session", "Reset session files"),
      selectInput("files", "Select input file:",
                  c("Example file" = "exp",
                    "User file" = "usr")),
      conditionalPanel(
        condition = "input.files == 'usr'",
      fileInput("upload_file", "Load user population data (compressed .zip or .gz)", 
                accept = c(".zip", ".gz")),
      textInput("popID", "Assign population identifier (use a three-letter code like 'USR')", "USR"),
      textInput("lat", "Assign latitude"),
      textInput("lon", "Assign longitude")
      ),
      h4("Process input File"),
      radioButtons("selector", "Select option:", 
              choices = list("Merge with Woldwide populations" = "wwp","Only loaded population" = "usr")),
      actionButton("process_file", "Decompress and label"),
      conditionalPanel(
         condition = "input.selector == 'wwp'",
      actionButton("run_merge", "Merge files")
      ),
      h4("Perform population analysis"),
      uiOutput("pop1_select"),  # Dynamic SelectInput for POP1
      shinyjs::hidden(actionButton("select_all_pop1", "Select all")),
      tags$hr(style="border-color: grey;"),
 
      h4("PCA analysis"),
      actionButton("run_plink", "Run PCA analysis"),
      h4("Admixture analysis"),
      sliderInput("sliderK", "Set K values:", min = 1, max = 10, value = 1),
      sliderInput("sliderM", "Set MAF:", min = 0.05, max = 0.50, value = 0.05, step = 0.05),
      
      actionButton("run_admix", "Run ADMIXTURE"),
      h4("Decoration plots"),
      radioButtons("popcolor", "Colored:", 
                   choices = list("Superpopulation" = "SUP", "Subpopulation" = "POP")),
      sliderInput("sliderR", "Set Pie radius:", min = 0, max = 10, value = 4, step = 0.5),
      sliderInput("sliderE", "Expand map:", min = 0, max = 2, value = 0.2, step = 0.05),
      uiOutput("sorter"),
      #selectInput("sorter", "Sort admixture plot:",
      #                 c("By first component" = "V1",
      #                   "By second component" = "V2")),
      tags$hr(style="border-color: grey;"),
      h4("Compare two populations"),
      h4("FST analysis"),
      uiOutput("pop2_checkbox"),  # Dynamic Checkbox for POP2
      shinyjs::hidden(actionButton("select_all_pop2", "Remove selected")),
      actionButton("run_fst", "Run FST analysis"),
      
      h4("Report"),
      downloadButton("downloadReport", "Download Report")
    ),
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

server <- function(input, output, session) {
  processed_file_path <- reactiveVal(NULL)  # Store processed file path
  # Disable the PCA, ADMIXTURE, and FST buttons by default
  shinyjs::disable("run_plink")
  shinyjs::disable("run_admix")
  shinyjs::disable("run_fst")
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
      shinyjs::enable("run_plink")
    }

    if (input$files == "exp"){                   ############### if the example file is selected
 
    # go to folder with the example file
    file_path <- "example/example.zip"
    file_ext <- tools::file_ext(file_path)
    base_name <- "example"
 
    if(input$selector == "wwp") {                 ####### if example file is merged with WWW

    unzip(file_path, exdir = "decompressed_files") #>>>>>> send to "decompressed_files" folder
    decompressed_file_path <- paste0("decompressed_files/", base_name)
    # Store the processed file path
    processed_file_path(decompressed_file_path)

    # rename files as INPUT mantenint l'extensió original
    file_extension <- tools::file_ext(decompressed_file_path)  # Obtenim l'extensió
    new_name <- file.path("decompressed_files", paste0("INPUT", file_extension))

    # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
    file.rename(paste0(decompressed_file_path, ".bed"), paste0(new_name, ".bed"))
    file.rename(paste0(decompressed_file_path, ".fam"), paste0(new_name, ".fam"))
    file.rename(paste0(decompressed_file_path, ".bim"), paste0(new_name, ".bim")) 
    # Retorna el nou camí del fitxer processat
    processed_file_path <- new_name

    # recover codes from first column
    df2_reactive <- reactive({
      fam <- fread(paste0("decompressed_files/INPUT.fam"), header = F)
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
      selectInput("pop1", "Select Superpopulation: (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
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
      
      unzip(file_path, exdir = "temporal")   #>>>>>>> send to "temporal" folder
      temporal_file_path <- paste0("temporal/", base_name)

      processed_file_path(temporal_file_path)

      # rename file as "MERGED" mantenint l'extensió original
      file_extension <- tools::file_ext(temporal_file_path)  # Obtenim l'extensió
      new_name <- file.path("temporal", paste0("MERGED", file_extension))

      # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
      file.rename(paste0(temporal_file_path, ".bed"), paste0(new_name, ".bed"))
      file.rename(paste0(temporal_file_path, ".fam"), paste0(new_name, ".fam"))
      file.rename(paste0(temporal_file_path, ".bim"), paste0(new_name, ".bim"))
      
      # Retorna el nou camí del fitxer processat
        processed_file_path <- new_name

      # recover codes from first column
      df2_reactive <- reactive({
        fam <- fread(paste0("temporal/MERGED.fam"), header = F)
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
        selectInput("pop1", "Select Superpopulation: (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
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
    
    } else {                                 ################ selection of file loaded by the user file
    req(input$upload_file)

    file_path <- input$upload_file$datapath
    file_ext <- tools::file_ext(file_path)
    
    base_name <- NULL
    
    if(input$selector == "wwp") {             ########### if loaded file is merged with WWW 

    if (file_ext == "zip") {
      unzip(file_path, exdir = "decompressed_files") #>>>>>> send to "decompressed_files" folder
      base_name <- tools::file_path_sans_ext(input$upload_file$name)
    } else if (file_ext == "gz") {
      base_name <- tools::file_path_sans_ext(basename(file_path))
      system(paste("gunzip", file_path))
    } else {
      stop("Unsupported file type.")
    }
    decompressed_file_path <- paste0("decompressed_files/", base_name)
    processed_file_path(decompressed_file_path)
    
    # to assign new code to .fam file
    user.fam <- fread(paste0(decompressed_file_path,".fam"), header = F)
    colnames(user.fam) <- c("V1","V2","V3","V4","V5","V6")
    user.id <- paste0(input$popID,"_",input$popID)
    user.fam$V1 <- user.id
    write.table(user.fam, paste0(decompressed_file_path,".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
   
    # Renomena amb el nom "INPUT" mantenint l'extensió original
    file_extension <- tools::file_ext(decompressed_file_path)  # Obtenim l'extensió
    new_name <- file.path("decompressed_files", paste0("INPUT", file_extension))
    
    # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
    file.rename(paste0(decompressed_file_path, ".bed"), paste0(new_name, ".bed"))
    file.rename(paste0(decompressed_file_path, ".fam"), paste0(new_name, ".fam"))
    file.rename(paste0(decompressed_file_path, ".bim"), paste0(new_name, ".bim"))
    
     } else {                                ############  if loaded file is processed alone
      if (file_ext == "zip") {
        unzip(file_path, exdir = "temporal")   #>>>>>> send to "temporal" folder
        base_name <- tools::file_path_sans_ext(input$upload_file$name)
      } else if (file_ext == "gz") {
        base_name <- tools::file_path_sans_ext(basename(file_path))
        system(paste("gunzip", file_path))
      } else {
        stop("Unsupported file type.")
      }
       temporal_file_path <- paste0("temporal/", base_name)
     processed_file_path(temporal_file_path)
     
     # to assign new code to .fam file
     user.fam <- fread(paste0(temporal_file_path,".fam"), header = F)
     colnames(user.fam) <- c("V1","V2","V3","V4","V5","V6")
     user.id <- paste0(input$popID,"_",input$popID)
     user.fam$V1 <- user.id
     write.table(user.fam, paste0(temporal_file_path,".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
     
    # Defineix el nou camí amb el nom "MERGED" mantenint l'extensió original
    file_extension <- tools::file_ext(temporal_file_path)  # Obtenim l'extensió
    new_name <- file.path("temporal", paste0("MERGED", file_extension))
    
    # Canvia el nom del fitxer descomprimit (amb totes les extensions relacionades)
    file.rename(paste0(temporal_file_path, ".bed"), paste0(new_name, ".bed"))
    file.rename(paste0(temporal_file_path, ".fam"), paste0(new_name, ".fam"))
    file.rename(paste0(temporal_file_path, ".bim"), paste0(new_name, ".bim"))
    
    # Retorna el nou camí del fitxer processat
  #  processed_file_path <- new_name
     
    # recover codes from first column
    df2_reactive <- reactive({
      fam <- fread(paste0("temporal/MERGED.fam"), header = F)
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
        selectInput("pop1", "Select Superpopulation: (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
      })
      
      # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
      output$pop2_checkbox <- renderUI({
        req(input$pop1)  # Assegura que POP1 està seleccionat
        df2 <- df2_reactive()
        pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
        selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
        # checkboxGroupInput("pop2", "Select POP2 Groups:", choices = unique(pop2_choices))
      })
    }
    # Final step - Update status and progress bar
    output$status <- renderText("File processed successfully.")
    progress$inc(1)  # Ensure progress is fully completed 
    }
  })
###############################################################################  
  # Function to merge files
  observeEvent(input$run_merge, {
  #  req(input$upload_file)  # Ensure the file is uploaded before proceeding  
    # Obtenir el fitxer processat
    decompressed_file_path <- paste0("decompressed_files/INPUT")
    processed_file_path(decompressed_file_path)
    bfile_path <- processed_file_path()
   
  # Initialize progress bar
  progress <- shiny::Progress$new()
  progress$set(message = "Merging files...", value = 0)
  on.exit(progress$close())  # Ensure progress bar closes on exit
    
    # Function to update progress
    update_progress <- function(amount) {
      progress$inc(amount)
    }
    
    # Set status text to running
    output$status <- renderText("Files merged")
    
    user.bim <- fread(paste0(bfile_path,".bim"))
    # extract from reference those comon to user
    system(paste0("www/plink19 --bfile REFERENCE/Reference --extract ",bfile_path,".bim --aec --make-bed --out temporal/to_merge1"))
    # extract from user those comon to reference-extracted 
    system(paste0("www/plink19 --bfile ",bfile_path," --extract temporal/to_merge1.bim --aec --make-bed --out temporal/to_merge2"))
    
    # merge step 1
    system(paste0("www/plink19 --bfile temporal/to_merge1 --bmerge temporal/to_merge2 --aec --make-bed --out temporal/MERGE1"))
   # Check if ".missnp" file exist
     # Set the file path
    file_path <- "temporal/MERGE1-merge.missnp"
    
    # Check if the file exists
    if (file.exists(file_path)) { # if exist procede as:
      # merge step 1
      system(paste0("www/plink19 --bfile temporal/to_merge1 --bmerge temporal/to_merge2 --aec --make-bed --out temporal/MERGE1"))
      # exclude step
      system(paste0("www/plink19 --bfile temporal/to_merge1 --exclude temporal/MERGE1-merge.missnp --aec --make-bed --out temporal/preMERGE_a"))
      system(paste0("www/plink19 --bfile temporal/to_merge2 --exclude temporal/MERGE1-merge.missnp --aec --make-bed --out temporal/preMERGE_b"))
      system(paste0("www/plink19 --bfile temporal/preMERGE_a --bmerge temporal/preMERGE_b --aec --make-bed --out temporal/MERGED"))
      # merge step 2 
    } else {                      # if not exist procede as:
      # merge step 2    
      system(paste0("www/plink19 --bfile temporal/to_merge1 --bmerge temporal/to_merge2 --aec --make-bed --out temporal/MERGED"))
    }
    
    # generate dataframe to selector populations
    df2_reactive <- reactive({
      fam <- fread(paste0("temporal/MERGED.fam"), header = F)
      fam.1 <- fam[, 1:2]
      colnames(fam.1) = c("IDI", "IDII")
      fam.2 <- fam.1[, 1]
      colnames(fam.2) <- "V1"
      fam2 <- fam.2 %>%
        separate(V1, into = c("POP1", "POP2"), sep = "_")
      df2 <- as.data.frame(cbind(fam2,fam.2))
      return(df2)
    })
    
 #   # UI dinàmic per a la selecció del admixture component
 #   output$sorter <- renderUI({
 #     selectInput("sorter", "Sort admixture plot by 'Ancestry component': ", choices = 1:input$sliderK, multiple = F,selected=1)
 #   }) 
    
    # UI dinàmic per a la selecció de POP1
    output$pop1_select <- renderUI({
      df2 <- df2_reactive()
      selectInput("pop1", "Select Superpopulation: (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
    })
    
    # UI dinàmic per als checkboxes de POP2 basats en el POP1 seleccionat
    output$pop2_checkbox <- renderUI({
      req(input$pop1)  # Assegura que POP1 està seleccionat
      df2 <- df2_reactive()
      pop2_choices <- df2[df2$POP1 == input$pop1, "POP2"]
      selectInput("pop2", "Select a pair of Subpopulations: (FST analysis)", choices = unique(pop2_choices), multiple = TRUE)
     # checkboxGroupInput("pop2", "Select POP2 Groups:", choices = unique(pop2_choices))
    }) 
    shinyjs::enable("run_plink")
    showNotification("Files merged successfully! You can now perform PCA analyses by selecting target populations.", type = "message", duration = 10)
  })
  
  # Executar PLINK quan es fa clic al botó
  observeEvent(input$run_plink, {
    req(processed_file_path())  # Assegurar que el fitxer es va processar
    req(input$pop1)  # Assegurar que hi ha alguna superpoblació seleccionada
    
    # Inicialitzar el bar de progrés
    progress <- shiny::Progress$new()
    progress$set(message = "Running PLINK...", value = 0)
    on.exit(progress$close())  # Assegurar que es tanca el bar
    
    # Funció per actualitzar el bar de progrés
    update_progress <- function(amount) {
      progress$inc(amount)
    }
    
    # Obtenir el fitxer processat
    bfile_path <- "temporal/MERGED"
    
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
    write.table(pop2,"temporal/pop_data.txt", quote = FALSE, row.names = FALSE, col.names = T)
    
    to_keep <- subset(pop2, POP1 %in% input$pop1)
    write.table(to_keep[, 3:4], "temporal/to_keep.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(to_keep[, 1:2], "temporal/to_color.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

    tryCatch({
      # Step 0: Filter the data by selected superpopulation(s) unless "Select all" is checked
      if (input$select_all_pop1=="Select all") {
        # If "Select All" is checked, skip the --keep step
        system(paste("www/plink19 --bfile", bfile_path, "--make-bed --out PCA/step1"))
        update_progress(0.20)  # Update progress to 20%
      } else {
        pop.pre1 <- fread("temporal/to_keep.txt", header = F)
        system(paste("www/plink19 --bfile", bfile_path, "--keep temporal/to_keep.txt --make-bed --out PCA/step1"))
        update_progress(0.20)  # Update progress to 20%
      } 
      # Step 1: Run PLINK with MAF filter
      system(paste("www/plink19 --bfile PCA/step1 --maf 0.05 --make-bed --out PCA/step2"))
      update_progress(0.20)  # Update progress to 40%
      
      # Step 2: Perform linkage disequilibrium (LD) pruning
      system(paste("www/plink19 --bfile PCA/step2 --indep-pairwise 200 25 0.3 --out PCA/prune1"))
      update_progress(0.20)  # Update progress to 60%
      
      # Step 3: Extract pruned SNPs and create a new binary file
      system(paste("www/plink19 --bfile PCA/step2 --extract PCA/prune1.prune.in --make-bed --out PCA/step3"))
      update_progress(0.20)  # Update progress to 80%
      
      # Step 4: Perform PCA
      system(paste("www/plink19 --bfile PCA/step3 --pca 10 header tabs var-wts --out PCA/input_PCA"))
      update_progress(0.20)  # Update progress to 100%
      
      output$status <- renderText("PCA analysis completed successfully.")

    }, error = function(e) {
      showNotification("PLINK analysis failed. Check the logs.", type = "error", duration = 10)
      output$status <- renderText(paste("Error:", e$message))
    })
    
    update_progress(1)  # Finalitzar el progrés
    
    # Read the PCA results from the PLINK output file
    PCA_results <- fread("PCA/input_PCA.eigenvec", header = TRUE)
    
    
    # Render the PCA plot
    output$PCA_plot <- renderPlot({
      to_color <- fread("temporal/to_color.txt", header = F)
      colnames(to_color) <- c("POP1", "POP2")
      
      # Get unique items in the order they first appear
      unique_names <- unique(unique(to_color$POP1))
      # Reorder the Name column based on unique values
      to_color$POP1 <- factor(to_color$POP1, levels = unique_names)
      # Get unique items in the order they first appear
      unique_names2 <- unique(unique(to_color$POP2))
      # Reorder the Name column based on unique values
      to_color$POP2 <- factor(to_color$POP2, levels = unique_names2)
      
      if (input$popcolor == "SUP") {
        PCA <- ggplot(PCA_results, aes(x = PC1, y = PC2, color = factor(to_color$POP1))) +
          geom_point(size = 3, shape=1, stroke = 1.5) +
          labs(title = "PCA Plot", x = "PC1", y = "PC2") +
          guides(color = guide_legend(title = "Superpopulation")) 
      } else if (input$popcolor == "POP") {
        PCA <-  ggplot(PCA_results, aes(x = PC1, y = PC2, color = factor(to_color$POP2), shape = factor(to_color$POP1))) +
          geom_point(size = 3, stroke = 1.5) +
          labs(title = "PCA Plot", x = "PC1", y = "PC2") +
          guides(color = guide_legend(title = "Subpopulation"),
          shape = guide_legend(title = "Superpulation")) 
      }
      plot(PCA)
      ggsave("plots/PCA.png")
    })
    
    shinyjs::enable("run_admix")
    shinyjs::enable("run_fst")
    showNotification("PCA analysis finished. You can now run ADMIXTURE or FST analysis.", type = "message", duration = 10)
    
  })  

  ###########################################
  observeEvent(input$run_admix, {
    req(input$run_plink)  # Ensure PLINK is run first
    
    # Initialize progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Running ADMIXTURE analysis (be patient. This can take a lot of time)...", value = 0)
    on.exit(progress$close())  # Ensure progress bar closes on exit
    
    # Function to update progress
    update_progress <- function(amount) {
      progress$inc(amount)
    }
    
    # Set status text to running
    output$status <- renderText("Running ADMIXTURE analysis...")

    # adjust MAF
    m = input$sliderM
        system(paste("www/plink19 --bfile PCA/step3 --maf", m, "--make-bed --out PCA/input_admx"))  
      update_progress(0.05 / (input$sliderM - 0.05))  # Update progress incrementally based on M

    # Run ADMIXTURE for the selected K value
    for (i in 1:input$sliderK) {
      system(paste("www/admixture --cv PCA/input_admx.bed", i, "| tee log", i, ".out"))
      update_progress(1.0 / (input$sliderK - 1))  # Update progress incrementally based on K
      }
    
      # UI dinàmic per a la selecció del admixture component
      output$sorter <- renderUI({
        selectInput("sorter", "Sort admixture plot by 'Ancestry component': ", choices = 1:input$sliderK, multiple = F,selected=1)
      }) 
      
    # Render the ADMIXTURE plot
    output$ADMIXTURE_plot <- renderPlot({
      # Load data
      runs <- fread(paste0("input_admx.", input$sliderK, ".Q"))
      to_color <- fread("temporal/to_color.txt", header = FALSE)
      colnames(to_color) <- c("POP1", "POP2")

      # Combine population and ancestry data
      admixture_data <- cbind(to_color, runs)
      # Calculate mean values of first (V1) component by POP2
      # create a new column with mean values assigned
      for (i in 1:input$sorter) { 
       ii <- paste0("V",i)
       MEAN_POP2 <- aggregate(admixture_data[[ii]], by = list(admixture_data$POP2), FUN = mean)
      }
      #order by mean value
      MEAN_POP2 = MEAN_POP2[order(MEAN_POP2$x),]
      MEAN_POP2$Group.1
      level_orderX <- c(MEAN_POP2$Group.1)
      admixture_data <- merge(admixture_data, MEAN_POP2, by.x = "POP2", by.y = "Group.1", all.x = TRUE)
      # Determine which population column to use based on input$popcolor
      if (input$popcolor == "SUP") {
        population_column <- "POP1"
      } else {
        population_column <- "POP2"
      }
      
      # Sort admixture data based on the selected population column
      admixture_data <- admixture_data[order(admixture_data[["x"]]), ]
      admixture_data$x <- 1:nrow(admixture_data)  # Add observation column
    #  runs_sorted <- admixture_data[, -(1:2)]  # Remove the first two columns (POP1 and POP2)
      admixture_data <- admixture_data %>%
        relocate(x, .before = V1)

      # Reshape the data for ggplot (long format)
      admixture_long <- pivot_longer(admixture_data, cols = 4:ncol(admixture_data))

      # Calculate the center position for each unique population group based on the observation column
      pop_center_positions <- admixture_data[, .(center = mean(x)), by = list(get(population_column))]
      # Create individual-level admixture plot (stacked bar plot)
      ADM <- ggplot(admixture_long, aes(x = x, y = value, fill = name)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = rainbow(input$sliderK), name = "Ancestry\nComponent") +  # Color scale based on K
        labs(x = "Populations", y = "Ancestry Proportion", fill = "Ancestry Component") +
      #  theme_minimal() +
        scale_x_continuous(breaks = pop_center_positions$center, 
                           labels = pop_center_positions$get)  # Add population labels at center positions
      
      # Plot the graph
      print(ADM)
      
      # Save the plot as PNG
      ggsave("plots/ADM.png", plot = ADM)
    })
      
      output$box_plot <- renderPlot({
        # Load data
        runs <- fread(paste0("input_admx.", input$sliderK, ".Q"))
        # runs <- fread(paste0("input_admx.", input$sliderK, ".Q"))
        to_color <- fread("temporal/to_color.txt", header = FALSE)
        colnames(to_color) <- c("POP1", "POP2")
     #   to_color$observation <- 1:nrow(to_color)  # Add observation column
        
        
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
        #  admixture_data <- admixture_data[order(admixture_data[[population_column]]), ]
      #  runs_sorted <- admixture_data[, -(1:2)]  # Remove the first two columns (POP1 and POP2)
        admixture_data <- admixture_data %>%
          relocate(x, .before = V1)
        # Reshape the data for ggplot (long format)
        admixture_long <- pivot_longer(admixture_data, cols = 4:ncol(admixture_data))
        
        if (input$popcolor == "SUP") {
          p <- ggplot(admixture_long, aes(x = name, y = value, fill = POP1)) + 
            geom_boxplot() +
            labs(x = "Ancestry component", y = "Ancestry Proportion", fill = "Population") +
            theme_minimal()
        } else {

          p <- ggplot(admixture_long, aes(x = name, y = value, fill = factor(POP2, level = level_orderX))) + 
            geom_boxplot() +
            facet_grid(cols = vars(POP1),scales = "free", space = "free") +
            scale_x_discrete(expand = c(0, 0.5)) +
            labs(x = "Ancestry component", y = "Ancestry Proportion", fill = "Population") +
            ggforce::facet_row(vars(POP1), scales = 'free', space = 'free')
          #+
           # theme_minimal()
        }    
      
      print(p)
      # Save the plot as PNG
      ggsave("plots/box_ADM.png", plot = p)
      })
################################################################################      
      output$map_plot <- renderPlot({
        # Load data
        runs <- fread(paste0("input_admx.", input$sliderK, ".Q"))
        to_color <- fread("temporal/to_color.txt", header = FALSE)
        colnames(to_color) <- c("POP1", "POP2")
        coord <- fread(paste0("www/coord.txt"), header = TRUE)
        
        if (input$files == "usr") {
          user_data <-  data.frame(POP1 = input$popID,POP2 = input$popID, Pop_name = "user population", LAT = as.numeric(input$lat),LON = as.numeric(input$lon))
          coord <- rbind(coord,user_data)
        } else {
          coord < coord
        }
        
        pop_counts <-to_color %>% 
          count(POP2)
        
        to_color2 <- to_color %>%
          inner_join(pop_counts, by = "POP2")
        
        # Combine population and ancestry data
        admixture_data <- cbind(to_color, runs)
        admixture_data <- merge(admixture_data, coord, by = c("POP1", "POP2"))
       
        admixture_data <- admixture_data %>%
          relocate(c(Pop_name, LAT, LON), .before = V1)

        cordi <- select(admixture_data, c("POP1", "POP2", "LAT", "LON"))
        
        Npop <- length(unique(admixture_data$POP2))
        qpop <- matrix(NA, ncol = input$sliderK, nrow = Npop)
        
        # Calculate mean ancestry proportions for each population
        for (i in seq_len(Npop)) {  
          pop_subset <- admixture_data$POP2 == unique(admixture_data$POP2)[i]
          qpop[i, ] <- apply(runs[pop_subset, ], 2, mean)
        }
        
        # Calculate mean coordinates for each population
        coord_pop <- matrix(NA, ncol = 2, nrow = Npop) # Assuming LAT and LON as the two columns
        for (i in seq_len(Npop)) {
          pop_subset <- admixture_data$POP2 == unique(admixture_data$POP2)[i]
          coord_pop[i, ] <- apply(cordi[pop_subset, c("LON", "LAT")], 2, mean)
        }
        
        
        # Expand x and y limits by 20%
        xlim <- c(
          min(coord_pop[, 1]) - input$sliderE * (max(coord_pop[, 1]) - min(coord_pop[, 1])), 
          max(coord_pop[, 1]) + input$sliderE * (max(coord_pop[, 1]) - min(coord_pop[, 1]))
        )
        ylim <- c(
          min(coord_pop[, 2]) - input$sliderE * (max(coord_pop[, 2]) - min(coord_pop[, 2])), 
          max(coord_pop[, 2]) + input$sliderE * (max(coord_pop[, 2]) - min(coord_pop[, 2]))
        )

                # Create base map with adjusted limits
        basemap(xlim, ylim)
        plot(coord_pop, xlab = "LON", ylab = "LAT", type = "n", xlim = xlim, ylim = ylim)
        map(add = TRUE, col = "grey90", fill = TRUE)
        
        # Example population sizes for dynamic pie chart radius
        lon_range <- max(coord_pop[, 1]) - min(coord_pop[, 1])  # Longitude range
        lat_range <- max(coord_pop[, 2]) - min(coord_pop[, 2])  # Latitude range
        plot_area <- lon_range * lat_range
        
        R = input$sliderR/200
        radius_scaling_factor <- sqrt(plot_area) * R  # Adjust the divisor based on desired size
        
        pop_size <- runif(Npop, min = 1, max = 5)  # Random example
        
        # Add pie charts with scaled radius
        for (i in seq_len(Npop)) {
          scaled_radius <- pop_size[i] * radius_scaling_factor
          add.pie(z = qpop[i, ], x = coord_pop[i, 1], y = coord_pop[i, 2], labels = "", 
                  radius = scaled_radius, col = rainbow(input$sliderK))
        }
############################## save plot as png ##################################
        lon_range <- max(coord_pop[, 1]) - min(coord_pop[, 1])
        lat_range <- max(coord_pop[, 2]) - min(coord_pop[, 2])
        aspect_ratio <- lon_range / lat_range
        width <- 2400  # Example width in pixels
        height <- width / aspect_ratio
        
        png("plots/map_ADM.png", width = width, height = height, res = 150)
        
        runs <- fread(paste0("input_admx.", input$sliderK, ".Q"))
        to_color <- fread("temporal/to_color.txt", header = FALSE)
        colnames(to_color) <- c("POP1", "POP2")
        coord <- fread(paste0("www/coord.txt"), header = TRUE)
        
        if (input$files == "usr") {
          user_data <-  data.frame(POP1 = input$popID,POP2 = input$popID, Pop_name = "user population", LAT = as.numeric(input$lat),LON = as.numeric(input$lon))
          coord <- rbind(coord,user_data)
        } else {
          coord < coord
        }
        
        pop_counts <-to_color %>% 
          count(POP2)
        
        to_color2 <- to_color %>%
          inner_join(pop_counts, by = "POP2")
        
        # Combine population and ancestry data
        admixture_data <- cbind(to_color, runs)
        admixture_data <- merge(admixture_data, coord, by = c("POP1", "POP2"))
        
        admixture_data <- admixture_data %>%
          relocate(c(Pop_name, LAT, LON), .before = V1)
        
        cordi <- select(admixture_data, c("POP1", "POP2", "LAT", "LON"))
        
        Npop <- length(unique(admixture_data$POP2))
        qpop <- matrix(NA, ncol = input$sliderK, nrow = Npop)
        
        # Calculate mean ancestry proportions for each population
        for (i in seq_len(Npop)) {  
          pop_subset <- admixture_data$POP2 == unique(admixture_data$POP2)[i]
          qpop[i, ] <- apply(runs[pop_subset, ], 2, mean)
        }
        
        # Calculate mean coordinates for each population
        coord_pop <- matrix(NA, ncol = 2, nrow = Npop) # Assuming LAT and LON as the two columns
        for (i in seq_len(Npop)) {
          pop_subset <- admixture_data$POP2 == unique(admixture_data$POP2)[i]
          coord_pop[i, ] <- apply(cordi[pop_subset, c("LON", "LAT")], 2, mean)
        }
        
        # Expand x and y limits by 20%
        xlim <- c(
          min(coord_pop[, 1]) - 0.2 * (max(coord_pop[, 1]) - min(coord_pop[, 1])), 
          max(coord_pop[, 1]) + 0.2 * (max(coord_pop[, 1]) - min(coord_pop[, 1]))
        )
        ylim <- c(
          min(coord_pop[, 2]) - 0.2 * (max(coord_pop[, 2]) - min(coord_pop[, 2])), 
          max(coord_pop[, 2]) + 0.2 * (max(coord_pop[, 2]) - min(coord_pop[, 2]))
        )
        print(admixture_data)
        print(coord_pop)
        # Create base map with adjusted limits
        basemap(xlim, ylim)
        plot(coord_pop, xlab = "LON", ylab = "LAT", type = "n", xlim = xlim, ylim = ylim)
        map(add = TRUE, col = "grey90", fill = TRUE)
        
        # Example population sizes for dynamic pie chart radius
        lon_range <- max(coord_pop[, 1]) - min(coord_pop[, 1])  # Longitude range
        lat_range <- max(coord_pop[, 2]) - min(coord_pop[, 2])  # Latitude range
        plot_area <- lon_range * lat_range
        
        R = input$sliderR/200
        radius_scaling_factor <- sqrt(plot_area) * R  # Adjust the divisor based on desired size
        
        pop_size <- runif(Npop, min = 1, max = 5)  # Random example
        
        # Add pie charts with scaled radius
        for (i in seq_len(Npop)) {
          scaled_radius <- pop_size[i] * radius_scaling_factor
          add.pie(z = qpop[i, ], x = coord_pop[i, 1], y = coord_pop[i, 2], labels = "", 
                  radius = scaled_radius, col = rainbow(input$sliderK))
        }
        
        dev.off()
      })
      
       
  })
  ###########################################
  # Run GCTA Fst test 
  observeEvent(input$run_fst, {
    req(input$run_plink)  # Ensure PLINK is run first
    
    # Step 0: Filter the data by selected superpopulation(s) unless "Select All" is checked
    if (input$select_all_pop2=="Select all") {
      pop2 <- fread("temporal/pop_data.txt", header = TRUE)
      to_fst <- to_fst[, c(3,4,2)]
      write.table(pop2, "temporal/to_fst.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
      system(paste("www/gcta64 --bfile PCA/step3 --fst --sub-popu temporal/to_fst.txt --out FST/FST_output"))
    } else {
      pop2 <- fread("temporal/pop_data.txt", header = TRUE)
      to_fst <- subset(pop2, POP2 %in% input$pop2)
      to_fst <- to_fst[, c(3,4,2)]
      write.table(to_fst, "temporal/to_fst.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
      system(paste("www/gcta64 --bfile PCA/step3 --fst --sub-popu temporal/to_fst.txt --out FST/FST_output"))
    } 
    
    # Generate and render the Manhattan plot
    output$manhattan_plot <- renderPlot({
      pop_names <- fread("temporal/to_fst.txt", header = F)
      names_pop <- unique(pop_names$V3)
      FST_output <- fread("FST/FST_output.fst", header = TRUE)
      FST_output <- FST_output[complete.cases(FST_output[, c("Chr", "bp", "Fst")]), ]
      FST_output <- FST_output[is.finite(FST_output$Fst), ]
      
      # Convert necessary columns to numeric if they aren't already
      FST_output$Chr <- as.numeric(as.character(FST_output$Chr))
      FST_output$bp <- as.numeric(as.character(FST_output$bp))
      FST_output$Fst <- as.numeric(as.character(FST_output$Fst))
      FST_output <- as.data.frame(FST_output)
      # Find the row with the maximum Fst value
      max_fst_row <- FST_output[which.max(FST_output$Fst), ]
      
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
      png("plots/FST.png", width = 2400, height = 1200, res = 150)
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
  ###########################################
#  observeEvent(input$remove, {
#    # Remove previous files in specified directories
#    unlink("temporal/*", recursive = TRUE)
#    unlink("temporal/*.log", recursive = TRUE)
#    
#    unlink("PCA/*", recursive = TRUE)
#    unlink("decompressed_files/*.bim", recursive = TRUE)
#    unlink("decompressed_files/*.fam", recursive = TRUE)
#    unlink("decompressed_files/*.bed", recursive = TRUE)
#    #if exist
#    unlink("1")
#    unlink("2")
#    unlink("3")
#    unlink("4")
#    unlink("5")
#    unlink("6")
#    unlink("7")
#    unlink("8")
#    unlink("9")
#    unlink("10")
#    unlink("log")
#    
#    unlink("ADM/*.Q")
#    unlink("ADM/*.P")
#    unlink("*.Q")
#    unlink("*.P")
#    unlink("log*")
#    unlink(paste0("log", input$sliderK, ".out"))
#    unlink(input$sliderK)
#    unlink("FST/FST_output.log")
#    unlink("FST/FST_output.fst")
#    unlink("plots/PCA.png")
#    unlink("plots/FST.png")
#    unlink("plots/ADM.png")
#    unlink("plots/box_ADM.png")
#    unlink("plots/map_ADM.png")
#    
#   # unlink(temp_folder, recursive = TRUE)
#    
#    # Update status to indicate removal is complete
#    output$status <- renderText("Previous files removed.")
#  })
  ###########################################
  
  observeEvent(input$reset_all, {
    # Reset all inputs
    shinyjs::reset("files")
    shinyjs::reset("upload_file")
    shinyjs::reset("sliderK")
    shinyjs::reset("sliderM")
    shinyjs::reset("sliderR")
    shinyjs::reset("popID")
    shinyjs::reset("lat")
    shinyjs::reset("lon")
    
    # Reset plot outputs
    output$PCA_plot <- renderPlot(NULL)
    output$ADMIXTURE_plot <- renderPlot(NULL)
    output$box_plot <- renderPlot(NULL)
    output$map_plot <- renderPlot(NULL)
    output$manhattan_plot <- renderPlot(NULL)
    
    # Optionally, reset status text or other outputs
    output$status <- renderText("All inputs and plots have been reset.")
    unlink("temporal/*", recursive = TRUE)
    unlink("temporal/*.log", recursive = TRUE)
    
    unlink("PCA/*", recursive = TRUE)
    unlink("decompressed_files/*.bim", recursive = TRUE)
    unlink("decompressed_files/*.fam", recursive = TRUE)
    unlink("decompressed_files/*.bed", recursive = TRUE)
    #if exist
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
    
    unlink("ADM/*.Q")
    unlink("ADM/*.P")
    unlink("*.Q")
    unlink("*.P")
    unlink("log*")
    unlink(paste0("log", input$sliderK, ".out"))
    unlink(input$sliderK)
    unlink("FST/FST_output.log")
    unlink("FST/FST_output.fst")
    unlink("plots/PCA.png")
    unlink("plots/FST.png")
    unlink("plots/ADM.png")
    unlink("plots/box_ADM.png")
    unlink("plots/map_ADM.png")
    
    # Hide progress bar if needed
    output$progress_bar_ui <- renderUI(NULL)
  })
  
  observeEvent(input$reset_session, {
    # Reset only session files
   # shinyjs::reset("files")
   # shinyjs::reset("upload_file")
    shinyjs::reset("sliderK")
    shinyjs::reset("sliderM")
    shinyjs::reset("sliderR")
    shinyjs::reset("popID")
    shinyjs::reset("lat")
    shinyjs::reset("lon")
    
    # Reset plot outputs
    output$PCA_plot <- renderPlot(NULL)
    output$ADMIXTURE_plot <- renderPlot(NULL)
    output$box_plot <- renderPlot(NULL)
    output$map_plot <- renderPlot(NULL)
    output$manhattan_plot <- renderPlot(NULL)
    
    # Optionally, reset status text or other outputs
    output$status <- renderText("All inputs and plots on this session have been reset.")
   # unlink("temporal/*", recursive = TRUE)
   # unlink("temporal/*.log", recursive = TRUE)
    
    unlink("PCA/*", recursive = TRUE)
   # unlink("decompressed_files/*.bim", recursive = TRUE)
   # unlink("decompressed_files/*.fam", recursive = TRUE)
   # unlink("decompressed_files/*.bed", recursive = TRUE)
    #if exist
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
    
    unlink("ADM/*.Q")
    unlink("ADM/*.P")
    unlink("*.Q")
    unlink("*.P")
    unlink("log*")
    unlink(paste0("log", input$sliderK, ".out"))
    unlink(input$sliderK)
    unlink("FST/FST_output.log")
    unlink("FST/FST_output.fst")
    unlink("plots/PCA.png")
    unlink("plots/FST.png")
    unlink("plots/ADM.png")
    unlink("plots/box_ADM.png")
    unlink("plots/map_ADM.png")
    
    # Hide progress bar if needed
    output$progress_bar_ui <- renderUI(NULL)
  })

  ###########################################
  # Download report
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste("Analysis_Report_", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # Generate a report using Rmarkdown
      rmarkdown::render("report_template.Rmd", output_file = file,
                        params = list(
                          pca_plot = "plots/PCA.png",
                          admixture_plot = "plots/ADM.png",
                          box_ADM = "plots/box_ADM.png",
                          map_ADM = "plots/map_ADM.png",
                          manhattan_plot = "plots/FST.png"
                        ))
    }
  )
}

# Reset All button functionality


shinyApp(ui = ui, server = server)

