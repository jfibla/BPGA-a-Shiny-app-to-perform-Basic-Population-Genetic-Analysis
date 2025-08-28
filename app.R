# BPGA — Shiny UI for Basic Population Genetic Analysis
# Revised and structured 08-28-2025
# .libPaths(c("/home/joan_fibla/R/x86_64-pc-linux-gnu-library/4.3", .libPaths()))

# Use a temp folder outside /tmp to avoid disk-full errors
custom_tmp <- "/home/ubuntu/tmp_alt"
dir.create(custom_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = custom_tmp)

options(bitmapType = "cairo")

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
library(shinycssloaders)
library(sf)
library(leaflet)
library(htmlwidgets)
library(leaflet.minicharts)
library(mapview)
library(plotly)
library(Cairo)
library(pandoc)
library(orca)

# Prefer specific functions on name conflicts
conflict_prefer("map", "maps")
conflict_prefer("count", "dplyr")
conflict_prefer("validate", "shiny")

# Increase max upload size to 600 MB
options(shiny.maxRequestSize = 600 * 1024^2)

ui <- fluidPage(
  useShinyjs(),
  
  # Global styles and error suppression in browser console
  tags$head(
    tags$style(HTML("
      body {
        background-image: url('fons.jpg');
        background-size: cover;
        background-repeat: no-repeat;
        background-attachment: fixed;
      }
      .shiny-output-error { visibility: hidden; }
      .shiny-output-error:before { visibility: hidden; }
    "))
  ),
  
  titlePanel("Basic Population Genetic Analysis (BPGA)"),
  
  sidebarLayout(
    sidebarPanel(
      
      # === App info panel shown before starting ===
      conditionalPanel(
        condition = "output.panel_visible === false",
        absolutePanel(
          id = "project_description",
          class = "panel panel-default",
          fixed = FALSE,
          draggable = TRUE,
          top = 100, left = 40, right = 40, bottom = "auto",
          width = "auto", height = "auto",
          style = "
            background-color: #E8E8E8;
            opacity: 0.85;
            z-index: 1000;
            padding: 20px;
            margin: auto;
            border-radius: 4px;
            box-shadow: 0 0 10px rgba(0,0,0,0.2);
            padding-bottom: 2mm;
            padding-top: 1mm;
          ",
          h3("Project description"),
          fluidRow(
            column(
              11,
              p(
                style = "text-align: justify;",
                "This Shiny app provides an interactive platform for visualizing population genetic structure.
                 Users can upload their own datasets or use example files to explore population relationships,
                 genetic diversity, and ancestral components. The app reads PLINK binary files and integrates
                 them with worldwide reference populations from 1000 Genomes (1000G) and the Human Genome
                 Diversity Project (HGDP), enabling PCA, ADMIXTURE, and FST analyses."
              )
            ),
            column(
              11,
              p(
                style = "text-align: justify;",
                "Designed for researchers and students, the app offers a user-friendly interface and
                 produces publication-ready plots."
              )
            ),
            column(
              11,
              strong("Click "),
              actionButton("start_analysis", "Start analysis"),
              strong(" to begin and enjoy.")
            )
          )
        )
      ),
      
      # Always visible controls
      actionButton("reload_button", "Reload session"),
      actionButton("info_00", "ℹ"), # Info about reloading session at start
      verbatimTextOutput("folder_path"),
      
      # === Main controls (visible after starting analysis) ===
      conditionalPanel(
        condition = "output.panel_visible === true",
        
        actionButton("info_0", "ℹ"), # Info about the input file selection
        
        selectInput(
          "files", "Select input file:",
          choices = c(
            "Example single population" = "exp1",
            "Example multi-population"  = "exp2",
            "User file"                 = "usr"
          )
        ),
        
        # ----- User file inputs -----
        conditionalPanel(
          condition = "input.files == 'usr'",
          
          h3("Load user genotype data"),
          fileInput(
            "upload_file",
            "Upload .bed/.bim/.fam or a compressed .zip/.gz",
            accept = c(".zip", ".gz")
          ),
          
          h3("Load input file metadata"),
          tags$hr(style = "border-color: grey;"),
          
          actionButton("info_23", "ℹ"),
          
          radioButtons(
            inputId = "fam_formated",
            label   = "FAM first column format:",
            choices = list(
              "Preformatted as POP1_POP2"           = "pop1pop2",
              "Single population (all samples same)" = "unic",
              "Multiple populations (subpopulations)"= "subpop"
            ),
            selected = character(0)
          ),
          
          # USR code (used when 'unic' or 'subpop')
          conditionalPanel(
            condition = "input.fam_formated == 'unic' || input.fam_formated == 'subpop'",
            textInput("popID", "Assigned code (default 'USR')", "USR")
          ),
          
          # Lat/Lon for 'unic'
          conditionalPanel(
            condition = "input.fam_formated == 'unic'",
            textInput("lat", "Assign latitude"),
            textInput("lon", "Assign longitude")
          ),
          
          # Coordinates file for 'pop1pop2' or 'subpop'
          conditionalPanel(
            condition = "input.fam_formated == 'pop1pop2' || input.fam_formated == 'subpop'",
            actionButton("info_23", "ℹ"),
            fileInput(
              "upload_coord",
              "Upload population coordinates (.txt)",
              accept = ".txt"
            )
          )
        ),
        
        # ----- Processing: build & liftover -----
        conditionalPanel(
          condition = "input.files == 'usr'",
          tags$hr(style = "border-color: grey;"),
          h3("Process input file"),
          actionButton("info_25", "ℹ"),
          radioButtons(
            "build_select", "User genome build:",
            choices  = c("hg18", "hg19", "hg38"),
            selected = "hg38", inline = TRUE
          ),
          conditionalPanel(
            condition = "input.build_select == 'hg18' || input.build_select == 'hg19'",
            actionButton("liftover_to_hg38", "Liftover to hg38")
          )
        ),
        
        tags$hr(style = "border-color: grey;"),
        
        # ----- Process alone vs merge with worldwide -----
        actionButton("info_01", "ℹ"),
        radioButtons(
          "selector", "Select option:",
          choices  = list(
            "Only the loaded population"     = "only",
            "Merge with worldwide populations" = "wwp"
          ),
          selected = character(0)
        ),
        
        # Only-process path
        conditionalPanel(
          condition = "input.selector == 'only'",
          actionButton("decompress_file", "Process input file")
        ),
        
        verbatimTextOutput("wait_ref"),
        
        # Merge path
        conditionalPanel(
          condition = "input.selector == 'wwp'",
          verbatimTextOutput("reference_path"),
          verbatimTextOutput("info_ref"),
          
          # Harmonization only when user file is being merged with reference
          conditionalPanel(
            condition = "input.files == 'usr'",
            h5("Harmonization step"),
            actionButton("info_26", "ℹ"),
            actionButton("harmonize_to_ref", "Harmonize alleles to reference"),
            verbatimTextOutput("harmonization_summary")
          ),
          
          h5("Merging step"),
          actionButton("run_merge", "Merge files")
        ),
        
        verbatimTextOutput("wait_merge"),
        verbatimTextOutput("merged_output"),
        
        # ----- Analyses (enabled after processing or merging) -----
        conditionalPanel(
          condition = "input.decompress_file > 0 || input.run_merge > 0",
          
          h3("Perform population analysis"),
          
          # Single info button (removed duplicated IDs)
          actionButton("info_2", "ℹ"),
          
          # Population selectors
          uiOutput("pop1_select"),                 # Dynamic SelectInput for POP1
          shinyjs::hidden(actionButton("select_all_pop1", "Select all")),
          tags$hr(style = "border-color: grey;"),
          
          # PCA
          h4("PCA analysis"),
          actionButton("info_3", "ℹ"),
          sliderInput(
            "sliderM1", "MAF threshold for pruning:",
            min = 0.05, max = 0.450, value = 0.05, step = 0.05
          ),
          actionButton("run_pca", "Run PCA"),
          conditionalPanel(
            condition = "input.run_pca > 0",
            verbatimTextOutput("pca_log")
          ),
          
          tags$hr(style = "border-color: grey;"),
          
          # FST
          h4("FST analysis: compare two populations"),
          actionButton("info_5", "ℹ"),
          uiOutput("pop2_checkbox"),               # Dynamic Checkbox for POP2
          shinyjs::hidden(actionButton("select_all_pop2", "Remove selected")),
          actionButton("run_fst", "Run FST"),
          conditionalPanel(
            condition = "input.run_fst > 0",
            tags$h5("Top 5 FST values"),
            verbatimTextOutput("top5_fst")
          ),
          
          tags$hr(style = "border-color: grey;"),
          
          # ADMIXTURE
          h4("ADMIXTURE analysis"),
          actionButton("info_4", "ℹ"),
          sliderInput(
            "sliderM2", "MAF threshold:",
            min = 0.05, max = 0.450, value = 0.05, step = 0.05
          ),
          verbatimTextOutput("maf_output"),
          actionButton("info_24", "ℹ"),
          radioButtons(
            "Kselect", "K selection:",
            choices = list(
              "Single K value" = "one",
              "Range of K values" = "range"
            )
          ),
          sliderInput("sliderK", "Set K values:", min = 1, max = 10, value = 1),
          
          actionButton("run_admix", "Run ADMIXTURE"),
          verbatimTextOutput("cvErrors"),
          verbatimTextOutput("logAdmixture"),
          
          tags$hr(style = "border-color: grey;"),
          
          # Plot decorations
          h4("Decorate plots"),
          radioButtons(
            "popcolor", "Color by:",
            choices = list("Superpopulation (POP1)" = "POP1",
                           "Subpopulation (POP2)"  = "POP2")
          ),
          sliderInput(
            "sliderR", "Pie radius on map:",
            min = 20, max = 200, value = 100, step = 10
          ),
          
          uiOutput("sorter"),
          
          tags$hr(style = "border-color: grey;"),
          
          # Report / downloads
          h4("Report"),
          # downloadButton("downloadReport", "Download report"),
          downloadButton("downloadPlotsZip", "Download plots (.zip)")
        )
      )
    ),
    
    # ===================== MAIN PANEL =====================
    mainPanel(
      
      uiOutput("image_ui"),
      uiOutput("image2_ui"),
      
      textOutput("status"),
      uiOutput("progress_bar_ui"),
      
      # PCA plot
      conditionalPanel(
        condition = "input.run_pca > 0",
        h4("PCA plot"),
        withSpinner(plotOutput("PCA_plot"), type = 1, color = "#007BFF")
      ),
      
      # FST Manhattan plot
      conditionalPanel(
        condition = "input.run_fst > 0",
        h4("FST values by chromosome"),
        withSpinner(plotlyOutput("manhattan_plot"), type = 1, color = "#007BFF"),
        uiOutput("dbsnp_link")
      ),
      
      # ADMIXTURE barplot
      conditionalPanel(
        condition = "input.run_admix > 0",
        h4("ADMIXTURE plot"),
        div(
          style = "text-align: center;",
          withSpinner(plotOutput("ADMIXTURE_plot"), type = 1, color = "#007BFF"),
          div(
            id = "loading_message1",
            HTML("Loading data, running ADMIXTURE, and drawing barplot.<br>Please be patient…"),
            style = "color: #007BFF; font-size: 16px; display: block; margin-top: 10px;"
          )
        )
      ),
      
      # ADMIXTURE boxplot
      conditionalPanel(
        condition = "input.run_admix > 0",
        h4("Boxplot of ADMIXTURE components"),
        div(
          style = "text-align: center;",
          withSpinner(plotOutput("box_plot"), type = 1, color = "#007BFF"),
          div(
            id = "loading_message2",
            HTML("Loading data, running ADMIXTURE, and drawing boxplot.<br>Good time for a coffee…"),
            style = "color: #007BFF; font-size: 16px; display: block; margin-top: 10px;"
          )
        )
      ),
      
      # ADMIXTURE map
      conditionalPanel(
        condition = "input.run_admix > 0",
        h4("Geographic distribution of ancestry components"),
        div(
          style = "text-align: center;",
          withSpinner(leafletOutput("map_plot", height = "600px"), type = 1, color = "#007BFF"),
          div(
            id = "loading_message3",
            HTML("Loading data, running ADMIXTURE, and mapping results.<br>Back in a moment…"),
            style = "color: #007BFF; font-size: 16px; display: block; margin-top: 10px;"
          )
        )
      )
    )
  )
)


################################################################################  
########################### SERVER #############################################  
################################################################################  

server <- function(input, output, session) {

  
  # --- Utilities ---------------------------------------------------------------
  source("utils_liftover.R")
  source("utils_harmonize.R", local = FALSE)
  
  # 600 MB upload limit (matches UI)
  options(shiny.maxRequestSize = 600 * 1024^2)
  
  # Use Chromium for webshot2 (e.g., for plot saving)
  options(webshot2.chrome = "/usr/bin/chromium-browser")
  
  # Clean temporary files from previous sessions (if any)
  unlink(file.path(tempdir(), "user_output_*"), recursive = TRUE, force = TRUE)
  
  # --- Reactive state / triggers ----------------------------------------------
  # FST reactivation/tracking at server start
  fst_trigger    <- reactiveVal(0)
  fst_out_files  <- reactiveVal()
  fst_html_path  <- reactiveVal(NULL)
  
  # Temporary folder to store plots and serve them as static resources
  temp_plots <- file.path(tempdir(), "plots")
  dir.create(temp_plots, showWarnings = FALSE, recursive = TRUE)
  addResourcePath("plots", temp_plots)
  ploted_files <- temp_plots
  
  # --- UI loading messages during ADMIXTURE -----------------------------------
  # These three observers show/hide loading messages when ADMIXTURE starts.
  # NOTE: Sys.sleep() here simulates a long process and will block the R session.
  observeEvent(input$run_admix, {
    shinyjs::show("loading_message1")
    Sys.sleep(10)
    shinyjs::hide("loading_message1")
  })
  
  observeEvent(input$run_admix, {
    shinyjs::show("loading_message2")
    Sys.sleep(10)
    shinyjs::hide("loading_message2")
  })
  
  observeEvent(input$run_admix, {
    shinyjs::show("loading_message3")
    Sys.sleep(10)
    shinyjs::hide("loading_message3")
  })
  
  # --- Intro images ------------------------------------------------------------
  # Show infographic before starting the analysis
  output$image_ui <- renderUI({
    if (input$start_analysis == 0) {
      img(src = "infograph.png", height = "700px", width = "800px")
    } else {
      NULL
    }
  })
  
  # Show workflow image after starting and before PCA has run
  output$image2_ui <- renderUI({
    if (input$start_analysis != 0 && input$run_pca == 0) {
      tags$img(
        src = "workflow.png",
        width = "720px",
        height = "405px",
        alt = "BPGA workflow diagram"
      )
    } else {
      NULL
    }
  })
  
  # --- Info modals -------------------------------------------------------------
  observeEvent(input$info_00, {
    showModal(modalDialog(
      title = "Info: Start a new analysis",
      p(style = "text-align: justify;",
        "Reload the session to start a new analysis from scratch."
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_0, {
    showModal(modalDialog(
      title = "Info: Select input files",
      p(style = "text-align: justify;",
        "Select example files (single population or multi-population) or upload your own input file (a compressed ",
        strong(".zip"), " containing ", strong(".bed, .bim, and .fam"), " files). ",
        "Upload size is limited to 600 MB."
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_01, {
    showModal(modalDialog(
      title = "Info: Process files",
      p(style = "text-align: justify;",
        "You can process only the input file or merge it with worldwide reference populations. ",
        "Then click the ", strong("'Decompress file'"), " button to continue."
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_1, {
    showModal(modalDialog(
      title = "Info: Select reference files",
      p(style = "text-align: justify;",
        "Select ", strong("Small"), " or ", strong("Large"), " reference files. Choosing ",
        strong("Large"), " will significantly increase processing time due to a larger (≈27 GB) archive. ",
        strong("Please be patient!")
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_2, {
    showModal(modalDialog(
      title = "Info: Select populations to analyze",
      p(style = "text-align: justify;",
        "From the expanded list above, select ",
        strong("one"), " or ", strong("multiple"), " populations to include in the PCA analysis."
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_3, {
    showModal(modalDialog(
      title = "Info: Select MAF threshold",
      p(style = "text-align: justify;",
        "Select a ", strong("Minor Allele Frequency (MAF)"), " threshold to filter SNPs. ",
        "Higher MAF values reduce processing time but may reduce accuracy.",
        br(),
        em("(Note: Be cautious when selecting MAF = 0.5, as it could remove all SNPs from the analysis.)")
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_4, {
    showModal(modalDialog(
      title = "Info: ADMIXTURE analysis",
      p(style = "text-align: justify;",
        "Choose a ", strong("MAF"), " threshold and specify ", strong("K values"),
        " (number of ancestral populations/clusters) from 1 to 10. ",
        "Higher MAF can speed processing but may compromise accuracy.",
        br(),
        em("Note: Selecting MAF = 0.5 may exclude all SNPs. Processing time also increases with larger K."),
        br(),
        em(span(style = "color: red;",
                "ADMIXTURE may take several minutes (or hours). Enjoy a coffee while you wait!"))
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_24, {
    showModal(modalDialog(
      title = "Info: K selection for ADMIXTURE",
      p(style = "text-align: justify;",
        "You can select a ", strong("single K value"), " or a ",
        strong("range of K values (from 1 up to your chosen K)"), ".",
        br(),
        em("Note: Processing time increases substantially when using a range of K values."),
        br(),
        em(span(style = "color: red;",
                "ADMIXTURE may take several minutes (or hours). Enjoy a coffee while you wait!"))
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_5, {
    showModal(modalDialog(
      title = "Info: Select populations for FST",
      p(style = "text-align: justify;",
        "From the expanded list above, select a pair of populations for the FST analysis."
      ),
      em("Note: You may select more than two populations; in that case, the FST reflects genetic differentiation across all selected populations."),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_21, {
    showModal(modalDialog(
      title = "Info: FAM file formats",
      HTML(
        '<p style="text-align: justify;">
         The FAM file must be formatted as:
         <strong>Preformatted POP1_POP2</strong>
       </p>
       <table style="border-collapse: collapse; width: 100%;" border="1">
         <thead>
           <tr>
             <th>#ID</th><th>IID</th><th>Father</th><th>Mother</th><th>Gender</th><th>Status</th>
           </tr>
         </thead>
         <tbody>
           <tr><td>AFR_AWS</td><td>ASW</td><td>African</td><td>-3.82</td><td>12.93</td></tr>
           <tr><td>AFR_LWK</td><td>LWK</td><td>Luhya_Kenya</td><td>-1.27</td><td>36.61</td></tr>
           <tr><td>AFR_MKK</td><td>MKK</td><td>Maasai</td><td>-1.00</td><td>35.00</td></tr>
           <tr><td>AFR_YRI</td><td>YRI</td><td>Yoruba</td><td>7.40</td><td>3.92</td></tr>
           <tr><td>AMR_MEX</td><td>MEX</td><td>Mexican</td><td>22.50</td><td>-100.00</td></tr>
         </tbody>
       </table>'
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_22, {
    showModal(modalDialog(
      title = "Info: Load coordinate data for user populations",
      HTML(
        '<p style="text-align: justify;">
         Coordinate file should be a data frame with 5 columns:
         <strong>POP1, POP2, Pop_name, LAT, LON</strong>
       </p>
       <table style="border-collapse: collapse; width: 100%;" border="1">
         <thead>
           <tr>
             <th>POP1</th><th>POP2</th><th>Pop_name</th><th>LAT</th><th>LON</th>
           </tr>
         </thead>
         <tbody>
           <tr><td>AFR</td><td>ASW</td><td>African</td><td>-3.82</td><td>12.93</td></tr>
           <tr><td>AFR</td><td>LWK</td><td>Luhya_Kenya</td><td>-1.27</td><td>36.61</td></tr>
           <tr><td>AFR</td><td>MKK</td><td>Maasai</td><td>-1.00</td><td>35.00</td></tr>
           <tr><td>AFR</td><td>YRI</td><td>Yoruba</td><td>7.40</td><td>3.92</td></tr>
           <tr><td>AMR</td><td>MEX</td><td>Mexican</td><td>22.50</td><td>-100.00</td></tr>
         </tbody>
       </table>'
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_23, {
    showModal(modalDialog(
      HTML('
      <div style="font-size: 12px;">
        <p style="text-align: justify;">
          FAM file examples and the assigned coordinate file for each format option,
          following the PLINK binary ".fam" convention.
        </p>

        <!-- Pair 1 -->
        <h4>First column preformatted as <code>POP1_POP2</code></h4>
        <div style="display: flex; gap: 20px;">
          <div style="flex: 1;">
            <table style="border-collapse: collapse; width: 100%;" border="1">
              <thead>
                <tr><th>Family ID</th><th>Individual ID</th><th>Father</th><th>Mother</th><th>Sex</th><th>Phenotype</th></tr>
              </thead>
              <tbody>
                <tr><td>AFR_ASW</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>AFR_ASW</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>EUR_CEU</td><td>NA12341</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
                <tr><td>EUR_CEU</td><td>NA06984</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>EAS_CHB</td><td>NA18532</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
                <tr><td>EAS_CHB</td><td>NA18561</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
              </tbody>
            </table>
          </div>
          <div style="flex: 1;">
            <table style="border-collapse: collapse; width: 100%;" border="1">
              <thead>
                <tr><th>POP1</th><th>POP2</th><th>Pop_name</th><th>LAT</th><th>LON</th></tr>
              </thead>
              <tbody>
                <tr><td>AFR</td><td>ASW</td><td>African</td><td>-3.82</td><td>12.93</td></tr>
                <tr><td>EAS</td><td>CHD</td><td>Chinese</td><td>40.00</td><td>115.00</td></tr>
                <tr><td>EUR</td><td>CEU</td><td>European</td><td>51.30</td><td>11.66</td></tr>
              </tbody>
            </table>
            <p><em>Note:</em> Coordinates are taken from <code>POP2</code> populations; <code>POP1</code> coordinates are estimated as the mean of their <code>POP2</code> subpopulations.</p>
          </div>
        </div>

        <!-- Pair 2 -->
        <h4>First column describes one single population</h4>
        <div style="display: flex; gap: 20px;">
          <div style="flex: 1;">
            <table style="border-collapse: collapse; width: 100%;" border="1">
              <thead>
                <tr><th>Family ID</th><th>Individual ID</th><th>Father</th><th>Mother</th><th>Sex</th><th>Phenotype</th></tr>
              </thead>
              <tbody>
                <tr><td>AFR</td><td>NA19818</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>AFR</td><td>NA20346</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>AFR</td><td>NA19921</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
                <tr><td>AFR</td><td>NA20281</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>AFR</td><td>NA20301</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
                <tr><td>AFR</td><td>NA20294</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
                <tr><td>AFR</td><td>NA20357</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
              </tbody>
            </table>
          </div>
          <div style="flex: 1;">
            <table style="border-collapse: collapse; width: 100%;" border="1">
              <thead><tr><th>Instructions</th></tr></thead>
              <tbody>
                <tr><td>User-assigned population code</td></tr>
                <tr><td>User-provided LAT, LON coordinates</td></tr>
              </tbody>
            </table>
          </div>
        </div>

        <!-- Pair 3 -->
        <h4>First column describes several populations</h4>
        <div style="display: flex; gap: 20px;">
          <div style="flex: 1;">
            <table style="border-collapse: collapse; width: 100%;" border="1">
              <thead>
                <tr><th>Family ID</th><th>Individual ID</th><th>Father</th><th>Mother</th><th>Sex</th><th>Phenotype</th></tr>
              </thead>
              <tbody>
                <tr><td>ASW</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>ASW</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>CEU</td><td>NA12341</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
                <tr><td>CEU</td><td>NA06984</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
                <tr><td>CHB</td><td>NA18532</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
                <tr><td>CHB</td><td>NA18561</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
              </tbody>
            </table>
          </div>
          <div style="flex: 1;">
            <table style="border-collapse: collapse; width: 100%;" border="1">
              <thead>
                <tr><th>POP1</th><th>POP2</th><th>Pop_name</th><th>LAT</th><th>LON</th></tr>
              </thead>
              <tbody>
                <tr><td>ASW</td><td>ASW</td><td>Africa_sample</td><td>9.3</td><td>19.3</td></tr>
                <tr><td>CEU</td><td>CEU</td><td>Europe_sample</td><td>50</td><td>15</td></tr>
                <tr><td>CHB</td><td>CHB</td><td>Asia_sample</td><td>38</td><td>83</td></tr>
              </tbody>
            </table>
            <p><em>Note:</em> Coordinates are taken from <code>POP2</code> populations; <code>POP1</code> coordinates are estimated as the mean of all <code>POP2</code> populations.</p>
          </div>
        </div>
      </div>
    '),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_25, {
    showModal(modalDialog(
      title = "Info: Assign user genome build",
      p(style = "text-align: justify;",
        "Select the genome build of the user genotype data. If hg18 or hg19 is selected, coordinates will be lifted over to hg38 using the appropriate chain file (e.g., 'hg19ToHg38.over.chain')."
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_26, {
    showModal(modalDialog(
      title = "Info: Harmonization step",
      p(style = "text-align: justify;",
        "The harmonization step aligns alleles and SNP identifiers in USER.bim to a reference BIM file: ",
        "it fixes strand orientation, applies complements when needed, removes mismatches and palindromic SNPs, ",
        "and produces a filtered, standardized BIM for downstream analysis."
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  # ----------------------- end info modals -------------------------------------
  # Keep sliderM2's minimum equal to sliderM1
  observeEvent(input$sliderM1, {
    updateSliderInput(session, "sliderM2", min = input$sliderM1)
  })
  
  # Tracks whether the main panel should be visible
  panel_visible <- reactiveVal(FALSE)
  
  # Expose panel_visible to the client (used by conditionalPanels in UI)
  output$panel_visible <- reactive({
    panel_visible()
  })
  outputOptions(output, "panel_visible", suspendWhenHidden = FALSE)
  
  # Folder browser (debug/helper) — shows contents of a given path
  observeEvent(input$show_content, {
    req(input$folder)  # Ensure a folder path is provided
    
    output$folder_content <- renderText({
      folder_path <- input$folder
      
      if (!dir.exists(folder_path)) {
        return(paste("Folder does not exist:", folder_path))
      }
      
      folder_content <- list.files(folder_path, full.names = TRUE)
      if (length(folder_content) == 0) {
        return(paste("The folder is empty:", folder_path))
      }
      
      paste("Content of", folder_path, ":\n", paste(folder_content, collapse = "\n"))
    })
  })
  
  # Path to the processed file (if needed later)
  processed_file_path <- reactiveVal(NULL)
  
  # Disable analysis buttons by default; they will be enabled at the right steps
  # shinyjs::disable("decompress_file")
  shinyjs::disable("run_merge")
  shinyjs::disable("harmonize_to_ref")
  shinyjs::disable("run_pca")
  shinyjs::disable("run_admix")
  shinyjs::disable("run_fst")
  
  ################################################################################
  # Create a per-session working folder with subfolders
  random_folder <- reactiveVal()
  
  # ====== Create temp workspace at analysis start ======
  observeEvent(input$start_analysis, {
    # Toggle visibility of the main analysis panel
    panel_visible(!panel_visible())
    outputOptions(output, "panel_visible", suspendWhenHidden = FALSE)
    
    # Clean custom TMP (if configured) BEFORE creating the new temp folder
    tmp_base <- Sys.getenv("TMPDIR")
    if (nzchar(tmp_base)) {
      unlink(tmp_base, recursive = TRUE, force = TRUE)
      dir.create(tmp_base, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Create a unique session folder under tempdir()
    temp_folder <- tempfile(pattern = "user_output_", tmpdir = tempdir())
    dir.create(temp_folder, recursive = TRUE, showWarnings = FALSE)
    
    # Assign it to your external reactiveVal (defined elsewhere in your server)
    folder_path(temp_folder)  # <- keep as-is; used by the rest of your app
    
    # Create standard subfolders
    subfolders <- c("example_files", "pca", "merged", "reference", "admx", "fst", "plots", "user")
    lapply(subfolders, function(subfolder) {
      path <- file.path(temp_folder, subfolder)
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      cat("✅ Created:", path, "\n")
    })
    
    # Store paths in a single reactive structure
    random_folder(list(
      main_folder = temp_folder,
      subfolders  = subfolders,
      paths       = setNames(lapply(subfolders, function(sub) file.path(temp_folder, sub)), subfolders)
    ))
  })
  
  # ====== Show created folders ======
  output$folder_path <- renderText({
    rf <- random_folder()
    req(rf)
    paste("Temporary user folder:", rf$main_folder)
  })
  
  output$subfolder_paths <- renderText({
    rf <- random_folder()
    req(rf)
    paste("Subfolders created:\n", paste(file.path(rf$main_folder, rf$subfolders), collapse = "\n"))
  })
  
  # ====== Load and prepare reference set when "wwp" (worldwide populations) is selected ======
  # Reference (reactive) from Dropbox
  reference_file <- reactive({
    req(input$selector)
    rf <- random_folder()
    req(rf)
    
    if (input$selector == "wwp") {
      dropbox_url <- "https://www.dropbox.com/scl/fi/73i2iaoybsb7prafp6j4w/Reference_M.zip?rlkey=e5qkyujmuznzszol06s3rtcyu&st=7muwbu5f&dl=1"
      temp_dir <- file.path(rf$main_folder, "reference")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
      
      zip_local_path <- file.path(temp_dir, "Reference_M.zip")
      ref_bed <- file.path(temp_dir, "Reference.bed")
      ref_bim <- file.path(temp_dir, "Reference.bim")
      ref_fam <- file.path(temp_dir, "Reference.fam")
      
      if (!file.exists(zip_local_path)) {
        withProgress(message = "Downloading reference file from repository...", value = 0, {
          download.file(dropbox_url, destfile = zip_local_path, mode = "wb")
          incProgress(1)
        })
      }
      
      if (!all(file.exists(ref_bed, ref_bim, ref_fam))) {
        withProgress(message = "Decompressing reference file...", value = 0, {
          if (Sys.info()[["sysname"]] %in% c("Linux", "Darwin")) {
            system2("unzip", args = c("-o", zip_local_path, "-d", temp_dir), wait = TRUE)
          } else {
            unzip(zip_local_path, exdir = temp_dir)
          }
          incProgress(1, detail = "Finished decompressing.")
        })
      }
      
      shinyjs::enable("run_merge")
      shinyjs::enable("harmonize_to_ref")
      shiny::showNotification("Reference file ready.", duration = 5, type = "message")
      
      return(file.path(temp_dir, "Reference"))
    }
    
    return(NULL)
  })
  
  # ====== Reference file info/status ======
  output$reference_path <- renderText({
    ref_file <- reference_file()
    if (is.null(ref_file)) {
      return("Please select the reference option to proceed.")
    }
    paste("Reference file at:", ref_file)
  })
  
  output$wait_ref <- renderText({
    req(input$selector)
    if (input$selector == "wwp") {
      "Please wait while the reference file is decompressed…"
    } else {
      NULL
    }
  })
  
  merge_status <- reactiveVal("")
  output$wait_merge <- renderText({ merge_status() })
  
  harmonize_status <- reactiveVal("")
  output$wait_harmonize <- renderText({ harmonize_status() })
  
  output$info_ref <- renderText({
    req(input$selector)
    if (input$selector == "wwp") {
      "Reference file:\nTotal genotyping rate is 0.997.\n3,000,000 variants and 3,258 individuals pass filters and QC."
    } else {
      NULL
    }
  })
  
  ########################################################################################
  # ====================== Load, decompress, and normalize user file ====================
  ########################################################################################
  
  observeEvent(input$upload_file, {
    req(input$files == "usr")
    rf <- random_folder()
    req(rf)
    folder_path <- rf$main_folder
    
    file_ext <- tools::file_ext(input$upload_file$datapath)
    
    if (file_ext == "zip") {
      # Decompress ZIP into session's 'user' subfolder
      withProgress(message = "Decompressing user file...", value = 0, {
        decompress_dir <- file.path(folder_path, "user")
        if (Sys.info()[["sysname"]] %in% c("Linux", "Darwin")) {
          system2("unzip", args = c("-o", input$upload_file$datapath, "-d", decompress_dir), wait = TRUE)
        } else {
          unzip(input$upload_file$datapath, exdir = decompress_dir)
        }
        incProgress(1, detail = "Finished decompressing.")
        
        # Remove macOS metadata folder if present
        macosx_path <- file.path(decompress_dir, "__MACOSX")
        if (dir.exists(macosx_path)) {
          unlink(macosx_path, recursive = TRUE, force = TRUE)
        }
      })
    } else if (file_ext == "gz") {
      # Decompress .gz in place (note: this occurs in Shiny's temp upload dir)
      system(paste("gunzip", input$upload_file$datapath))
    } else {
      stop("Unsupported file type.")
    }
    
    # Find the .bim inside the session 'user' folder and derive base path
    bim_path <- list.files(file.path(folder_path, "user"), pattern = "\\.bim$", full.names = TRUE)
    stopifnot(length(bim_path) > 0)
    base_name <- tools::file_path_sans_ext(basename(bim_path))
    user_file_path <- file.path(folder_path, "user", base_name)
    
    # Keep only autosomal SNPs (CHR 1..22)
    bim <- fread(paste0(user_file_path, ".bim"), header = FALSE)
    bim[[1]] <- as.numeric(bim[[1]])
    bim_auto <- bim[bim[[1]] %in% 1:22, ]
    
    # Write list of SNP IDs to keep
    snps_to_keep_file <- file.path(folder_path, "keep_snps.txt")
    fwrite(data.frame(SNP = bim_auto[[2]]), snps_to_keep_file, col.names = FALSE)
    
    # Filter to autosomes with PLINK and write a new BED/BIM/FAM
    filtered_output_path <- file.path(folder_path, "user", "filtered_tmp")
    system(paste0(
      "www/plink19 --bfile ", user_file_path,
      " --extract ", snps_to_keep_file,
      " --allow-extra-chr --allow-no-sex --make-bed --out ", filtered_output_path
    ))
    
    # Ensure PLINK created outputs
    if (!file.exists(paste0(filtered_output_path, ".bed"))) {
      stop("❌ Error: PLINK did not generate filtered files. Is keep_snps.txt empty?")
    }
    
    # Normalize allele codes: convert 1/2/3/4 -> A/C/G/T (retain 0)
    library(data.table)
    
    in_bim  <- paste0(filtered_output_path, ".bim")
    out_bim <- paste0(filtered_output_path, "_alleles_fix.bim")  # (kept for reference)
    
    bim <- fread(in_bim, header = FALSE)
    setnames(bim, c("CHR", "SNP", "REC", "POS", "A1", "A2"))
    
    # Normalize types (avoid factors/spaces); upper-case alleles
    bim[, A1 := toupper(trimws(as.character(A1)))]
    bim[, A2 := toupper(trimws(as.character(A2)))]
    
    # Lookup table applied only to 0/1/2/3/4
    lut <- c("0" = "0", "1" = "A", "2" = "C", "3" = "G", "4" = "T")
    bim[A1 %chin% names(lut), A1 := lut[A1]]
    bim[A2 %chin% names(lut), A2 := lut[A2]]
    
    # Validation
    stopifnot(all(bim$A1 %in% c("A", "C", "G", "T", "0")))
    stopifnot(all(bim$A2 %in% c("A", "C", "G", "T", "0")))
    
    # Overwrite BIM with normalized alleles (keeps same basename)
    fwrite(bim, paste0(filtered_output_path, ".bim"), sep = "\t", col.names = FALSE)
    
    # Store path for downstream steps
    processed_file_path(filtered_output_path)
    
    output$status <- renderText("✅  User file loaded")
  })
  
  #############################################################################
  # =========================== Load user metadata ============================
  
  observeEvent(input$upload_coord, {
    req(input$files == "usr")
    rf <- random_folder()
    req(rf)
    folder_path <- rf$main_folder
    
    # Ensure 'user' subfolder exists
    user_dir <- file.path(folder_path, "user")
    dir.create(user_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Read uploaded coordinate file
    coord_df <- fread(input$upload_coord$datapath, header = TRUE)
    
    # Add "*" to first two columns without altering header
    coord_df <- coord_df %>%
      mutate(across(1:2, ~ paste0(.x, "*")))
    colnames(coord_df) <- c("POP1", "POP2", "Pop_name", "LAT", "LON")
    
    # Save as 'user_coord.txt' inside the 'user' folder
    coord_path <- file.path(user_dir, "user_coord.txt")
    fwrite(coord_df, coord_path, sep = "\t", col.names = TRUE)
    
    output$status <- renderText("✅ Loaded user metadata")
  })
  
  #################################################################################
  # ============================== FORMAT USER.fam ================================
  #################################################################################
  observeEvent(input$fam_formated, {
    req(processed_file_path())
    req(input$fam_formated)
    
    fam_path  <- paste0(processed_file_path(), ".fam")
    user.fam  <- fread(fam_path, header = FALSE)
    colnames(user.fam) <- c("V1","V2","V3","V4","V5","V6")
    
    if (input$fam_formated == "unic") {
      # Single population: set FID to "<USR*>_<USR*>"
      user.id      <- paste0(input$popID, "*_", input$popID, "*")
      user.fam$V1  <- user.id
    } else if (input$fam_formated == "subpop") {
      # Multiple populations: prefix existing FID with "<USR*>_"
      user.fam$V1 <- paste0(input$popID, "*_", user.fam$V1, "*")
    } else if (input$fam_formated == "pop1pop2") {
      # Preformatted POP1_POP2: inject asterisks around both parts
      user.fam$V1 <- sub("(.*)_(.*)", "\\1*_\\2*", user.fam$V1)
    }
    
    # Save modified .fam, then rename filtered files to USER.*
    write.table(user.fam, fam_path, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    file.rename(paste0(processed_file_path(), ".bed"), file.path(dirname(processed_file_path()), "USER.bed"))
    file.rename(paste0(processed_file_path(), ".bim"), file.path(dirname(processed_file_path()), "USER.bim"))
    file.rename(paste0(processed_file_path(), ".fam"), file.path(dirname(processed_file_path()), "USER.fam"))
    
    # Build POP1 / POP2 selectors from USER.fam
    fam   <- fread(file.path(dirname(processed_file_path()), "USER.fam"), header = FALSE)
    fam.1 <- fam[, 1:2]
    colnames(fam.1) <- c("FID","IID")
    
    fam.2 <- fam.1[, 1]           # first column (FID) with POP1*_POP2*
    colnames(fam.2) <- "V1"
    
    fam2 <- fam.2 %>%              # split into POP1 / POP2
      separate(V1, into = c("POP1","POP2"), sep = "_", extra = "merge")
    
    df2 <- as.data.frame(cbind(fam2, fam.2))
    
    output$pop1_select <- renderUI({
      selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)",
                  choices = unique(df2$POP1), multiple = TRUE)
    })
    
    observe({
      req(df2)
      updateSelectInput(session, "pop1", selected = unique(df2$POP1)[1])
    })
    
    output$pop2_checkbox <- renderUI({
      req(input$pop1)
      pop2_choices <- df2[df2$POP1 %in% input$pop1, "POP2"]
      selectInput("pop2", "Select a pair of populations:",
                  choices = unique(pop2_choices), multiple = TRUE)
    })
    
    output$status <- renderText("✅ FAM file formatted and saved successfully.")
  })
  
  #################################################################################
  # =========================== LIFTOVER to hg38 (if needed) ======================
  #################################################################################
  observeEvent(input$liftover_to_hg38, {
    req(processed_file_path())
    showNotification("Liftover running...", type = "message")
    
    if (input$build_select %in% c("hg18","hg19")) {
      showNotification(paste0("🔁 Liftover: ", input$build_select, " → hg38"), type = "message")
      
      user_base <- file.path(dirname(processed_file_path()), "USER")
      bim_file  <- paste0(user_base, ".bim")
      
      validate(need(file.exists(bim_file), "❌ USER.bim not found. Load data first."))
      
      bim <- data.table::fread(bim_file, header = FALSE)
      data.table::setnames(bim, c("CHR","SNP","REC","POS","A1","A2"))
      
      withProgress(message = "Converting USER to hg38...", value = 0, {
        temp_lifted <- file.path(dirname(user_base), "lifted_USER")
        
        # Pass source genome build to liftover utility
        run_liftover(bim, temp_lifted, from_build = input$build_select)
        incProgress(0.30, detail = "✔️ LiftOver running")
        
        lifted_bim <- data.table::fread(paste0(temp_lifted, ".bim"), header = FALSE)
        validate(need(nrow(lifted_bim) > 0, "❌ LiftOver produced no mapped SNPs."))
        data.table::setnames(lifted_bim, c("CHR","SNP","REC","POS","A1","A2"))
        
        # SNP list to extract
        snp_extract <- file.path(dirname(user_base), "extract_snps.txt")
        data.table::fwrite(lifted_bim[, .(SNP)], snp_extract, col.names = FALSE)
        
        # Subset USER to lifted SNPs with PLINK
        subset_base <- file.path(dirname(user_base), "USER_subset")
        system(paste(
          "www/plink19 --bfile", user_base,
          "--extract", snp_extract,
          "--make-bed --out", subset_base,
          "--allow-extra-chr --allow-no-sex"
        ))
        incProgress(0.60, detail = "✔️ Subset created")
        
        # Replace USER.* with lifted coordinates and subset genotypes
        file.copy(paste0(subset_base, ".bed"), paste0(user_base, ".bed"), overwrite = TRUE)
        file.copy(paste0(temp_lifted, ".bim"), paste0(user_base, ".bim"), overwrite = TRUE)
        file.copy(paste0(subset_base, ".fam"), paste0(user_base, ".fam"), overwrite = TRUE)
        
        incProgress(1.00, detail = "✔️ USER.bim at hg38")
      })
      
      showNotification("✅ LiftOver completed", type = "message")
      output$status <- renderText("✅  LiftOver completed")
      
    } else {
      showNotification("✅ Already on hg38. LiftOver not required.", type = "message")
    }
  })
  
  ###################################################################################################################
  # ========================================= HARMONIZE ============================================================
  ###################################################################################################################
  # Reactive to store harmonization summary
  harmonization_result <- reactiveVal(NULL)
  
  observeEvent(input$harmonize_to_ref, {
    on.exit({
      # Re-enable UI, etc., if you disable anything above.
    }, add = TRUE)
    
    harmonize_status("⏳ Harmonizing alleles to reference...")
    req(processed_file_path())
    req(reference_file())
    
    user_base <- file.path(dirname(processed_file_path()), "USER")
    user_bim  <- paste0(user_base, ".bim")
    user_bed  <- paste0(user_base, ".bed")
    user_fam  <- paste0(user_base, ".fam")
    
    temp_harmo <- file.path(dirname(user_base), "harmo_USER")
    harmo_base <- file.path(dirname(user_base), "USER_harmo")
    
    # Reference from reactive
    ref_base <- reference_file()
    ref_bim  <- paste0(ref_base, ".bim")
    
    withProgress(message = "Harmonizing...", value = 0, {
      setProgress(0.1, detail = "📄 Loading reference and user files...")
      
      result <- tryCatch({
        harmonize_bim_to_reference(
          user_bim_path = user_bim,
          ref_bim_path  = ref_bim,
          out_bim_path  = temp_harmo
        )
      }, error = function(e) {
        showNotification(paste("❌ Error:", e$message), type = "error", duration = 10)
        harmonize_status("❌ Harmonization failed.")
        return(NULL)
      })
      req(result)
      
      print("✅ Harmonization result:")
      print(result)
      setProgress(0.5, detail = sprintf("✅ Harmonized: %s SNPs retained", result$n_kept))
      
      # Read harmonized BIM from temp_harmo
      harmo_bim <- data.table::fread(paste0(temp_harmo, ".bim"), header = FALSE)
      data.table::setnames(harmo_bim, c("CHR","SNP","CM","POS","A1","A2"))
      
      # SNP list to extract
      snp_harmo <- file.path(dirname(user_base), "harmo_snps.txt")
      data.table::fwrite(harmo_bim[, .(SNP)], snp_harmo, col.names = FALSE)
      
      # Rebuild binary set with PLINK using only harmonized SNPs
      setProgress(0.65, detail = "🔄 Rebuilding binary files with PLINK...")
      plink_bin  <- "www/plink19"  # or simply "plink" if on PATH
      plink_args <- c(
        "--bfile", shQuote(user_base),
        "--extract", shQuote(snp_harmo),
        "--make-bed",
        "--out", shQuote(harmo_base),
        "--allow-extra-chr",
        "--allow-no-sex"
      )
      plink_out <- system2(plink_bin, args = plink_args, stdout = TRUE, stderr = TRUE)
      
      if (!all(file.exists(paste0(harmo_base, c(".bed",".bim",".fam"))))) {
        showNotification("❌ PLINK failed. Check the console log.", type = "error", duration = 10)
        print(plink_out)
        harmonize_status("❌ Harmonization failed during PLINK step.")
        return(NULL)
      }
      setProgress(0.8, detail = "✔️ Subset created")
      
      # Ensure SNP order matches before copying harmonized BIM onto USER.*
      hb <- data.table::fread(paste0(harmo_base, ".bim"), header = FALSE)
      if (!identical(hb[[2]], harmo_bim$SNP)) {
        # Reorder temp_harmo.bim to match harmo_base.bim SNP order
        setkey(harmo_bim, SNP)
        harmo_bim <- harmo_bim[J(hb[[2]])]
        data.table::fwrite(harmo_bim, paste0(temp_harmo, ".bim"), col.names = FALSE, sep = "\t")
      }
      
      # Copy back to USER.*
      ok <- file.copy(paste0(harmo_base, ".bed"), paste0(user_base, ".bed"), overwrite = TRUE) &
        file.copy(paste0(temp_harmo, ".bim"),  paste0(user_base, ".bim"), overwrite = TRUE) &
        file.copy(paste0(harmo_base, ".fam"), paste0(user_base, ".fam"), overwrite = TRUE)
      validate(need(ok, "❌ Could not update USER.*"))
      
      setProgress(1.0, detail = "✅ USER.* files updated")
      harmonization_result(result)
    })
    
    harmonize_status("✅ Harmonization complete.")
  })
  
  output$harmonization_summary <- renderText({
    res <- harmonization_result()
    if (is.null(res)) return(NULL)
    
    reason_text <- ""
    if (!is.null(res$removed)) {
      reasons <- table(res$removed$remove_reason)
      reason_text <- paste(names(reasons), reasons, sep = ": ", collapse = "\n")
    }
    
    paste0(
      "📊 Harmonization results\n",
      "---------------------------\n",
      "SNPs input:       ", res$n_input, "\n",
      "SNPs harmonized:  ", res$n_kept, "\n",
      "SNPs discarded:   ", res$n_input - res$n_kept, "\n\n",
      if (reason_text != "") paste0("🔍 Reasons for discard:\n", reason_text) else ""
    )
  })
 
  ########################################################################################################
  # =================== Decompress example files and rename them as INPUT (exp1 / exp2) =================
  ########################################################################################################
  observeEvent(input$decompress_file, {
    rf <- random_folder()
    req(rf)
    folder_path <- rf$main_folder
    
    # Allow PCA once example data is prepared
    shinyjs::enable("run_pca")
    
    # Only act when an example dataset is selected
    if (input$files == "exp1" || input$files == "exp2") {
      
      # Choose example archive and base name
      if (input$files == "exp1") {
        ex_file_path <- "example/example_single.zip"  # single-population example
        base_name    <- "example"
      } else {
        ex_file_path <- "example/example_multi.zip"   # multi-population example
        base_name    <- "example"
      }
      
      # Decompress into session's example_files folder
      withProgress(message = "Decompressing example file...", value = 0, {
        ex_dir <- file.path(folder_path, "example_files")
        if (Sys.info()[["sysname"]] %in% c("Linux", "Darwin")) {
          system2("unzip", args = c("-o", ex_file_path, "-d", ex_dir), wait = TRUE)
        } else {
          unzip(ex_file_path, exdir = ex_dir)
        }
        incProgress(1, detail = "Finished decompressing.")
        
        # Remove macOS metadata folder if present
        macosx_path <- file.path(ex_dir, "__MACOSX")
        if (dir.exists(macosx_path)) {
          unlink(macosx_path, recursive = TRUE, force = TRUE)
        }
      })
      
      # Path to decompressed trio (without extension)
      decompressed_file_path <- file.path(folder_path, "example_files", base_name)
      
      # Register processed path (reactive) for downstream steps
      processed_file_path(decompressed_file_path)
      
      # Rename decompressed PLINK trio to INPUT.bed/.bim/.fam (keep extensions explicit)
      new_name <- file.path(folder_path, "example_files", "INPUT")
      file.rename(paste0(decompressed_file_path, ".bed"), paste0(new_name, ".bed"))
      file.rename(paste0(decompressed_file_path, ".fam"), paste0(new_name, ".fam"))
      file.rename(paste0(decompressed_file_path, ".bim"), paste0(new_name, ".bim"))
      
      # Update reactive to the new INPUT base path (important!)
      processed_file_path(new_name)
      
      # Build small helper reactive to parse POP1/POP2 from INPUT.fam
      df2_reactive <- reactive({
        fam   <- fread(file.path(folder_path, "example_files", "INPUT.fam"), header = FALSE)
        fam.1 <- fam[, 1:2]
        colnames(fam.1) <- c("FID", "IID")
        fam.2 <- fam.1[, 1]
        colnames(fam.2) <- "V1"
        fam2 <- fam.2 %>%
          tidyr::separate(V1, into = c("POP1", "POP2"), sep = "_")
        as.data.frame(cbind(fam2, fam.2))
      })
      
      # Dynamic UI: POP1 selector
      output$pop1_select <- renderUI({
        df2 <- df2_reactive()
        selectInput(
          "pop1",
          "Select Population(s): (PCA & Admixture analysis)",
          choices = unique(df2$POP1),
          multiple = TRUE
        )
      })
      
      # Dynamic UI: POP2 selector depends on selected POP1
      output$pop2_checkbox <- renderUI({
        req(input$pop1)
        df2 <- df2_reactive()
        pop2_choices <- df2[df2$POP1 %in% input$pop1, "POP2"]
        selectInput(
          "pop2",
          "Select a pair of populations:",
          choices = unique(pop2_choices),
          multiple = TRUE
        )
      })
      
      # Status message
      output$status <- renderText("✅  File processed successfully.")
    }
  })
  
  #####################################################################################
  # ################################### MERGE STEP ####################################
  # ########################## when 'wwp' option is selected ##########################
  observeEvent(input$run_merge, {
    merge_status("⏳ Merging reference and user populations...")
    
    rf <- random_folder()
    req(rf)
    
    folder_path      <- rf$main_folder
    path_decompressed <- rf$paths$example_files  # example files
    path_user         <- rf$paths$user           # user files
    path_merged       <- rf$paths$merged
    ref_file_path     <- file.path(folder_path, "reference")  # ✅ ensure Reference path
    
    if (input$files %in% c("exp1", "exp2")) {
      bfile_prefix <- file.path(path_decompressed, "INPUT")
      
      # 🧠 If INPUT.* does not exist yet, auto-decompress the corresponding example
      if (!file.exists(paste0(bfile_prefix, ".bed"))) {
        ex_file_path <- if (input$files == "exp1") {
          "example/example_single.zip"
        } else {
          "example/example_multi.zip"
        }
        
        withProgress(message = "Decompressing example file...", value = 0, {
          if (Sys.info()[["sysname"]] %in% c("Linux", "Darwin")) {
            system2("unzip", args = c("-o", ex_file_path, "-d", path_decompressed), wait = TRUE)
          } else {
            unzip(ex_file_path, exdir = path_decompressed)
          }
          
          base_name <- "example"
          file.rename(file.path(path_decompressed, paste0(base_name, ".bed")), file.path(path_decompressed, "INPUT.bed"))
          file.rename(file.path(path_decompressed, paste0(base_name, ".bim")), file.path(path_decompressed, "INPUT.bim"))
          file.rename(file.path(path_decompressed, paste0(base_name, ".fam")), file.path(path_decompressed, "INPUT.fam"))
        })
      }
      
    } else {
      bfile_prefix <- file.path(path_user, "USER")
      
      # 🧠 If USER.* does not exist, the user hasn't uploaded/processed their file yet
      if (!file.exists(paste0(bfile_prefix, ".bed"))) {
        showNotification("❌ You must upload and process the user file before merging.", type = "error")
        return(NULL)
      }
    }
    
    # Initialize a manual progress bar (note: not strictly needed if you prefer withProgress)
    progress <- shiny::Progress$new()
    progress$set(message = "Merging files...", value = 0)
    on.exit({
      progress$close()  # always close progress bar
    }, add = TRUE)
    
    # Helper (currently unused, kept to preserve original structure)
    update_progress <- function(amount) {
      progress$inc(amount)
    }
    
    # Status text (kept in original position)
    output$status <- renderText("✅  Files merged")
    
    # 1) Extract from reference those SNPs present in user/example BIM
    system(paste0(
      "www/plink19 --bfile ", ref_file_path, "/Reference",
      " --extract ", bfile_prefix, ".bim",
      " --aec --allow-no-sex --make-bed",
      " --out ", path_merged, "/to_merge1"
    ))
    
    # 2) Extract from user/example those SNPs present in the reference-extracted set
    system(paste0(
      "www/plink19 --bfile ", bfile_prefix,
      " --extract ", path_merged, "/to_merge1.bim",
      " --aec --allow-no-sex --make-bed",
      " --out ", path_merged, "/to_merge2"
    ))
    
    # 3) First merge attempt
    system(paste0(
      "www/plink19 --bfile ", path_merged, "/to_merge1",
      " --bmerge ", path_merged, "/to_merge2",
      " --aec --allow-no-sex --make-bed",
      " --out ", path_merged, "/MERGE1"
    ))
    
    # If a ".missnp" file exists, perform the exclude-and-merge fallback
    miss_file_path <- paste0(path_merged, "/MERGE1-merge.missnp")
    
    if (file.exists(miss_file_path)) {
      # Repeat merge, then exclude problematic SNPs on both sides and merge again
      system(paste0(
        "www/plink19 --bfile  ", path_merged, "/to_merge1",
        " --bmerge ", path_merged, "/to_merge2",
        " --aec --allow-no-sex --make-bed",
        " --out ", path_merged, "/MERGE1"
      ))
      system(paste0(
        "www/plink19 --bfile  ", path_merged, "/to_merge1",
        " --exclude  ", path_merged, "/MERGE1-merge.missnp",
        " --aec --allow-no-sex --make-bed",
        " --out ", path_merged, "/preMERGE_a"
      ))
      system(paste0(
        "www/plink19 --bfile  ", path_merged, "/to_merge2",
        " --exclude  ", path_merged, "/MERGE1-merge.missnp",
        " --aec --allow-no-sex --make-bed",
        " --out ", path_merged, "/preMERGE_b"
      ))
      system(paste0(
        "www/plink19 --bfile  ", path_merged, "/preMERGE_a",
        " --bmerge  ", path_merged, "/preMERGE_b",
        " --aec --allow-no-sex --make-bed",
        " --out ", path_merged, "/MERGED"
      ))
    } else {
      # Straight merge if there were no problematic SNPs
      system(paste0(
        "www/plink19 --bfile  ", path_merged, "/to_merge1",
        " --bmerge  ", path_merged, "/to_merge2",
        " --aec --allow-no-sex --make-bed",
        " --out ", path_merged, "/MERGED"
      ))
    }
    
    # Activate processed_file_path pointing to the merged dataset
    processed_file_path(file.path(path_merged, "MERGED"))
    
    # Parse POP1/POP2 codes from the first column (FID) of MERGED.fam
    df2_reactive <- reactive({
      fam   <- fread(paste0(path_merged, "/MERGED.fam"), header = FALSE)
      fam.1 <- fam[, 1:2]
      colnames(fam.1) <- c("FID", "IID")
      fam.2 <- fam.1[, 1]
      colnames(fam.2) <- "V1"
      fam2 <- fam.2 %>%
        separate(V1, into = c("POP1", "POP2"), sep = "_")
      as.data.frame(cbind(fam2, fam.2))
    })
    
    # Dynamic UI for POP1 selection
    output$pop1_select <- renderUI({
      df2 <- df2_reactive()
      selectInput(
        "pop1",
        "Select Population(s): (PCA & Admixture analysis)",
        choices = unique(df2$POP1),
        multiple = TRUE
      )
    })
    
    # Dynamic UI for POP2 selection (depends on chosen POP1)
    output$pop2_checkbox <- renderUI({
      req(input$pop1)  # ensure POP1 is selected
      df2 <- df2_reactive()
      pop2_choices <- df2[df2$POP1 %in% input$pop1, "POP2"]
      selectInput(
        "pop2",
        "Select a pair of populations:",
        choices = unique(pop2_choices),
        multiple = TRUE
      )
    })
    
    merge_status("✅ Merging complete.")
    shinyjs::enable("run_pca")
  })
  
  ###########################################################################################
  # ====================================== PCA STEP ========================================
  ###########################################################################################
  
  # Global reactives
  pca_files <- reactiveVal(NULL)
  PCA_data  <- reactiveVal(NULL)
  
  observeEvent(input$run_pca, {
    req(processed_file_path())
    req(input$pop1)
    
    cat("📁 folder_path():", folder_path(), "\n")
    cat("📄 processed_file_path():", processed_file_path(), "\n")
    cat("🌍 input$pop1:", input$pop1, "\n")
    
    print("🔄 Updating `PCA_data()` with new values...")
    
    # Progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Running PCA with PLINK...", value = 0)
    on.exit({ progress$close() }, add = TRUE)
    update_progress <- function(amount) progress$inc(amount)
    
    # Session paths
    rf          <- random_folder()
    base_folder <- rf$main_folder
    pca_dir     <- file.path(base_folder, "pca")
    pca_files(pca_dir)  # store PCA folder path for later use
    
    path_example <- rf$paths$example_files
    path_user    <- rf$paths$user
    path_merged  <- rf$paths$merged
    
    # Choose dataset to analyze
    if (input$files %in% c("exp1", "exp2") && input$selector != "wwp") {
      bfile_path     <- file.path(path_example, "INPUT")
      required_files <- file.path(path_example, paste0("INPUT", c(".bed", ".bim", ".fam")))
    } else if (input$files == "usr" && input$selector != "wwp") {
      bfile_path     <- file.path(path_user, "USER")
      required_files <- file.path(path_user, paste0("USER", c(".bed", ".bim", ".fam")))
    } else {
      bfile_path     <- file.path(path_merged, "MERGED")
      required_files <- file.path(path_merged, paste0("MERGED", c(".bed", ".bim", ".fam")))
    }
    
    if (!all(file.exists(required_files))) {
      showNotification("Missing required PLINK files (.bed, .bim, .fam).", type = "error", duration = 10)
      return(NULL)
    }
    
    # Build POP1/POP2 table from FAM
    pop.pre1    <- fread(paste0(bfile_path, ".fam"), header = FALSE)
    pop.pre1.1  <- pop.pre1[, 1:2]
    colnames(pop.pre1.1) <- c("IDI", "IDII")
    pop.pre2    <- pop.pre1[, 1]
    colnames(pop.pre2) <- "V1"
    
    pop2 <- pop.pre2 %>%
      separate(V1, into = c("POP1", "POP2"), sep = "_") %>%
      cbind(pop.pre1.1)
    
    cat("Current state ➤ selector:", input$selector, ", fam_formated:", input$fam_formated, "\n")
    
    # Persist helper files for PLINK filtering and coloring
    write.table(pop2, file.path(pca_dir, "pop_data.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    to_keep <- subset(pop2, POP1 %in% input$pop1)
    write.table(to_keep[, 3:4], file.path(pca_dir, "to_keep.txt"),  quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(to_keep[, 1:2], file.path(pca_dir, "to_color.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Run PCA pipeline with PLINK
    tryCatch({
      # NOTE: input$select_all_pop1 is an actionButton (numeric). The equality
      # check against "Select all" is kept to preserve original behavior.
      if (input$select_all_pop1 == "Select all") {
        system(paste0(
          "www/plink19 --bfile ", bfile_path,
          " --make-bed --out ", file.path(pca_dir, "step1")
        ))
      } else {
        system(paste0(
          "www/plink19 --bfile ", bfile_path,
          " --keep ", file.path(pca_dir, "to_keep.txt"),
          " --make-bed --out ", file.path(pca_dir, "step1")
        ))
      }
      
      update_progress(0.20)
      system(paste0(
        "www/plink19 --bfile ", file.path(pca_dir, "step1"),
        " --maf ", input$sliderM1,
        " --make-bed --out ", file.path(pca_dir, "step2")
      ))
      
      update_progress(0.20)
      system(paste0(
        "www/plink19 --bfile ", file.path(pca_dir, "step2"),
        " --indep-pairwise 200 25 0.3 --out ", file.path(pca_dir, "prune1")
      ))
      
      update_progress(0.20)
      system(paste0(
        "www/plink19 --bfile ", file.path(pca_dir, "step2"),
        " --extract ", file.path(pca_dir, "prune1.prune.in"),
        " --make-bed --out ", file.path(pca_dir, "step3")
      ))
      
      update_progress(0.20)
      system(paste0(
        "www/plink19 --bfile ", file.path(pca_dir, "step3"),
        " --pca 10 header tabs var-wts --out ", file.path(pca_dir, "input_PCA")
      ))
      
      update_progress(0.20)
      output$status <- renderText("✅  PCA analysis completed successfully.")
    }, error = function(e) {
      showNotification("PLINK analysis failed. Check the logs.", type = "error", duration = 10)
      output$status <- renderText(paste("Error:", e$message))
    })
    
    # Read PCA eigenvectors and prepare data for plotting
    PCA_results <- fread(file.path(pca_dir, "input_PCA.eigenvec"), header = TRUE)
    
    to_color <- fread(file.path(pca_dir, "to_color.txt"), header = FALSE)
    colnames(to_color) <- c("POP1", "POP2")
    
    PCA_to_plot <- cbind(to_color, PCA_results %>% select(PC1, PC2))
    PCA_data(PCA_to_plot)
    
    print("✅ `PCA_data()` has been updated successfully.")
    shinyjs::enable("run_admix")
    shinyjs::enable("run_fst")
    showNotification("PCA finished. You can now run ADMIXTURE or FST.", type = "message", duration = 10)
  })
  
  # Show a useful line from PLINK PCA log
  output$pca_log <- renderPrint({
    fp <- folder_path()
    if (is.null(fp) || is.na(fp) || fp == "") {
      cat("⚠️ Sample size and variants.")
      return()
    }
    
    log_file <- file.path(fp, "pca", "input_PCA.log")
    if (!file.exists(log_file)) {
      cat("❌ Log file not found at:", log_file)
      return()
    }
    
    log_lines <- readLines(log_file)
    matched   <- grep("people pass filters and QC", log_lines, value = TRUE)
    
    if (length(matched) == 0) {
      cat("⚠️ Line not found in log.")
    } else {
      cat(matched)
    }
  })
  
  # PCA plot
  output$PCA_plot <- renderPlot({
    req(PCA_data())
    req(pca_files())
    
    PCA_to_plot <- PCA_data()
    pca_dir     <- pca_files()
    
    message("✅ `renderPlot` is generating the plot correctly.")
    print(head(PCA_to_plot))
    
    to_color <- fread(file.path(pca_dir, "to_color.txt"), header = FALSE)
    colnames(to_color) <- c("POP1", "POP2")
    
    to_color$POP1 <- factor(to_color$POP1, levels = unique(to_color$POP1))
    to_color$POP2 <- factor(to_color$POP2, levels = unique(to_color$POP2))
    
    PCA_to_plot <- cbind(to_color, PCA_to_plot %>% select(PC1, PC2))
    
    base_palette <- c(
      "red", "blue", "orange", "green", "yellow", "purple", "coral", "cyan",
      "salmon", "skyblue", "tomato", "magenta", "sandybrown", "gray", "chocolate"
    )
    
    # If there are many groups, expand palette via interpolation
    too_many_groups <- length(unique(PCA_to_plot$POP1)) > 15 || length(unique(PCA_to_plot$POP2)) > 15
    selected_palette <- if (too_many_groups) colorRampPalette(base_palette)(74) else base_palette
    
    color_groups <- if (input$popcolor == "POP1") {
      sort(unique(as.character(PCA_to_plot$POP1)))
    } else {
      sort(unique(as.character(PCA_to_plot$POP2)))
    }
    
    palette_named <- setNames(selected_palette[seq_along(color_groups)], color_groups)
    
    if (input$popcolor == "POP1") {
      PCA_to_plot$POP1 <- factor(PCA_to_plot$POP1, levels = color_groups)
      PCA <- ggplot(PCA_to_plot, aes(x = PC1, y = PC2, color = POP1)) +
        geom_point(size = 3, shape = 1, stroke = 1.5) +
        labs(title = "PCA Plot", x = "PC1", y = "PC2") +
        guides(color = guide_legend(title = "Population")) +
        scale_color_manual(values = palette_named) +
        theme_minimal(base_size = 14)
    } else {
      PCA_to_plot$POP2 <- factor(PCA_to_plot$POP2, levels = color_groups)
      PCA_to_plot$POP1 <- factor(PCA_to_plot$POP1, levels = sort(unique(as.character(PCA_to_plot$POP1))))
      PCA <- ggplot(PCA_to_plot, aes(x = PC1, y = PC2, color = POP2, shape = POP1)) +
        geom_point(size = 3, stroke = 1.5) +
        labs(title = "PCA Plot", x = "PC1", y = "PC2") +
        guides(color = guide_legend(title = "Subpopulation"),
               shape  = guide_legend(title = "Population")) +
        scale_color_manual(values = palette_named) +
        theme_minimal(base_size = 14)
    }
    
    message("Colors used:")
    print(palette_named)
    
    # Save PNG
    plot_dir <- file.path(dirname(pca_dir), "plots")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(plot_dir, "PCA.png"), plot = PCA)
    
    return(PCA)
  })
  
  ##########################################################################################
  # ===================================== FST STEP ========================================
  ##########################################################################################
  
  # Global reactives for FST
  folder_path   <- reactiveVal(NULL)
  fst_out_files <- reactiveVal(NULL)
  fst_trigger   <- reactiveVal(0)
  
  # Run FST analysis
  observeEvent(input$run_fst, {
    req(input$run_pca)     # ensure PCA was run (button pressed)
    req(folder_path())     # ensure working folder exists
    
    cat("📁 folder_path():", folder_path(), "\n")
    cat("📄 Expected pop file:", file.path(folder_path(), "pca", "pop_data.txt"), "\n")
    cat("📄 Expected output:",   file.path(folder_path(), "fst", "FST_output.fst"), "\n")
    
    progress <- shiny::Progress$new()
    progress$set(message = "Running FST analysis", value = 0)
    on.exit(progress$close(), add = TRUE)
    
    pca_files   <- file.path(folder_path(), "pca")
    b_files     <- file.path(pca_files, "step3")
    fst_files   <- file.path(folder_path(), "fst")
    fst_in_path <- file.path(fst_files, "to_fst.txt")
    fst_out_path <- file.path(fst_files, "FST_output")
    
    # Build sub-population file: (FID IID POP2?) based on original structure
    pop2   <- fread(file.path(pca_files, "pop_data.txt"), header = TRUE)
    to_fst <- subset(pop2, POP2 %in% input$pop2)[, c(3, 4, 2)]
    
    write.table(to_fst, fst_in_path, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cmd <- paste("www/gcta64 --bfile", b_files,
                 "--fst --sub-popu", fst_in_path,
                 "--out", fst_out_path)
    cat("📣 Running FST:", cmd, "\n")
    system(cmd)
    
    fst_out_files(fst_out_path)
    fst_trigger(fst_trigger() + 1)
  })
  
  # Load FST results
  fst_data <- reactive({
    req(fst_trigger())
    req(fst_out_files())
    
    df <- fread(paste0(fst_out_files(), ".fst"), header = TRUE)
    df <- df[complete.cases(df[, c("Chr", "bp", "Fst")]), ]
    df <- df[is.finite(df$Fst), ]
    df$Chr <- as.numeric(as.character(df$Chr))
    df$bp  <- as.numeric(as.character(df$bp))
    df$Fst <- as.numeric(as.character(df$Fst))
    df <- df[order(df$Chr, df$bp), ]
    
    # Cumulative genomic position for Manhattan-like plotting
    chr_lengths <- df %>%
      dplyr::group_by(Chr) %>%
      dplyr::summarise(chr_len = max(bp)) %>%
      dplyr::arrange(Chr) %>%
      dplyr::mutate(offset = dplyr::lag(cumsum(chr_len), default = 0))
    
    df <- df %>%
      dplyr::left_join(chr_lengths, by = "Chr") %>%
      dplyr::mutate(pos_cum = bp + offset)
    
    # dbSNP links for hover/click
    df$dbSNP_link <- paste0("https://www.ncbi.nlm.nih.gov/snp/", df$SNP)
    
    df
  })
  
  # Interactive Manhattan plot
  output$manhattan_plot <- renderPlotly({
    req(fst_data())
    fst_in_path <- file.path(folder_path(), "fst", "to_fst.txt")
    df <- fst_data()
    
    pop_names <- fread(fst_in_path, header = FALSE)
    names_pop <- unique(pop_names$V3)
    
    df$color_group <- as.factor(df$Chr %% 2)
    
    chr_centers <- df %>%
      dplyr::group_by(Chr) %>%
      dplyr::summarise(center = mean(pos_cum)) %>%
      dplyr::arrange(Chr)
    
    p <- plot_ly(
      data  = df,
      x     = ~pos_cum,
      y     = ~Fst,
      type  = "scatter",
      mode  = "markers",
      color = ~color_group,
      colors = c("black", "red"),
      text  = ~paste0(
        "SNP: ", SNP,
        "<br>Chr: ", Chr,
        "<br>Position: ", bp,
        "<br>Fst: ", round(Fst, 4),
        "<br><a href='", dbSNP_link, "' target='_blank'>View in dbSNP</a>"
      ),
      hoverinfo = "text",
      marker = list(size = 6)
    ) %>%
      plotly::layout(
        title = paste0("Interactive Manhattan plot of FST scores — populations: ", names_pop),
        xaxis = list(
          title   = "Chromosome",
          tickmode = "array",
          tickvals = chr_centers$center,
          ticktext = as.character(chr_centers$Chr)
        ),
        yaxis = list(title = "FST"),
        showlegend = FALSE
      )
    
    # Ensure plots directory exists (optional save disabled as in original)
    temp_dir <- file.path(folder_path(), "plots")
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    # htmlwidgets::saveWidget(p, file.path(temp_dir, "FST_interactive.html"), selfcontained = TRUE)
    
    p
  })
  
  ###########################################################################
  
  # Top 5 FST table
  output$top5_fst <- renderPrint({
    req(fst_data())
    table <- fst_data() %>% dplyr::select(1, 2, 3, 4, 5, 6, 7, 11)
    head(table %>% dplyr::arrange(dplyr::desc(Fst)), 5)
  })
  
  # UCSC and dbSNP links when clicking a SNP in the plot
  observeEvent(event_data("plotly_click"), {
    ed <- event_data("plotly_click")
    df <- fst_data()
    if (!is.null(ed)) {
      # Nearest point by cumulative position
      nearest_row <- df[which.min(abs(df$pos_cum - ed$x)), ]
      chr <- nearest_row$Chr
      pos <- nearest_row$bp
      start_ucsc <- max(0, pos - 10000)
      end_ucsc   <- pos + 10000
      
      ucsc_url <- paste0(
        "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr",
        chr, ":", start_ucsc, "-", end_ucsc
      )
      
      output$dbsnp_link <- renderUI({
        tagList(
          tags$a(href = nearest_row$dbSNP_link, target = "_blank",
                 paste("View SNP", nearest_row$SNP, "in dbSNP")),
          tags$br(),
          tags$a(href = ucsc_url, target = "_blank",
                 paste("View ±10kb region in UCSC (chr", chr, ":", start_ucsc, "-", end_ucsc, ")"))
        )
      })
    }
  })
  
  ##########################################################################################
  # ================================== ADMIXTURE STEP =====================================
  ##########################################################################################
  
  observeEvent(input$run_admix, {
    req(input$run_pca)  # Ensure PCA ran first
    
    # Progress UI
    progress <- shiny::Progress$new()
    progress$set(message = "Running ADMIXTURE analysis. This may take a while...", value = 0)
    on.exit(progress$close(), add = TRUE)
    update_progress <- function(amount) { progress$inc(amount) }
    
    output$status <- renderText("✅  Running ADMIXTURE analysis...")
    
    admx_files    <- file.path(folder_path(), "admx")
    pca_files     <- file.path(folder_path(), "pca")
    # admixture_bin <- file.path(getwd(), "www", "admixture")
    # admixture_bin <- "/home/jfibla/miniconda3/bin/admixture"
    admixture_bin <- normalizePath("www/admixture")
    
    # Adjust MAF and build ADMIXTURE input from PCA step3
    m <- input$sliderM2
    system(paste0(
      "www/plink19 --bfile ", pca_files, "/step3",
      " --maf ", m, " --mind 1.0 --make-bed --out ", admx_files, "/input_admx"
    ))
    # system(paste0("/srv/shiny-server/bpga_VM/www/plink19 --bfile ",pca_files,"/step3 --maf ", m, " --mind 1.0 --make-bed  --out ",admx_files,"/input_admx"))
    
    # Keep original progress behavior
    update_progress(0.05 / (input$sliderM2 - 0.05))
    
    # Show MAF log summary
    log_file <- file.path(admx_files, "input_admx.log")
    if (file.exists(log_file)) {
      keyword      <- "variants and"
      grep_command <- paste("grep", shQuote(keyword), shQuote(log_file))
      log_lines    <- tryCatch(system(grep_command, intern = TRUE), error = function(e) "No matching lines found.")
      output$maf_output <- renderText({
        paste0("✅ Number of SNPs and sample size:\n", paste(log_lines, collapse = "\n"))
      })
    } else {
      output$maf_output <- renderText({ "Log file not found. Ensure PLINK ran successfully." })
    }
    
    # Use all but one CPU core
    Np <- parallel::detectCores() - 1
    
    ##################
    # (Intentionally repeated to preserve original behavior)
    admixture_bin <- normalizePath("www/admixture")
    
    # Run ADMIXTURE for a single K or for K=1..sliderK
    if (input$Kselect == "one") {
      r <- input$sliderK
      system(paste0(
        "cd ", admx_files, " && ", admixture_bin,
        " --cv -j", Np, " input_admx.bed ", r, " | tee log", r, ".out"
      ))
      update_progress(1.0 / (input$sliderK - 1))
    } else {
      for (i in 1:input$sliderK) {
        system(paste0(
          "cd ", admx_files, " && ", admixture_bin,
          " --cv -j", Np, " input_admx.bed ", i, " | tee log", i, ".out"
        ))
        update_progress(1.0 / (input$sliderK - 1))
      }
    }
    #########
    
    # Read all logs from admx_files
    all_logs <- list()
    if (input$Kselect == "one") {
      r <- input$sliderK
      log_file <- file.path(admx_files, paste0("log", r, ".out"))
      if (file.exists(log_file)) {
        all_logs[[paste0("K", r)]] <- readLines(log_file)
      } else {
        warning(paste("Log file not found:", log_file))
      }
    } else {
      for (r in 1:input$sliderK) {
        log_file <- file.path(admx_files, paste0("log", r, ".out"))
        if (file.exists(log_file)) {
          all_logs[[paste0("K", r)]] <- readLines(log_file)
        } else {
          warning(paste("Log file not found:", log_file))
        }
      }
    }
    
    # Extract lines containing CV error
    cv_error_lines <- unlist(lapply(all_logs, function(log) grep("CV error", log, value = TRUE)))
    
    output$cvErrors <- renderPrint({
      if (length(cv_error_lines) > 0) {
        cat(cv_error_lines, sep = "\n")
      } else {
        cat("No CV error lines found.")
      }
    })
    
    # Update the cluster selector (first definition kept)
    output$sorter <- renderUI({
      selectInput("sorter", "Sort admixture boxplot by 'Ancestry component': ",
                  choices = 1:input$sliderK, selected = 1)
    })
    
    # Remove residual files if needed (kept as in original)
    unlink(paste0(
      getwd(), "/", c("1","2","3","4","5","6","7","8","9","10","log","*.P","log*","*.out")
    ), recursive = FALSE)
    
    ########################################
    # Dynamic UI for selecting the ancestry component (second definition kept)
    output$sorter <- renderUI({
      selectInput("sorter", "Sort admixture boxplot by 'Ancestry component': ",
                  choices = 1:input$sliderK, multiple = FALSE, selected = 1)
    })
    ##############################################################################
    
    # Render the stacked ADMIXTURE plot
    output$ADMIXTURE_plot <- renderPlot({
      # Load data
      runs     <- fread(paste0(admx_files, "/input_admx.", input$sliderK, ".Q"))
      to_color <- fread(paste0(pca_files, "/to_color.txt"), header = FALSE)
      colnames(to_color) <- c("POP1", "POP2")
      
      # Combine population and ancestry data
      admixture_data <- cbind(to_color, runs)
      data <- admixture_data
      
      # 10-color palette for K components
      k_palette <- c(
        "#4682B4", "#FFD700", "#00FF00", "#FF69B4", "#556B2F",
        "#7FFFD4", "#9932CC", "#DEB887", "#B22222", "#008080"
      )
      
      # Order by the first ancestry component (V1) within each population group
      if (input$popcolor == "POP1") {
        data_sorted <- data %>%
          dplyr::group_by(POP1) %>%
          dplyr::arrange(V1, .by_group = TRUE) %>%
          dplyr::ungroup()
        
        # Add per-group index
        data_sorted <- data_sorted %>%
          dplyr::group_by(POP1) %>%
          dplyr::mutate(Index = dplyr::row_number()) %>%
          dplyr::ungroup()
        
        # Long format
        data_long <- data_sorted %>%
          tidyr::pivot_longer(cols = starts_with("V"), names_to = "Variable", values_to = "Value")
        
        # Sort POP1 alphabetically
        data_long$POP1 <- factor(data_long$POP1, levels = sort(unique(as.character(data_long$POP1))))
        
        # Plot
        ADM <- ggplot(data_long, aes(x = Index, y = Value, fill = Variable)) +
          geom_bar(stat = "identity", position = "stack") +
          ggh4x::facet_nested(~ POP1, scales = "free_x") +
          labs(x = "Individuals", y = "Ancestry component", title = "Admixture Plot", fill = "Components") +
          scale_fill_manual(values = k_palette) +
          theme(
            plot.title   = element_text(size = 18, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text    = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text  = element_text(size = 12)
          )
        
      } else {
        data_sorted <- data %>%
          dplyr::group_by(POP2) %>%
          dplyr::arrange(V1, .by_group = TRUE) %>%
          dplyr::ungroup()
        
        # Add per-group index
        data_sorted <- data_sorted %>%
          dplyr::group_by(POP2) %>%
          dplyr::mutate(Index = dplyr::row_number()) %>%
          dplyr::ungroup()
        
        # Long format
        data_long <- data_sorted %>%
          tidyr::pivot_longer(cols = starts_with("V"), names_to = "Variable", values_to = "Value")
        
        # Sort POP2 alphabetically
        data_long$POP2 <- factor(data_long$POP2, levels = sort(unique(as.character(data_long$POP2))))
        
        # Plot
        ADM <- ggplot(data_long, aes(x = Index, y = Value, fill = Variable)) +
          geom_bar(stat = "identity", position = "stack") +
          ggh4x::facet_nested(~ POP1 + POP2, scales = "free_x", space = "free") +
          labs(x = "Individuals", y = "Ancestry component", title = "Admixture Plot", fill = "Components") +
          scale_fill_manual(values = k_palette) +
          theme(
            plot.margin  = margin(10, 10, 10, 10),
            panel.spacing = unit(0, "lines"),
            strip.text.x = element_text(size = 14, face = "bold", angle = 0),
            strip.text.y = element_text(size = 14, face = "bold", angle = 45),
            plot.title   = element_text(size = 18, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text    = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text  = element_text(size = 12)
          )
      }
      
      print(ADM)
      
      # Save PNG
      temp_dir <- file.path(folder_path(), "plots")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
      ggsave(paste0(temp_dir, "/ADM.png"), plot = ADM)
    })
    
    ##############################################################################
    # Render the ADMIXTURE boxplot
    output$box_plot <- renderPlot({
      # Load data
      runs     <- fread(paste0(admx_files, "/input_admx.", input$sliderK, ".Q"))
      to_color <- fread(paste0(pca_files, "/to_color.txt"), header = FALSE)
      colnames(to_color) <- c("POP1", "POP2")
      
      # Combine population and ancestry data
      admixture_data <- cbind(to_color, runs)
      
      # Compute mean values of the selected component across POP2
      for (i in 1:input$sorter) {
        ii <- paste0("V", i)
        MEAN_POP2 <- aggregate(admixture_data[[ii]], by = list(admixture_data$POP2), FUN = mean)
      }
      # Order by mean and merge back
      MEAN_POP2 <- MEAN_POP2[order(MEAN_POP2$x), ]
      MEAN_POP2$Group.1
      level_orderX  <- c(MEAN_POP2$Group.1)
      admixture_data <- merge(admixture_data, MEAN_POP2, by.x = "POP2", by.y = "Group.1", all.x = TRUE)
      
      # Sort by the auxiliary mean column
      admixture_data <- admixture_data[order(admixture_data[["x"]]), ]
      admixture_data <- admixture_data %>% dplyr::relocate(x, .before = V1)
      
      # Long format for ggplot
      admixture_long <- tidyr::pivot_longer(admixture_data, cols = 4:ncol(admixture_data))
      colnames(admixture_long) <- c("POP2", "POP1", "x", "K", "value")
      
      # Color palette
      if (length(unique(admixture_long$POP2)) <= 15) {
        selected_palete <- c(
          "red", "blue", "orange", "green",
          "yellow", "purple", "coral", "cyan",
          "salmon", "skyblue", "tomato", "magenta",
          "sandybrown", "gray", "chocolate"
        )
      } else {
        selected_palete <- colorRampPalette(c("red", "blue", "orange", "green", "yellow", "purple"))(74)
      }
      
      admixture_long$POP1 <- factor(admixture_long$POP1, levels = sort(unique(as.character(admixture_long$POP1))))
      admixture_long$POP2 <- factor(admixture_long$POP2, levels = sort(unique(as.character(admixture_long$POP2))))
      print(admixture_long)
      
      if (input$popcolor == "POP1") {
        p <- ggplot(admixture_long, aes(x = K, y = value, fill = POP1)) +
          geom_boxplot() +
          labs(x = "Ancestry component", y = "Ancestry Proportion", fill = "Population") +
          scale_fill_manual(values = selected_palete) +
          theme(
            plot.title   = element_text(size = 18, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text    = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text  = element_text(size = 12)
          )
      } else {
        p <- ggplot(admixture_long, aes(x = K, y = value, fill = POP2)) +
          geom_boxplot() +
          ggh4x::facet_nested(~ POP1, scales = "free", space = "free") +
          # facet_grid(cols = vars(POP1), scales = "free", space = "free")
          scale_fill_manual(values = selected_palete) +
          labs(x = "Ancestry component", y = "Ancestry Proportion", fill = "Population") +
          theme(
            plot.margin  = margin(10, 10, 10, 10),
            panel.spacing = unit(0, "lines"),
            strip.text.x = element_text(size = 14, face = "bold", angle = 0),
            strip.text.y = element_text(size = 14, face = "bold", angle = 45),
            plot.title   = element_text(size = 18, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text    = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text  = element_text(size = 12)
          )
      }
      
      print(p)
      
      # Save PNG
      temp_dir <- file.path(folder_path(), "plots")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
      ggsave(paste0(temp_dir, "/box_ADM.png"), plot = p)
    })
    
    ################################################################################
    # ================================== GEO MAP STEP ==============================
    ################################################################################
    map_with_save <- function() {
      req(input$sliderK)
      print("map_with_save defined correctly")
      
      # Paths
      admx_files <- file.path(folder_path(), "admx")
      pca_files  <- file.path(folder_path(), "pca")
      
      # Read and prepare coordinate file
      coord <- fread("www/coord.txt", header = TRUE)
      
      if (!is.null(input$selector)) {
        if (input$files %in% c("exp1", "exp2")) {
          coord
        } else if (input$files == "usr") {
          if (!is.null(input$fam_formated) && input$fam_formated == "unic") {
            user_data <- data.frame(
              POP1     = paste0(input$popID, "*"),
              POP2     = paste0(input$popID, "*"),
              Pop_name = "user population",
              LAT      = as.numeric(input$lat),
              LON      = as.numeric(input$lon)
            )
            coord <- rbind(coord, user_data)
            
          } else if (!is.null(input$fam_formated) && input$fam_formated %in% c("pop1pop2", "subpop")) {
            user_file <- file.path(folder_path(), "user", "user_coord.txt")
            user_data <- fread(user_file, sep = "\t", header = TRUE)
            colnames(user_data)[1:2] <- c("POP1", "POP2")
            coord <- rbind(coord, user_data)
          }
        }
      }
      
      print("coord preview:")
      print(coord)
      
      # Read Q-matrix
      runs <- fread(paste0(admx_files, "/input_admx.", input$sliderK, ".Q"))
      value_columns <- paste0("V", 1:input$sliderK)
      colnames(runs) <- value_columns
      
      # Population labels
      to_color <- fread(paste0(pca_files, "/to_color.txt"), header = FALSE)
      colnames(to_color) <- c("POP1", "POP2")
      
      # Merge Q and population data
      admix_data <- cbind(to_color, runs)
      print("Admixture input file:")
      print(admix_data)
      
      # Mean ancestry proportions by subpopulation (POP2) and superpopulation (POP1)
      qpop <- aggregate(. ~ POP2, data = as.data.frame(admix_data)[, c("POP2", value_columns)], mean)
      QPOP <- aggregate(. ~ POP1, data = as.data.frame(admix_data)[, c("POP1", value_columns)], mean)
      print("Admixture mean qpop (by POP2):")
      print(qpop)
      print("Admixture mean QPOP (by POP1):")
      print(QPOP)
      
      # Count individuals by POP2 to control pie size
      pop_counts <- to_color %>% dplyr::count(POP2, name = "pop_size")
      pop_counts$pop_size <- as.numeric(pop_counts$pop_size)
      
      # Merge means with counts
      admixture_data <- merge(qpop, pop_counts, by = "POP2")
      print("admixture_data -> qpop + pop_counts:")
      print(admixture_data)
      
      # Merge with coordinates (POP2-level)
      coord_pop2_df <- merge(admixture_data, coord, by = "POP2")
      print("POP2 means + coordinates:")
      print(coord_pop2_df)
      
      # Aggregate to POP1-level for superpop map
      coord_pop1_df <- coord_pop2_df %>%
        dplyr::group_by(POP1) %>%
        dplyr::summarise(
          Pop_name = POP1,
          dplyr::across(starts_with("V"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
          LAT      = mean(LAT, na.rm = TRUE),
          LON      = mean(LON, na.rm = TRUE),
          pop_size = sum(pop_size, na.rm = TRUE),
          .groups  = "drop"
        )
      print("QPOP means + coordinates:")
      print(coord_pop1_df)
      
      # Choose which table to map
      if (input$popcolor == "POP1") {
        coord_pop_df <- coord_pop1_df
      } else {
        coord_pop_df <- coord_pop2_df
      }
      
      # Convert to sf and compute scaled marker sizes
      df_sf <- coord_pop_df %>%
        sf::st_as_sf(coords = c("LON", "LAT"), crs = 4326, remove = FALSE) %>%
        dplyr::mutate(
          LON = sf::st_coordinates(.)[, 1],
          LAT = sf::st_coordinates(.)[, 2]
        ) %>%
        dplyr::mutate(scaled_size = scales::rescale(pop_size, to = c(40, input$sliderR)))
      
      print("Map input (sf):")
      print(df_sf)
      
      # Color palette for pie slices
      k_palette <- c(
        "#4682B4", "#FFD700", "#00FF00", "#FF69B4", "#556B2F",
        "#7FFFD4", "#9932CC", "#DEB887", "#B22222", "#008080"
      )
      
      chartdata    <- coord_pop_df %>% dplyr::select(dplyr::starts_with("V")) %>% as.matrix()
      num_pieces   <- ncol(chartdata)
      dynamic_cols <- k_palette[1:num_pieces]
      
      leaflet_map <- leaflet::leaflet(data = df_sf) %>%
        leaflet::addTiles() %>%
        leaflet.minicharts::addMinicharts(
          lng        = df_sf$LON,
          lat        = df_sf$LAT,
          chartdata  = chartdata,
          type       = "pie",
          width      = df_sf$scaled_size,
          colorPalette = dynamic_cols
        ) %>%
        leaflet::addLabelOnlyMarkers(
          lng = df_sf$LON,
          lat = df_sf$LAT,
          label = lapply(df_sf$Pop_name, htmltools::HTML),
          labelOptions = leaflet::labelOptions(noHide = TRUE, direction = "top", textOnly = TRUE, textsize = "12px")
        )
      
      temp_dir <- file.path(folder_path(), "plots")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
      
      leaflet_map
    }
    
    # Render the leaflet map
    output$map_plot <- renderLeaflet({
      map_with_save()
    })
  })
  # ============================== Launch ADMIXTURE (manual) ===============================
  observeEvent(input$lancaAdmixture, {
    req(input$fitxer, input$k)
    
    admixture_bin <- "www/admixture"
    
    # Compose input paths and log path
    fitxer   <- paste0("data/", input$fitxer)
    k        <- input$k
    log_path <- paste0("logs/", tools::file_path_sans_ext(input$fitxer), ".", k, ".log")
    
    # Run command and capture stdout (errors may not be captured with system(..., intern=TRUE))
    cmd <- paste(admixture_bin, fitxer, k)
    log_output <- tryCatch(
      {
        system(cmd, intern = TRUE)
      },
      error = function(e) {
        paste("Execution error:", e$message)
      }
    )
    
    # Write output to log file
    writeLines(log_output, log_path)
  })
  
  # ============================== Show ADMIXTURE log ======================================
  output$logAdmixture <- renderText({
    req(input$fitxer, input$k)
    log_path <- paste0("logs/", tools::file_path_sans_ext(input$fitxer), ".", input$k, ".log")
    if (file.exists(log_path)) {
      paste(readLines(log_path), collapse = "\n")
    } else {
      "No log has been generated yet for this file/K."
    }
  })
  
  ##############################################################################
  ############################ DOWNLOAD RESULTS STEP ###########################
  ##############################################################################
  
  # ---- Helpers to generate artifacts into plots/ before zipping ----
  
  generate_admix_map <- function() {
    fp <- folder_path()
    validate(need(!is.null(fp) && nzchar(fp), "Base folder not set."))
    
    admx_files <- file.path(fp, "admx")
    pca_files  <- file.path(fp, "pca")
    out_dir    <- file.path(fp, "plots")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Inputs needed for map
    req(input$sliderK)
    value_columns <- paste0("V", 1:input$sliderK)
    
    # Coordinates
    coord <- data.table::fread("www/coord.txt", header = TRUE)
    
    # Add user coordinates when applicable
    if (!is.null(input$selector) && input$files == "usr") {
      if (!is.null(input$fam_formated) && input$fam_formated == "unic") {
        user_data <- data.frame(
          POP1     = paste0(input$popID, "*"),
          POP2     = paste0(input$popID, "*"),
          Pop_name = "user population",
          LAT      = as.numeric(input$lat),
          LON      = as.numeric(input$lon)
        )
        coord <- rbind(coord, user_data)
      } else if (!is.null(input$fam_formated) && input$fam_formated %in% c("pop1pop2", "subpop")) {
        user_file <- file.path(fp, "user", "user_coord.txt")
        if (file.exists(user_file)) {
          user_data <- data.table::fread(user_file, sep = "\t", header = TRUE)
          colnames(user_data)[1:2] <- c("POP1", "POP2")
          coord <- rbind(coord, user_data, fill = TRUE)
        }
      }
    }
    
    # Q-matrix (K components)
    runs <- data.table::fread(file.path(admx_files, paste0("input_admx.", input$sliderK, ".Q")))
    colnames(runs) <- value_columns
    
    # Population labels used during PCA
    to_color <- data.table::fread(file.path(pca_files, "to_color.txt"), header = FALSE)
    colnames(to_color) <- c("POP1", "POP2")
    
    # Merge and compute means
    admix_data <- cbind(to_color, runs)
    
    qpop <- aggregate(. ~ POP2, data = as.data.frame(admix_data)[, c("POP2", value_columns)], mean)
    QPOP <- aggregate(. ~ POP1, data = as.data.frame(admix_data)[, c("POP1", value_columns)], mean)
    
    # Counts by POP2 to scale pies
    pop_counts <- to_color %>% dplyr::count(POP2, name = "pop_size")
    pop_counts$pop_size <- as.numeric(pop_counts$pop_size)
    
    admixture_data <- merge(qpop, pop_counts, by = "POP2")
    
    # Merge with coordinates
    coord_pop2_df <- merge(admixture_data, coord, by = "POP2")
    
    # Aggregate to POP1 for superpopulation map
    coord_pop1_df <- coord_pop2_df %>%
      dplyr::group_by(POP1) %>%
      dplyr::summarise(
        Pop_name = POP1,
        dplyr::across(dplyr::starts_with("V"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
        LAT      = mean(LAT, na.rm = TRUE),
        LON      = mean(LON, na.rm = TRUE),
        pop_size = sum(pop_size, na.rm = TRUE),
        .groups  = "drop"
      )
    
    # Choose which table to map
    coord_pop_df <- if (input$popcolor == "POP1") coord_pop1_df else coord_pop2_df
    
    df_sf <- coord_pop_df %>%
      sf::st_as_sf(coords = c("LON", "LAT"), crs = 4326, remove = FALSE) %>%
      dplyr::mutate(
        LON = sf::st_coordinates(.)[, 1],
        LAT = sf::st_coordinates(.)[, 2],
        scaled_size = scales::rescale(pop_size, to = c(40, input$sliderR))
      )
    
    # Color palette for K slices
    k_palette <- c(
      "#4682B4", "#FFD700", "#00FF00", "#FF69B4", "#556B2F",
      "#7FFFD4", "#9932CC", "#DEB887", "#B22222", "#008080"
    )
    chartdata    <- coord_pop_df %>% dplyr::select(dplyr::starts_with("V")) %>% as.matrix()
    num_pieces   <- ncol(chartdata)
    dynamic_cols <- k_palette[1:num_pieces]
    
    m <- leaflet::leaflet(data = df_sf) %>%
      leaflet::addTiles() %>%
      leaflet.minicharts::addMinicharts(
        lng = df_sf$LON, lat = df_sf$LAT,
        chartdata = chartdata,
        type = "pie", width = df_sf$scaled_size,
        colorPalette = dynamic_cols
      ) %>%
      leaflet::addLabelOnlyMarkers(
        lng = df_sf$LON, lat = df_sf$LAT,
        label = lapply(df_sf$Pop_name, htmltools::HTML),
        labelOptions = leaflet::labelOptions(noHide = TRUE, direction = "top", textOnly = TRUE, textsize = "12px")
      )
    
    # Save as HTML (portable in ZIP)
    htmlwidgets::saveWidget(m, file.path(out_dir, "ADMIXTURE_map.html"), selfcontained = TRUE)
  }
  
  generate_fst_plot <- function() {
    fp <- folder_path()
    validate(need(!is.null(fp) && nzchar(fp), "Base folder not set."))
    
    out_dir <- file.path(fp, "plots")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Use the existing reactive FST data
    df <- tryCatch(isolate(fst_data()), error = function(e) NULL)
    validate(need(!is.null(df) && nrow(df) > 0, "No FST data available."))
    
    # Names for title (if available)
    fst_in_path <- file.path(fp, "fst", "to_fst.txt")
    names_pop <- ""
    if (file.exists(fst_in_path)) {
      pop_names <- data.table::fread(fst_in_path, header = FALSE)
      names_pop <- paste(unique(pop_names$V3), collapse = ", ")
    }
    
    df$color_group <- as.factor(df$Chr %% 2)
    
    chr_centers <- df %>%
      dplyr::group_by(Chr) %>%
      dplyr::summarise(center = mean(pos_cum)) %>%
      dplyr::arrange(Chr)
    
    p <- plotly::plot_ly(
      data = df,
      x = ~pos_cum, y = ~Fst,
      type = "scatter", mode = "markers",
      color = ~color_group, colors = c("black", "red"),
      text = ~paste0(
        "SNP: ", SNP,
        "<br>Chr: ", Chr,
        "<br>Position: ", bp,
        "<br>FST: ", round(Fst, 4),
        "<br><a href='", dbSNP_link, "' target='_blank'>View in dbSNP</a>"
      ),
      hoverinfo = "text",
      marker = list(size = 6)
    ) %>%
      plotly::layout(
        title = paste0("Interactive Manhattan plot of FST — populations: ", names_pop),
        xaxis = list(
          title = "Chromosome",
          tickmode = "array",
          tickvals = chr_centers$center,
          ticktext = as.character(chr_centers$Chr)
        ),
        yaxis = list(title = "FST"),
        showlegend = FALSE
      )
    
    # Save as HTML (portable in ZIP)
    htmlwidgets::saveWidget(p, file.path(out_dir, "FST_interactive.html"), selfcontained = TRUE)
  }
  
  output$downloadPlotsZip <- downloadHandler(
    filename = function() {
      paste0("All_Plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      temp_dir <- file.path(folder_path(), "plots")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)  # ensure exists
      
      # Try to generate ADMIXTURE map
      tryCatch({
        generate_admix_map()
      }, error = function(e) {
        showNotification(paste("⚠️ Could not generate ADMIXTURE map:", e$message), type = "error")
      })
      
      # Try to generate FST plot
      tryCatch({
        generate_fst_plot()
      }, error = function(e) {
        showNotification(paste("⚠️ Could not generate FST plot:", e$message), type = "error")
      })
      
      # Create ZIP from plots directory
      zipfile <- tempfile(fileext = ".zip")
      old_wd  <- setwd(temp_dir)
      on.exit(setwd(old_wd), add = TRUE)
      zip(zipfile, files = list.files(".", full.names = FALSE), flags = "-r9X")
      file.copy(zipfile, file, overwrite = TRUE)
    }
  )
 
  ################## Clear temp subfolders & show status (for logs/UI) ######################
  
  observeEvent(input$clear_folders, {
    base_folder <- file.path(folder_path) # Ensure folder_path is valid
    
    # Check if base_folder is valid
    if (!dir.exists(base_folder)) {
      output$clear_status <- renderText({
        paste("Base folder does not exist:", base_folder)
      })
      return() # Stop execution if base folder does not exist
    }
    
    subfolders <- c("admx","temporal", "merged", "pca", "plots", "fst", "example_files","temporal")
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
}

# Reset All button functionality
print(sapply(sort(sapply(ls(), function(x) object.size(get(x))), decreasing = TRUE), function(size) round(size / 1024^2, 2)))

##################################################################################
shinyApp(ui = ui, server = server)