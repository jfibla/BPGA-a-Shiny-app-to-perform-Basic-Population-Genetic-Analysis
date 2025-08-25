# BPGA-a-Shiny-app-to-perform-Basic-Population-Genetic-Analysis
# setwd("~/Library/CloudStorage/Dropbox/BPGA_app_folders/bpga_VM/")
# acces to all folders
.libPaths(c("/home/joan_fibla/R/x86_64-pc-linux-gnu-library/4.3", .libPaths()))
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

conflicts_prefer(
  maps::map(),
  dplyr::count())

# Increase max upload size to 600MB
options(shiny.maxRequestSize = 600 * 1024^2)  # 600MB limit


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
      conditionalPanel("output.panel_visible === false",                       #>>>>> cond1 Panel info app
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
      ),                                                                        #>>>>> /cond1 Panel info app
      
      #tags$a(href="https://github.com/jfibla/BPGA-a-Shiny-app-to-perform-Basic-Population-Genetic-Analysis.git", 
      #       h5( "Code available")),
      
      # to generate user folder
      #   actionButton("start_analysis", "Start analysis"),
      actionButton("reload_button", "Reload session"),
      actionButton("info_00", "ℹ"),                       # info about reload session each time to start
      #  h4("Main folders:"),
      verbatimTextOutput("folder_path"),                  # temporal folder activated
   #  verbatimTextOutput("subfolder_paths"),                  # temporal subfolder activated
      conditionalPanel(condition = "output.panel_visible === true",     #>>>>> cond2
                       actionButton("info_0", "ℹ"),             # info about file to be selected
                       selectInput("files", "Select input file:",
                                   c("Example single pop" = "exp1","Example multi pop" = "exp2",
                                     "User file" = "usr")),
                #>>>>> cond3 if user file is used
                conditionalPanel(condition = "input.files == 'usr'",
                                 
                                 fileInput("upload_file", "Load user population data (compressed .zip or .gz)", 
                                           accept = c(".zip", ".gz")),
                                 actionButton("info_23", "ℹ"),
                                 radioButtons(
                                   inputId = "fam_formated",
                                   label = "Format fam file:",
                                   choices = list(
                                     "First column preformated as POP1_POP2" = "pop1pop2",
                                     "First column describe one single population" = "unic",
                                     "First column describe several populations" = "subpop"
                                   ),
                                   selected = character(0)
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.fam_formated == 'unic' || input.fam_formated == 'subpop'",
                                   textInput("popID", "Assign population identifier (use a three-letter code like 'USR')", "USR")
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.fam_formated == 'unic'",
                                   textInput("lat", "Assign latitude"),
                                   textInput("lon", "Assign longitude")
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.fam_formated == 'pop1pop2' || input.fam_formated == 'subpop'",
                                   actionButton("info_22", "ℹ"),
                                   fileInput("upload_coord", "Load population coordinates (.txt)", accept = ".txt")
                                 )
                                 
                ),
                
                #>>>>> /cond3 
                       h3("Process input File"),
                       tags$hr(style="border-color: grey;"),
                       actionButton("info_01", "ℹ"),                                     # info about to procese file alone or merge
                       radioButtons("selector", "Select option:", 
                                    choices = list("Only loaded population" = "usr","Merge with Woldwide populations" = "wwp"),
                                    selected = character(0)
                       ),
                       conditionalPanel(condition = "input.selector == 'usr'",    #>>>>> cond3.1  if user file is selected then
                                        actionButton("decompress_file", "Decompress file"),
                       ),                                                         #>>>>> /cond3.1
                       verbatimTextOutput("wait_ref"),
                       
                       conditionalPanel(condition = "input.selector == 'wwp'",           #>>>>> cond4
                                        verbatimTextOutput("reference_path"),
                                        verbatimTextOutput("info_ref"),
                                        actionButton("run_merge", "Merge files")
                       ),                                                          #>>>>> /cond4
                       verbatimTextOutput("wait_merge"),
                       verbatimTextOutput("merged_output"),
                
    conditionalPanel(condition = "input.decompress_file > 0 || input.run_merge  > 0 " ,        
                       h3("Perform population analysis"),
                       conditionalPanel(
                         condition = "input.decompress_file == true && !(input.run_merge == true)",
                         actionButton("info_2", "ℹ")
                       ),
                       conditionalPanel(
                         condition = "input.run_merge == true && !(input.decompress_file == true)",
                         actionButton("info_2", "ℹ")
                       ),
                       uiOutput("pop1_select"),  # Dynamic SelectInput for POP1
                       shinyjs::hidden(actionButton("select_all_pop1", "Select all")),
                       tags$hr(style="border-color: grey;"),
                       
                       h4("PCA analysis"),
                       actionButton("info_3", "ℹ"),
                       sliderInput("sliderM1", "Set MAF for prunning:", min = 0.05, max = 0.450, value = 0.05, step = 0.05),
                       actionButton("run_pca", "Run PCA analysis"),
                       conditionalPanel(
                         condition = "input.run_pca > 0",
                         verbatimTextOutput("pca_log")
                       ),
                       tags$hr(style="border-color: grey;"),
                       
                       h4("FST analysis: compare two populations"),
                       actionButton("info_5", "ℹ"),
                       uiOutput("pop2_checkbox"),  # Dynamic Checkbox for POP2
                       shinyjs::hidden(actionButton("select_all_pop2", "Remove selected")),
                       actionButton("run_fst", "Run FST analysis"),
                       #verbatimTextOutput("max_fst_output"),  # Output for displaying the max FST rowtags$hr(style="border-color: grey;"),
                       conditionalPanel(condition = "input.run_fst > 0",
                                        tags$h5("Top 5 Fst values"),
                                        verbatimTextOutput("top5_fst"),
                       ),
                       tags$hr(style="border-color: grey;"),
                       h4("Admixture analysis"),
                       actionButton("info_4", "ℹ"),
                       sliderInput("sliderM2", "Set MAF:", min = 0.05, max = 0.450, value = 0.05, step = 0.05),
                       verbatimTextOutput("maf_output"),  
                       actionButton("info_24", "ℹ"),
                       radioButtons("Kselect", "Select K option:", 
                                  choices = list("Single K value" = "one", "Range of K values" = "range")),
                       sliderInput("sliderK", "Set K values:", min = 1, max = 10, value = 1),
                       
                       actionButton("run_admix", "Run ADMIXTURE"),
                       verbatimTextOutput("cvErrors"),
                       verbatimTextOutput("logAdmixture"),
                       tags$hr(style="border-color: grey;"),
                       
                       h4("Decorate plots"),
                       radioButtons("popcolor", "Colored:", 
                                    choices = list("Superpopulation" = "POP1", "Subpopulation" = "POP2")),
                       sliderInput("sliderR", "Set Pie radius:", min = 20, max = 200, value = 100, step = 10),
                       # sliderInput("sliderE", "Expand map:", min = 0, max = 2, value = 0.2, step = 0.05),
                       
                       uiOutput("sorter"),
                       #selectInput("sorter", "Sort admixture boxplot:",
                       #                 c("By first component" = "V1",
                       #                   "By second component" = "V2")),
                       tags$hr(style="border-color: grey;"),
                       h4("Report"),
                     #  downloadButton("downloadReport", "Download Report"),
                       downloadButton("downloadPlotsZip", "Download plots (.zip)"),
                       
      )                                                             #>>>>> /cond2
    )),
    
    ###################################################################### 
    mainPanel(
      
      uiOutput("image_ui"),  # The image UI will be dynamically controlled
      uiOutput("image2_ui"),  # The image UI will be dynamically controlled
      
      textOutput("status"),
      uiOutput("progress_bar_ui"),
      
      # Titles and plots
      
      conditionalPanel(
        h4("PCA plot"),
        condition = "input.run_pca > 0",
        withSpinner(plotOutput("PCA_plot"), type = 1, color = "#007BFF")),
      
      conditionalPanel(
        h4("FST values by chromosome"),
        condition = "input.run_fst > 0",
        withSpinner(plotlyOutput("manhattan_plot"), type = 1, color = "#007BFF"),
        uiOutput("dbsnp_link"),
      ),
      
      conditionalPanel(
        h4("ADMIXTURE plot"),
        condition = "input.run_admix > 0",
        div(
          style = "text-align: center;",  # Centrar contenido
          withSpinner(
            plotOutput("ADMIXTURE_plot"),
            type = 1,
            color = "#007BFF"
          ),
          div(
            id = "loading_message1",
            HTML("Loading data, running Admixture and ploting barplot.<br>It can take a long time, please be patient..."),
            style = "color: #007BFF; font-size: 16px; display: block; margin-top: 10px;"
          ))),
      
      conditionalPanel(
        h4("Boxplot of ADMIXTURE data"),
        condition = "input.run_admix  > 0",
        div(
          style = "text-align: center;",  # Centrar contenido
          withSpinner(
            plotOutput("box_plot"),
            type = 1,
            color = "#007BFF"
          ),
          div(
            id = "loading_message2",
            HTML("Loading data, running Admixture and ploting boxplot.<br>Is time to take a coffe..."),
            style = "color: #007BFF; font-size: 16px; display: block; margin-top: 10px;"
          ))),
      
      conditionalPanel(
        h4("Geographic distribution of ancestry components"),
        condition = "input.run_admix  > 0",
        div(
          style = "text-align: center;",  # Centrar contenido
          withSpinner(
            leafletOutput("map_plot", height = "600px"),  # Usa leafletOutput si usas leaflet
            type = 1,
            color = "#007BFF"
          ),
          div(
            id = "loading_message3",
            HTML("Loading data, running Admixture and ploting data on map.<br>As soon as possible we wil be back..."),
            style = "color: #007BFF; font-size: 16px; display: block; margin-top: 10px;"
          ))),
      
    )
  )
)

################################################################################  
########################### SERVER #############################################  
################################################################################  

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 600 * 1024^2)  # 600MB limit
  # reactivacio de processos
  # Al inicio del server
  fst_trigger <- reactiveVal(0)
  fst_out_files <- reactiveVal()
  
  fst_html_path <- reactiveVal(NULL)
  
  
 # # Crear carpeta temporal para almacenar los plots
 # dir.create("/tmp/plots", showWarnings = FALSE, recursive = TRUE)
 # addResourcePath("plots", "/tmp/plots")
 # ploted_files <- "/tmp/plots/"
  
  temp_plots <- file.path(tempdir(), "plots")
  dir.create(temp_plots, showWarnings = FALSE, recursive = TRUE)
  addResourcePath("plots", temp_plots)
  ploted_files <- temp_plots
  
  # message during admixture process
  observeEvent(input$run_admix, {
    shinyjs::show("loading_message1")  # Muestra el mensaje al iniciar el proceso
    Sys.sleep(10)  # Simula un proceso largo
    shinyjs::hide("loading_message1")  # Oculta el mensaje cuando el proceso termina
  })
  
  observeEvent(input$run_admix, {
    shinyjs::show("loading_message2")  # Muestra el mensaje al iniciar el proceso
    Sys.sleep(10)  # Simula un proceso largo
    shinyjs::hide("loading_message2")  # Oculta el mensaje cuando el proceso termina
  })
  
  observeEvent(input$run_admix, {
    shinyjs::show("loading_message3")  # Muestra el mensaje al iniciar el proceso
    Sys.sleep(10)  # Simula un proceso largo
    shinyjs::hide("loading_message3")  # Oculta el mensaje cuando el proceso termina
  })
  # link to image
  output$image_ui <- renderUI({
    if (input$start_analysis == 0) {
      img(src = "infograph.png", height = "700px", width = "800px")
    } else {
      NULL  # Remove the image after the button is clicked
    }
  })
  output$image2_ui <- renderUI({
    if (input$start_analysis != 0 && input$run_pca == 0) {
      tags$img(src = "workflow.jpg", height = "400px", width = "900px")
    } else {
      NULL
    }
  })
  #################### info mesages ############################################## 
  observeEvent(input$info_00, {
    showModal(modalDialog(
      title = "Info: Select input files",
      p(style="text-align: justify;","Reload session in order to perform a new analysis"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # Info claims 
  observeEvent(input$info_0, {
    showModal(modalDialog(
      title = "Info: Select input files",
      p(style="text-align: justify;","Select example files (single population or multipopulation) or load your own input file (zip compressed file containing",strong(".bed, .bim and .fam"), " files of your population. Upload size is limited to 1GB)"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_01, {
    showModal(modalDialog(
      title = "Info: Process files",
      "Process input file alone or merge it by woldwide reference population. Then click on 'Decompress file' button to continue",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_1, {
    showModal(modalDialog(
      title = "Info: Select reference files",
      p(style = "text-align: justify;", 
        "Select ", strong("Small"), " or ", strong("Large"), " reference files. Selecting ", 
        strong("Large"), " may significantly increase processing time (time to decompress a larger (27G) file). ", strong("Please be patient!")),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_2, {
    showModal(modalDialog(
      title = "Info: Select populations to analyze",
      p(style = "text-align: justify;", 
        "From the expanded list (above) select ",strong("one")," or ",strong("multiple")," populations to be included in the PCA analysis"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_3, {
    showModal(modalDialog(
      title = "Info: Select MAF threshold",
      p(style = "text-align: justify;", 
        "You can select a ", strong("MAF (Minor Allele Frequency)"), " threshold to filter SNPs. Larger MAF values reduce processing time but may result in a loss of accuracy. ",
        br(), em("(Note: Be cautious when selecting MAF = 0.5, as it could result in the total loss of SNPs for analysis.)")
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_4, {
    showModal(modalDialog(
      title = "Info: Admixture analysis",
      p(style = "text-align: justify;", 
        "You can select a ", strong("Minor Allele Frequency (MAF)"), 
        " threshold to filter SNPs and specify ", strong("K values"), 
        " (the number of ancestral populations or clusters hypothesized in the model) ranging from 1 to 10. Higher MAF values reduce processing time but may compromise accuracy.", 
        br(),
        em("Note: Selecting a MAF of 0.5 could result in the complete exclusion of SNPs from the analysis, so proceed with caution. Additionally, processing time increases significantly as K values rise."),
        br(),
        em(span(style = "color: red;","ADMIXTURE analysis may take several minutes (hours) to complete. Take the time to enjoy a coffee while you wait!"))
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_24, {
    showModal(modalDialog(
      title = "Info: Admixture analysis",
      p(style = "text-align: justify;", 
        "You can select a ", strong("single K value"), 
        " or a ", strong("range of K values (from 1 to K selected)"), 
        " (K corresponds to the number of ancestral populations or clusters hypothesized in the model) ranging from 1 to 10.", 
        br(),
        em("Note: Processing time increases significantly from single to range of K values."),
        br(),
        em(span(style = "color: red;","ADMIXTURE analysis may take several minutes (hours) to complete. Take the time to enjoy a coffee while you wait!"))
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  observeEvent(input$info_5, {
    showModal(modalDialog(
      title = "Info: Select populations to analyze",
      p(style = "text-align: justify;", 
        "From the expanded list (above) select a pair of populations to be included in the FST analysis"),
      # br(),
      em("Note: You can select more than two populations, then FST value reflects the degree of genetic differentiation across all populations considered."),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_21, {
    showModal(modalDialog(
      title = "Info: Fam file formats",
      HTML(
        '<p style="text-align: justify;">
    Fam file must be formated as: 
    <strong>Preformated POP1_POP2:</strong>
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
      title = "Info: Load coordinate data of user populations",
      HTML(
        '<p style="text-align: justify;">
    Coordinate file, data frame with 5 columns: 
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
<p style="text-align: justify;">
  Fam file examples for each format option, following the PLINK binary ".fam" file convention.
</p>

<table style="border-collapse: collapse; width: 100%;" border="1">
  <thead>
    <tr>
      <th>Family ID</th><th>Individual ID</th><th>Father</th><th>Mother</th><th>Sex</th><th>Phenotype</th>
    </tr>
  </thead>
</table>

<h4>First column preformatted as <code>POP1_POP2</code></h4>
<table style="border-collapse: collapse; width: 100%;" border="1">
  <tbody>
    <tr><td>AFR_ASW</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
    <tr><td>AFR_ASW</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
    <tr><td>EUR_CEU</td><td>NA12341</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
    <tr><td>EUR_CEU</td><td>NA06984</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
    <tr><td>EAS_CHB</td><td>NA18532</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
    <tr><td>EAS_CHB</td><td>NA18561</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
  </tbody>
</table>

<h4>First column describes one single population</h4>
<table style="border-collapse: collapse; width: 100%;" border="1">
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

<h4>First column describes several populations</h4>
<table style="border-collapse: collapse; width: 100%;" border="1">
  <tbody>
    <tr><td>ASW</td><td>NA19916</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
    <tr><td>ASW</td><td>NA19703</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
    <tr><td>CEU</td><td>NA12341</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
    <tr><td>CEU</td><td>NA06984</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
    <tr><td>CHB</td><td>NA18532</td><td>0</td><td>0</td><td>2</td><td>-9</td></tr>
    <tr><td>CHB</td><td>NA18561</td><td>0</td><td>0</td><td>1</td><td>-9</td></tr>
  </tbody>
</table>
'),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  ########################################### end info files ##############################  
  # value of slider M as min slider M1
  observeEvent(input$sliderM1, {
    updateSliderInput(session, "sliderM2", min = input$sliderM1)
  })
  
  # Reactive value to track button state
  panel_visible <- reactiveVal(FALSE)
  
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
 # shinyjs::disable("decompress_file")
  shinyjs::disable("run_merge")
  shinyjs::disable("run_pca")
  shinyjs::disable("run_admix")
  shinyjs::disable("run_fst")
  
  ################################################################################  
  # create a random folder and subfolders
  random_folder <- reactiveVal()
  ################################################## aqui empieza en fragmento a corregir 
  # ====== Crear carpeta temporal al iniciar análisis ======
  observeEvent(input$start_analysis, {
    panel_visible(!panel_visible())
    outputOptions(output, "panel_visible", suspendWhenHidden = FALSE)
    
    temp_folder <- tempfile(pattern = "user_output_", tmpdir = tempdir())
    dir.create(temp_folder)
    folder_path(temp_folder)  # ✅ aquí se asigna correctamente al reactiveVal
    
    subfolders <- c("decompressed_files", "pca", "merged", "reference", "admx", "fst", "plots", "user")
    lapply(subfolders, function(subfolder) {
      dir.create(file.path(temp_folder, subfolder))
    })
    
    random_folder(list(
      main_folder = temp_folder,
      subfolders = subfolders,
      paths = setNames(
        lapply(subfolders, function(sub) file.path(temp_folder, sub)),
        subfolders
      )
    ))
    
    for (subfolder in subfolders) {
      path <- file.path(temp_folder, subfolder)
      dir.create(path)
      cat("✅ Creado:", path, "\n")
    }
  })
  
  
  # ====== Mostrar carpetas creadas ======
  output$folder_path <- renderText({
    rf <- random_folder()
    req(rf)
    paste("Temporal user folder:", rf$main_folder)
  })
  
  output$subfolder_paths <- renderText({
    rf <- random_folder()
    req(rf)
    paste("Subfolders created:\n", paste(file.path(rf$main_folder, rf$subfolders), collapse = "\n"))
  })
  
  # load and processe reference file if wwp is selected
  # ====== Referencia reactiva desde Dropbox ======
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
        withProgress(message = "Downloading reference file from Dropbox...", value = 0, {
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
      shiny::showNotification("Reference file ready.", duration = 5, type = "message")
      
      return(file.path(temp_dir, "Reference"))
    }
    
    return(NULL)
  })
  
  # ====== Mostrar info sobre archivo de referencia ======
  output$reference_path <- renderText({
    ref_file <- reference_file()
    if (is.null(ref_file)) {
      return("Please select Minor or Major reference files.")
    }
    paste("Reference file at:", ref_file)
  })
  
  output$wait_ref <- renderText({
    req(input$selector)
    if (input$selector == "wwp") {
      "Wait until reference and target files are decompressed"
    } else {
      NULL
    }
  })
  
  
  output$wait_merge <- renderText({
    req(input$run_merge)
    if (input$run_merge == TRUE) {
      paste("Wait until files are merged")
    }
  })
  
  output$info_ref <- renderText({
    req(input$selector)
    if (input$selector == "wwp") {
      "Reference file:\nTotal genotyping rate is 0.997.\n3000000 variants and 3258 people pass filters and QC."
    } else {
      NULL
    }
  })
  

  ###################### load decompress and rename user file   #########################

  observeEvent(input$upload_file, {
    req(input$files == "usr")
    rf <- random_folder()
    req(rf)
    folder_path <- rf$main_folder
    
    file_ext <- tools::file_ext(input$upload_file$datapath)
    
    if (file_ext == "zip") {
      withProgress(message = "Decompressing user file...", value = 0, {
        decompress_dir <- file.path(folder_path, "user")
        if (Sys.info()[["sysname"]] %in% c("Linux", "Darwin")) {
          system2("unzip", args = c("-o", input$upload_file$datapath, "-d", decompress_dir), wait = TRUE)
        } else {
          unzip(input$upload_file$datapath, exdir = decompress_dir)
        }
        incProgress(1, detail = "Finished decompressing.")
        
        macosx_path <- file.path(decompress_dir, "__MACOSX")
        if (dir.exists(macosx_path)) {
          unlink(macosx_path, recursive = TRUE, force = TRUE)
        }
      })
    } else if (file_ext == "gz") {
      system(paste("gunzip", input$upload_file$datapath))
    } else {
      stop("Unsupported file type.")
    }
    
    bim_path <- list.files(file.path(folder_path, "user"), pattern = "\\.bim$", full.names = TRUE)
    stopifnot(length(bim_path) > 0)
    base_name <- tools::file_path_sans_ext(basename(bim_path))
    merged_file_path <- file.path(folder_path, "user", base_name)
    
    bim <- fread(paste0(merged_file_path, ".bim"), header = FALSE)
    bim[[1]] <- as.numeric(bim[[1]])
    bim_auto <- bim[bim[[1]] %in% 1:22, ]
    
    snps_to_keep_file <- file.path(folder_path, "keep_snps.txt")
    fwrite(data.frame(SNP = bim_auto[[2]]), snps_to_keep_file, col.names = FALSE)
    
    filtered_output_path <- file.path(folder_path, "user", "filtered_tmp")
    system(paste0(
      "www/plink19 --bfile ", merged_file_path,
      " --extract ", snps_to_keep_file,
      " --allow-extra-chr --make-bed --out ", filtered_output_path
    ))
    
    if (!file.exists(paste0(filtered_output_path, ".bed"))) {
      stop("❌ Error: PLINK no generó archivos filtrados. ¿Está vacío keep_snps.txt?")
    }
    
    # Guarda path para usar luego
    processed_file_path(filtered_output_path)
    
    output$status <- renderText("✅  User file loaded. Select FAM file format.")
  })
  
  # load file with user metadata
  observeEvent(input$upload_coord, {
    req(input$files == "usr")
    rf <- random_folder()
    req(rf)
    folder_path <- rf$main_folder
    
    # Crear subcarpeta 'user' si no existe
    user_dir <- file.path(folder_path, "user")
    dir.create(user_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Leer el archivo subido
    coord_df <- fread(input$upload_coord$datapath, header = FALSE)
    
    # Añadir "_u" al final de los elementos de la primera y segunda columna
    coord_df[[1]] <- paste0(coord_df[[1]], "_u")
    coord_df[[2]] <- paste0(coord_df[[2]], "_u")
    
    # Guardar como 'user_coord.txt' en la carpeta user
    coord_path <- file.path(user_dir, "user_coord.txt")
    fwrite(coord_df, coord_path, sep = "\t", col.names = FALSE)
    
    output$status <- renderText("✅ Loaded user metadata")
  })
  
  observeEvent(input$fam_formated, {
    req(processed_file_path())
    req(input$fam_formated)
    
    fam_path <- paste0(processed_file_path(), ".fam")
    user.fam <- fread(fam_path, header = FALSE)
    colnames(user.fam) <- c("V1", "V2", "V3", "V4", "V5", "V6")
    
    if (input$fam_formated == "unic") {
      user.id <- paste0(input$popID, "_", input$popID) # USR_USR
      user.fam$V1 <- user.id
    } else if (input$fam_formated == "subpop") {
      user.id <- paste0(input$popID, "_")
      user.fam[[1]] <- paste0(user.id, user.fam[[1]]) # USR_"as.famV1"
    } else if (input$fam_formated == "pop1pop2") {
      #user.fam
    }
    
    write.table(user.fam, paste0(processed_file_path(), ".fam"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    file.rename(paste0(processed_file_path(), ".bed"), file.path(dirname(processed_file_path()), "USER.bed"))
    file.rename(paste0(processed_file_path(), ".bim"), file.path(dirname(processed_file_path()), "USER.bim"))
    file.rename(paste0(processed_file_path(), ".fam"), file.path(dirname(processed_file_path()), "USER.fam"))
    
    # Procesar nombres
    fam <- fread(file.path(dirname(processed_file_path()), "USER.fam"), header = FALSE)
    fam.1 <- fam[, 1:2]
    colnames(fam.1) <- c("IDI", "IDII")
    fam.2 <- fam.1[, 1]
    colnames(fam.2) <- "V1"
    fam2 <- fam.2 %>%
      separate(V1, into = c("POP1", "POP2"), sep = "_", extra = "merge")
    if (input$fam_formated != "unic") {
      fam2[[1]] <- paste0(fam2[[1]], "_u")
      fam2[[2]] <- paste0(fam2[[2]], "_u")
    } else {
      fam2
    }
    
    df2 <- as.data.frame(cbind(fam2, fam.2))
    
    output$pop1_select <- renderUI({
      selectInput("pop1", "Select Population(s): (PCA & Admixture analysis)", choices = unique(df2$POP1), multiple = TRUE)
    })
    
    observe({
      req(df2)
      updateSelectInput(session, "pop1", selected = unique(df2$POP1)[1])
    })
    
    output$pop2_checkbox <- renderUI({
      req(input$pop1)
      pop2_choices <- df2[df2$POP1 %in% input$pop1, "POP2"]
      selectInput("pop2", "Select a pair of populations: ", choices = unique(pop2_choices), multiple = TRUE)
    })
    
    output$status <- renderText("✅ FAM file formatted and saved successfully.")
    
  })
  
################################################################################################################### 
################################################################################################################## 
 ####################### decompress and rename as INPUT example files ############################### 
  observeEvent(input$decompress_file, {
    
    rf <- random_folder()
    req(rf)
    folder_path <- rf$main_folder
    
    # Ensure only one progress object is created at the start
    progress <- shiny::Progress$new()
    progress$set(message = "Processing file...", value = 0)
    on.exit({
      progress$close()  # Ensure progress bar is closed at the end of processing
    }, add = TRUE)
    
  #  if (input$selector == "usr") {
      shinyjs::enable("run_pca")
  #  }
    
    if (input$files == "exp1" || input$files == "exp2"){           
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
      
        withProgress(message = "Decompressing file...", value = 0, {
          if (Sys.info()[["sysname"]] %in% c("Linux", "Darwin")) {
            system2("unzip", args = c("-o", ex_file_path, "-d", file.path(folder_path, "decompressed_files")), wait = TRUE)
          } else {
            unzip(ex_file_path, exdir = file.path(folder_path, "decompressed_files"))
          }
          
          incProgress(1, detail = "Finished decompressing.")
          
          # Eliminar carpeta __MACOSX si existe
          macosx_path <- file.path(folder_path, "decompressed_files", "__MACOSX")
          if (dir.exists(macosx_path)) {
            unlink(macosx_path, recursive = TRUE, force = TRUE)
          }
        })
        
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
          pop2_choices <- df2[df2$POP1 %in% input$pop1, "POP2"]
          selectInput("pop2", "Select a pair of populations: ", choices = unique(pop2_choices), multiple = TRUE)
        })
        
        # Update status to indicate file processed successfully
        output$status <- renderText("File processed successfully.")
    }
    })
 
  #####################################################################################  
  ##################################### MERGE STEP ####################################
  ################################# if wwp is selected ################################ 
  observeEvent(input$run_merge, {
    rf <- random_folder()
    req(rf)
    
    folder_path <- rf$main_folder
    path_decompressed <- rf$paths$decompressed_files
    path_user <- rf$paths$user
    path_merged <- rf$paths$merged
    ref_file_path <- file.path(folder_path, "reference")  # ✅ Asegura path a Reference
    
    if (input$files %in% c("exp1", "exp2")) {
      bfile_prefix <- file.path(path_decompressed, "INPUT")
      
      # 🧠 Si los archivos INPUT.* no existen, descomprime automáticamente
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
      
      # 🧠 Si los archivos USER.* no existen, el usuario no ha subido archivo
      if (!file.exists(paste0(bfile_prefix, ".bed"))) {
        showNotification("❌ Debes cargar y procesar el archivo del usuario antes de hacer el merge.", type = "error")
        return(NULL)
      }
    }
    
  
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
  
  # extract from reference those comon to user
  system(paste0("www/plink19 --bfile ",ref_file_path,"/Reference --extract ",bfile_prefix,".bim --aec --make-bed  --out ",path_merged,"/to_merge1"))
  
  # extract from user those comon to reference-extracted 
  system(paste0("www/plink19 --bfile ",bfile_prefix," --extract ",path_merged,"/to_merge1.bim --aec --make-bed  --out ",path_merged,"/to_merge2"))
  
  # merge step 1
  system(paste0("www/plink19 --bfile ",path_merged,"/to_merge1 --bmerge ",path_merged,"/to_merge2 --aec --make-bed --out ",path_merged,"/MERGE1"))
  
  # Check if ".missnp" file exist
  # Set the file path
  miss_file_path <- paste0(path_merged,"/MERGE1-merge.missnp")
  
  # Check if the file exists
  if (file.exists(miss_file_path)) { # if exist procede as:
    # merge step 1
    system(paste0("www/plink19 --bfile  ",path_merged,"/to_merge1 --bmerge ",path_merged,"/to_merge2 --aec --make-bed  --out ",path_merged,"/MERGE1"))
    # exclude step
    system(paste0("www/plink19 --bfile  ",path_merged,"/to_merge1 --exclude  ",path_merged,"/MERGE1-merge.missnp --aec --make-bed  --out ",path_merged,"/preMERGE_a"))
    system(paste0("www/plink19 --bfile  ",path_merged,"/to_merge2 --exclude  ",path_merged,"/MERGE1-merge.missnp --aec --make-bed  --out ",path_merged,"/preMERGE_b"))
    system(paste0("www/plink19 --bfile  ",path_merged,"/preMERGE_a --bmerge  ",path_merged,"/preMERGE_b --aec --make-bed  --out ",path_merged,"/MERGED"))
    
    # merge step 2 
  } else {                      # if not exist procede as:
    # merge step 2    
    system(paste0("www/plink19 --bfile  ",path_merged,"/to_merge1 --bmerge  ",path_merged,"/to_merge2 --aec --make-bed --out ",path_merged,"/MERGED"))
  }
  
  # activar processed_file_path
  processed_file_path(file.path(path_merged, "MERGED"))
  
  # recover codes from first column
  df2_reactive <- reactive({
    fam <- fread(paste0(path_merged,"/MERGED.fam"), header = F)
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
    pop2_choices <- df2[df2$POP1 %in% input$pop1, "POP2"]
    selectInput("pop2", "Select a pair of populations: ", choices = unique(pop2_choices), multiple = TRUE)
  })
  })
  
  ###########################################################################################
  ######################################## PCA STEP ###########################################
  ###########################################################################################
  # Variable reactiva global
  pca_files <- reactiveVal(NULL)
  PCA_data <- reactiveVal(NULL)
  
  observeEvent(input$run_pca, {
    req(processed_file_path())
    req(input$pop1)
    
    cat("📁 folder_path():", folder_path(), "\n")
    cat("📄 processed_file_path():", processed_file_path(), "\n")
    cat("🌍 input$pop1:", input$pop1, "\n")
    
    print("🔄 Actualizando `PCA_data()` con nuevos valores...")
    
    progress <- shiny::Progress$new()
    progress$set(message = "Running PCA by PLINK...", value = 0)
    on.exit({ progress$close() }, add = TRUE)
    
    update_progress <- function(amount) progress$inc(amount)
    
    rf <- random_folder()
    folder_path <- rf$main_folder
    pca_dir <- file.path(folder_path, "pca")
    pca_files(pca_dir)  # Guardamos path en variable reactiva
    
    path_example <- rf$paths$decompressed_files
    path_user <- rf$paths$user
    path_merged <- rf$paths$merged
    
    if (input$files %in% c("exp1", "exp2") && input$selector != 'wwp') {
      bfile_path <- file.path(path_example, "INPUT")
      required_files <- file.path(path_example, paste0("INPUT", c(".bed", ".bim", ".fam")))
    } else if (input$files == "usr" && input$selector != 'wwp') {
      bfile_path <- file.path(path_user, "USER")
      required_files <- file.path(path_user, paste0("USER", c(".bed", ".bim", ".fam")))
    } else {
      bfile_path <- file.path(path_merged, "MERGED")
      required_files <- file.path(path_merged, paste0("MERGED", c(".bed", ".bim", ".fam")))
    }
    
    if (!all(file.exists(required_files))) {
      showNotification("Missing required PLINK files (.bed, .bim, .fam).", type = "error", duration = 10)
      return(NULL)
    }
    
    pop.pre1 <- fread(paste0(bfile_path, ".fam"), header = FALSE)
    pop.pre1.1 <- pop.pre1[, 1:2]
    colnames(pop.pre1.1) <- c("IDI", "IDII")
    pop.pre2 <- pop.pre1[, 1]
    colnames(pop.pre2) <- "V1"
    pop2 <- pop.pre2 %>%
      separate(V1, into = c("POP1", "POP2"), sep = "_") %>%
      cbind(pop.pre1.1)
    
    cat("Estado actual ➤ selector:", input$selector, ", fam_formated:", input$fam_formated, "\n")
    
    if (input$selector == "usr") {
      if (!is.null(input$fam_formated) && input$fam_formated != "unic") {
        pop2[[1]] <- paste0(pop2[[1]], "_u")
        pop2[[2]] <- paste0(pop2[[2]], "_u")
      }
    }
    
    write.table(pop2, file.path(pca_dir, "pop_data.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    to_keep <- subset(pop2, POP1 %in% input$pop1)
    write.table(to_keep[, 3:4], file.path(pca_dir, "to_keep.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(to_keep[, 1:2], file.path(pca_dir, "to_color.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    tryCatch({
      if (input$select_all_pop1 == "Select all") {
        system(paste0("www/plink19 --bfile ", bfile_path, " --make-bed --out ", file.path(pca_dir, "step1")))
      } else {
        system(paste0("www/plink19 --bfile ", bfile_path, " --keep ", file.path(pca_dir, "to_keep.txt"),
                      " --make-bed --out ", file.path(pca_dir, "step1")))
      }
      
      update_progress(0.20)
      system(paste0("www/plink19 --bfile ", file.path(pca_dir, "step1"),
                    " --maf ", input$sliderM1, " --make-bed --out ", file.path(pca_dir, "step2")))
      
      update_progress(0.20)
      system(paste0("www/plink19 --bfile ", file.path(pca_dir, "step2"),
                    " --indep-pairwise 200 25 0.3 --out ", file.path(pca_dir, "prune1")))
      
      update_progress(0.20)
      system(paste0("www/plink19 --bfile ", file.path(pca_dir, "step2"),
                    " --extract ", file.path(pca_dir, "prune1.prune.in"), " --make-bed --out ", file.path(pca_dir, "step3")))
      
      update_progress(0.20)
      system(paste0("www/plink19 --bfile ", file.path(pca_dir, "step3"),
                    " --pca 10 header tabs var-wts --out ", file.path(pca_dir, "input_PCA")))
      
      update_progress(0.20)
      output$status <- renderText("PCA analysis completed successfully.")
    }, error = function(e) {
      showNotification("PLINK analysis failed. Check the logs.", type = "error", duration = 10)
      output$status <- renderText(paste("Error:", e$message))
    })
    
    PCA_results <- fread(file.path(pca_dir, "input_PCA.eigenvec"), header = TRUE)
    # Output the CV error lines to the verbatimTextOutput
    
    to_color <- fread(file.path(pca_dir, "to_color.txt"), header = FALSE)
    colnames(to_color) <- c("POP1", "POP2")
    
    PCA_to_plot <- cbind(to_color, PCA_results %>% select(PC1, PC2))
    PCA_data(PCA_to_plot)
    
    print("✅ `PCA_data()` se ha actualizado correctamente.")
    shinyjs::enable("run_admix")
    shinyjs::enable("run_fst")
    showNotification("PCA analysis finished. You can now run ADMIXTURE or FST analysis.", type = "message", duration = 10)
  })
  
  output$pca_log <- renderPrint({
    # Verifica que folder_path esté definido y válido
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
    matched <- grep("people pass filters and QC", log_lines, value = TRUE)
    
    if (length(matched) == 0) {
      cat("⚠️ Line not found in log.")
    } else {
      cat(matched)
    }
  })
  
  
  
  output$PCA_plot <- renderPlot({
    req(PCA_data())
    req(pca_files())
    
    PCA_to_plot <- PCA_data()
    pca_dir <- pca_files()
    
    message("✅ `renderPlot` está generando el gráfico correctamente.")
    print(head(PCA_to_plot))
    
    to_color <- fread(file.path(pca_dir, "to_color.txt"), header = FALSE)
    colnames(to_color) <- c("POP1", "POP2")
    
    to_color$POP1 <- factor(to_color$POP1, levels = unique(to_color$POP1))
    to_color$POP2 <- factor(to_color$POP2, levels = unique(to_color$POP2))
    
    PCA_to_plot <- cbind(to_color, PCA_to_plot %>% select(PC1, PC2))
    
    base_palette <- c(
      "red", "blue", "orange", "green", "yellow", "purple", "coral", "cyan", "salmon", "skyblue", "tomato",
      "magenta", "sandybrown", "gray", "chocolate"
    )
    
    selected_palette <- if (length(unique(PCA_to_plot$POP1)) > 15 || length(unique(PCA_to_plot$POP2)) > 15) {
      colorRampPalette(base_palette)(74)
    } else {
      base_palette
    }
    
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
        guides(color = guide_legend(title = "Subpopulation"), shape = guide_legend(title = "Population")) +
        scale_color_manual(values = palette_named) +
        theme_minimal(base_size = 14)
    }
    
    message("Colores usados:")
    print(palette_named)
    
    # Guardar imagen
    plot_dir <- file.path(dirname(pca_dir), "plots")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(plot_dir, "PCA.png"), plot = PCA)
    
    return(PCA)
  })
  
  ##########################################################################################
  ################################## FST STEP ##############################################
  ##########################################################################################
  
  # Variables reactivas globales para FST
  folder_path <- reactiveVal(NULL)
  fst_out_files <- reactiveVal(NULL)
  fst_trigger <- reactiveVal(0)
  
  # Ejecutar análisis FST
  observeEvent(input$run_fst, {
    req(input$run_pca)
    req(folder_path())
    
    cat("📁 folder_path():", folder_path(), "\n")
    cat("📄 Esperado pop file:", file.path(folder_path(), "pca", "pop_data.txt"), "\n")
    cat("📄 Esperado salida:", file.path(folder_path(), "fst", "FST_output.fst"), "\n")
    
    progress <- shiny::Progress$new()
    progress$set(message = "Running FST analysis", value = 0)
    on.exit(progress$close(), add = TRUE)
    
    pca_files <- file.path(folder_path(), "pca")
    b_files <- file.path(pca_files, "step3")
    fst_files <- file.path(folder_path(), "fst")
    fst_in_path <- file.path(fst_files, "to_fst.txt")
    fst_out_path <- file.path(fst_files, "FST_output")
    
    pop2 <- fread(file.path(pca_files, "pop_data.txt"), header = TRUE)
    to_fst <- subset(pop2, POP2 %in% input$pop2)[, c(3,4,2)]
    
    write.table(to_fst, fst_in_path, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cmd <- paste("www/gcta64 --bfile", b_files,
                 "--fst --sub-popu", fst_in_path,
                 "--out", fst_out_path)
    cat("📣 Ejecutando FST:", cmd, "\n")
    system(cmd)
    
    
    fst_out_files(fst_out_path)
    fst_trigger(fst_trigger() + 1)
  })
  
  # Cargar resultados FST
  fst_data <- reactive({
    req(fst_trigger())
    req(fst_out_files())
    
    df <- fread(paste0(fst_out_files(), ".fst"), header = TRUE)
    df <- df[complete.cases(df[, c("Chr", "bp", "Fst")]), ]
    df <- df[is.finite(df$Fst), ]
    df$Chr <- as.numeric(as.character(df$Chr))
    df$bp <- as.numeric(as.character(df$bp))
    df$Fst <- as.numeric(as.character(df$Fst))
    df <- df[order(df$Chr, df$bp), ]
    
    chr_lengths <- df %>%
      group_by(Chr) %>%
      summarise(chr_len = max(bp)) %>%
      arrange(Chr) %>%
      mutate(offset = dplyr::lag(cumsum(chr_len), default = 0))
    
    df <- df %>%
      left_join(chr_lengths, by = "Chr") %>%
      mutate(pos_cum = bp + offset)
    
    df$dbSNP_link <- paste0("https://www.ncbi.nlm.nih.gov/snp/", df$SNP)
    
    df
  })
  
  # Plot Manhattan interactivo
  output$manhattan_plot <- renderPlotly({
    req(fst_data())
    fst_in_path <- file.path(folder_path(), "fst", "to_fst.txt")
    df <- fst_data()
    
    pop_names <- fread(fst_in_path, header = FALSE)
    names_pop <- unique(pop_names$V3)
    
    df$color_group <- as.factor(df$Chr %% 2)
    
    chr_centers <- df %>%
      group_by(Chr) %>%
      summarise(center = mean(pos_cum)) %>%
      arrange(Chr)
    
    p <- plot_ly(
      data = df,
      x = ~pos_cum,
      y = ~Fst,
      type = 'scatter',
      mode = 'markers',
      color = ~color_group,
      colors = c("black", "red"),
      text = ~paste0(
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
        title = paste0("Interactive Manhattan plot of Fst-scores of population: ", names_pop),
        xaxis = list(
          title = "Chromosome",
          tickmode = "array",
          tickvals = chr_centers$center,
          ticktext = as.character(chr_centers$Chr)
        ),
        yaxis = list(title = "Fst"),
        showlegend = FALSE
      )
    
    temp_dir <- file.path(folder_path(), "plots")
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    htmlwidgets::saveWidget(p, file.path(temp_dir, "FST_interactive.html"), selfcontained = TRUE)
    
    p
  })
  
  # Tabla top 5 FST
  output$top5_fst <- renderPrint({
    req(fst_data())
    table <- fst_data() %>% dplyr::select(1, 2, 3, 4, 5, 6, 7,11)
    head(table %>% arrange(desc(Fst)), 5)
  })
  
  # Enlaces UCSC y dbSNP tras clic en SNP
  observeEvent(event_data("plotly_click"), {
    ed <- event_data("plotly_click")
    df <- fst_data()
    if (!is.null(ed)) {
      nearest_row <- df[which.min(abs(df$pos_cum - ed$x)), ]
      chr <- nearest_row$Chr
      pos <- nearest_row$bp
      start_ucsc <- max(0, pos - 10000)
      end_ucsc <- pos + 10000
      
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
  #################################### ADMIXTURE STEP  #####################################
  ##########################################################################################
 
  ########################################
  observeEvent(input$run_admix, {
    req(input$run_pca)  # Asegura que PCA se ha ejecutado primero
    
    progress <- shiny::Progress$new()
    progress$set(message = "Running ADMIXTURE analysis. This may take a while...", value = 0)
    on.exit(progress$close(), add = TRUE)
    
    update_progress <- function(amount) {
      progress$inc(amount)
    }
    
    output$status <- renderText("Running ADMIXTURE analysis...")
    
    admx_files <- file.path(folder_path(), "admx")
    pca_files <- file.path(folder_path(), "pca")
   # admixture_bin <- file.path(getwd(), "www", "admixture")
   # admixture_bin <- "/home/jfibla/miniconda3/bin/admixture"
    admixture_bin <- normalizePath("www/admixture")
    
    # Ajustar MAF
    m <- input$sliderM2
   # system(paste0("www/plink19 --bfile ", pca_files, "/step3 --maf ", m, " --mind 1.0 --make-bed --out ", admx_files, "/input_admx"))
    system(paste0("/srv/shiny-server/bpga_VM/www/plink19 --bfile ",pca_files,"/step3 --maf ", m, " --mind 1.0 --make-bed  --out ",admx_files,"/input_admx"))
    
    update_progress(0.05 / (input$sliderM2 - 0.05))
    
    
    
    
    # Mostrar salida del log MAF
    log_file <- file.path(admx_files, "input_admx.log")
    if (file.exists(log_file)) {
      keyword <- "variants and"
      grep_command <- paste("grep", shQuote(keyword), shQuote(log_file))
      log_lines <- tryCatch(system(grep_command, intern = TRUE), error = function(e) "No matching lines found.")
      output$maf_output <- renderText({
        paste0("✅ Number of SNPs and sample size:\n", paste(log_lines, collapse = "\n"))
      })
    } else {
      output$maf_output <- renderText({ "Log file not found. Ensure PLINK ran successfully." })
    }
    
    Np <- parallel::detectCores() - 1
    
  #  # Ejecutar ADMIXTURE (ya genera los .Q y log en admx_files directamente)
  #  if (input$Kselect == "one") {
  #    r <- input$sliderK
  #    system(paste0("cd ", admx_files, " && ", admixture_bin,
  #                  " --cv -j", Np, " input_admx.bed ", r, " | tee log", r, ".out"))
  #    update_progress(1.0 / (input$sliderK - 1))
  #  } else {
  #    for (i in 1:input$sliderK) {
  #      system(paste0("cd ", admx_files, " && ", admixture_bin,
  #                    " --cv -j", Np, " input_admx.bed ", i, " | tee log", i, ".out"))
  #      update_progress(1.0 / (input$sliderK - 1))
  #    }
  #  }
    ##################
   # admixture_bin <<- "/usr/local/bin/admixture_new"
   # admixture_bin <- "/home/jfibla/miniconda3/bin/admixture"
    admixture_bin <- normalizePath("www/admixture")
    
    if (input$Kselect == "one") {
      r <- input$sliderK
      system(paste0("cd ", admx_files, " && ", admixture_bin,
                    " --cv -j", Np, " input_admx.bed ", r, " | tee log", r, ".out"))
      update_progress(1.0 / (input$sliderK - 1))
    } else {
      for (i in 1:input$sliderK) {
        system(paste0("cd ", admx_files, " && ", admixture_bin,
                      " --cv -j", Np, " input_admx.bed ", i, " | tee log", i, ".out"))
        update_progress(1.0 / (input$sliderK - 1))
      }
    }
    #########
 
        # Leer todos los logs desde admx_files
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
    
    # Extraer líneas CV error
    cv_error_lines <- unlist(lapply(all_logs, function(log) grep("CV error", log, value = TRUE)))
    
    output$cvErrors <- renderPrint({
      if (length(cv_error_lines) > 0) {
        cat(cv_error_lines, sep = "\n")
      } else {
        cat("No CV error lines found.")
      }
    })
    
    # Actualizar input de cluster selector
    output$sorter <- renderUI({
      selectInput("sorter", "Sort admixture boxplot by 'Ancestry component': ", choices = 1:input$sliderK, selected = 1)
    })
    
    # Eliminar residuos si fuera necesario
    unlink(paste0(getwd(), "/", c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "log", "*.P", "log*", "*.out")), recursive = FALSE)

  
  ########################################
   # UI dinàmic per a la selecció del admixture component
    output$sorter <- renderUI({
      selectInput("sorter", "Sort admixture boxplot by 'Ancestry component': ", choices = 1:input$sliderK, multiple = F,selected=1)
    }) 
    ##############################################################################
    # Render the ADMIXTURE plot
    output$ADMIXTURE_plot <- renderPlot({
      # Load data
      runs <- fread(paste0(admx_files,"/input_admx.", input$sliderK, ".Q"))
      to_color <- fread(paste0(pca_files,"/to_color.txt"), header = F)
      colnames(to_color) <- c("POP1", "POP2")
      
      # Combine population and ancestry data
      admixture_data <- cbind(to_color, runs)
      data <- admixture_data
      # Paleta 10 colors per K
      k_palette <- c(
        "#4682B4", # steel blue"
        "#FFD700", # gold
        "#00FF00", # lime
        "#FF69B4", # hot pink
        "#556B2F", # dark olive green
        "#7FFFD4", # aquamarine
        "#9932CC", # dark orchid
        "#DEB887", # burlywood
        "#B22222",  # firebrick 
        "#008080" # teal
      )
      
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
        
        #sort POP1 alphabetic
        data_long$POP1 <- factor(data_long$POP1, levels = sort(unique(as.character(data_long$POP1))))
        
        # Crear el gràfic
        ADM <- ggplot(data_long, aes(x = Index, y = Value, fill = Variable)) +
          geom_bar(stat = "identity", position = "stack") +
          facet_nested(~POP1, scales = "free_x") +
          labs(x = "Individuals", y = "Ancestry component", title = "Admixture Plot",fill = "Components") +
          #  theme_minimal() +
          #  theme(legend.position = "bottom") +
          scale_fill_manual(values = k_palette) +  # Aplicar la paleta personalizada
          theme(
            plot.title = element_text(size = 18, face = "bold"),  # Title size and style
            axis.title.x = element_text(size = 16, face = "bold"),              # X-axis label size
            axis.title.y = element_text(size = 16, face = "bold"),              # Y-axis label size
            axis.text = element_text(size = 12),                 # Axis tick labels size
            legend.title = element_text(size = 14),              # Legend title size
            legend.text = element_text(size = 12)                # Legend text size
          )
        
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
        #sort POP2 alphabetic
        data_long$POP2 <- factor(data_long$POP2, levels = sort(unique(as.character(data_long$POP2))))
        
        ADM <- ggplot(data_long, aes(x = Index, y = Value, fill = Variable)) +
          geom_bar(stat = "identity", position = "stack") +
          facet_nested(~POP1 + POP2, scales = "free_x", space = "free") +
          labs(
            x = "Individuals",
            y = "Ancestry component",
            title = "Admixture Plot",
            fill = "Components"
          ) +
          scale_fill_manual(values = k_palette) +  # Aplicar la paleta personalizada
          theme(
            plot.margin = margin(10, 10, 10, 10), # Adjust as needed (top, right, bottom, left)
            panel.spacing = unit(0, "lines"), 
            strip.text.x = element_text(size = 14, face = "bold", angle = 0),  # Facet label size and style for columns
            strip.text.y = element_text(size = 14, face = "bold", angle = 45), # Facet label size and style for rows
            plot.title = element_text(size = 18, face = "bold"),  # Title size and style
            axis.title.x = element_text(size = 16, face = "bold"),              # X-axis label size
            axis.title.y = element_text(size = 16, face = "bold"),              # Y-axis label size
            axis.text = element_text(size = 12),                 # Axis tick labels size
            legend.title = element_text(size = 14),              # Legend title size
            legend.text = element_text(size = 12)                # Legend text size
          )
        
      }  
      print(ADM)
      temp_dir <- file.path(folder_path(), "plots")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
      
      # dir.create("/tmp/plots", showWarnings = FALSE, recursive = TRUE)
      # addResourcePath("plots", "/tmp/plots")
      # ploted_files <- "/tmp/plots/"
      # Save the plot as PNG
      ggsave(paste0(temp_dir, "/ADM.png"), plot = ADM)
    })
    ##############################################################################
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
      colnames(admixture_long) <-c("POP2","POP1","x","K","value")
      
      # Paleta  colores
      if (length(unique(admixture_long$POP2))<=15) {
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
            plot.title = element_text(size = 18, face = "bold"),  # Title size and style
            axis.title.x = element_text(size = 16, face = "bold"),              # X-axis label size
            axis.title.y = element_text(size = 16, face = "bold"),              # Y-axis label size
            axis.text = element_text(size = 12),                 # Axis tick labels size
            legend.title = element_text(size = 14),              # Legend title size
            legend.text = element_text(size = 12)                # Legend text size
          )
        
      } else {
        p <- ggplot(admixture_long, aes(x = K, y = value, fill = POP2)) + 
          geom_boxplot() +
          facet_nested(~POP1, scales = "free", space = "free") +
          # facet_grid(cols = vars(POP1), scales = "free", space = "free") +  # Manteniendo el facetado para POP1
          scale_fill_manual(values = selected_palete) +
          labs(x = "Ancestry component", y = "Ancestry Proportion", fill = "Population") +
          theme(
            plot.margin = margin(10, 10, 10, 10), # Adjust as needed (top, right, bottom, left)
            panel.spacing = unit(0, "lines"), 
            strip.text.x = element_text(size = 14, face = "bold", angle = 0),  # Facet label size and style for columns
            strip.text.y = element_text(size = 14, face = "bold", angle = 45), # Facet label size and style for rows 
            plot.title = element_text(size = 18, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
          )
      }    
      
      print(p)
      temp_dir <- file.path(folder_path(), "plots")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
      
      # dir.create("/tmp/plots", showWarnings = FALSE, recursive = TRUE)
      # addResourcePath("plots", "/tmp/plots")
      # ploted_files <- "/tmp/plots/"
      
      ggsave(paste0(temp_dir, "/box_ADM.png"), plot = p)
    })
    ################################################################################     
    ################################################################################     
    ################################################################################  
    map_with_save <- function() {
      req(input$sliderK)
      print("map_with_save definido correctamente")
      # define paths
      admx_files <- file.path(folder_path(), "admx")
      pca_files <- file.path(folder_path(), "pca")
      
      # read and format coordinate file
      coord <- fread("www/coord.txt", header = TRUE)

      if (!is.null(input$selector)) {
        if (input$selector %in% c("exp1", "exp2")) {
          coord
        } else if (input$selector == "usr") {
          if (!is.null(input$fam_formated) && input$fam_formated == "unic") {
            user_data <- data.frame(
              POP1 = input$popID,
              POP2 = input$popID,
              Pop_name = "user population",
              LAT = as.numeric(input$lat),
              LON = as.numeric(input$lon)
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
      
      # capture .Q runs
      runs <- fread(paste0(admx_files, "/input_admx.", input$sliderK, ".Q"))
      value_columns <- paste0("V", 1:input$sliderK)
      colnames(runs) <- value_columns
      
      # capture population data
      to_color <- fread(paste0(pca_files, "/to_color.txt"), header = FALSE)
      colnames(to_color) <- c("POP1", "POP2")
      
      # merge .Q and population data
      admix_data <- cbind(to_color, runs)
      
      # calculate .Q mean value by subpopulation (qpop) and superpopulation (QPOP)
      qpop <- aggregate(. ~ POP2, data = as.data.frame(admix_data)[, c("POP2", value_columns)], mean)
      QPOP <- aggregate(. ~ POP1, data = as.data.frame(admix_data)[, c("POP1", value_columns)], mean)
      print(qpop)
      print(QPOP)
      
      # calculate population POP2 counts to control cart size
      pop_counts <- to_color %>%
        count(POP2, name = "pop_size")
      pop_counts$pop_size <- as.numeric(pop_counts$pop_size)
      # merge file and prepare input for pie plot
      admixture_data <- merge(qpop, pop_counts, by = "POP2")
      coord_pop2_df <- merge(admixture_data, coord, by = "POP2")
      print(coord_pop2_df)
      
      # arrange file by POP1
      coord_pop1_df <- coord_pop2_df %>%
        group_by(POP1) %>%
        summarise(
          Pop_name = POP1,
          across(starts_with("V"), ~mean(.x, na.rm = TRUE), .names = "{.col}"),
          LAT = mean(LAT, na.rm = TRUE),
          LON = mean(LON, na.rm = TRUE),
          pop_size = sum(pop_size, na.rm = TRUE),
          .groups = "drop"
        )
      print(coord_pop1_df)
      
      
      # filter population to plot
       if (input$popcolor == "POP1") {
        coord_pop_df <- coord_pop1_df
      } else {
        coord_pop_df <- coord_pop2_df
      }
      
      df_sf <- coord_pop_df %>%
        st_as_sf(coords = c("LON", "LAT"), crs = 4326, remove = FALSE) %>%
        mutate(
          LON = st_coordinates(.)[, 1],
          LAT = st_coordinates(.)[, 2]
        )
      
      df_sf <- df_sf %>%
        mutate(scaled_size = scales::rescale(pop_size, to = c(40, input$sliderR)))
      
      k_palette <- c(
        "#4682B4", "#FFD700", "#00FF00", "#FF69B4", "#556B2F",
        "#7FFFD4", "#9932CC", "#DEB887", "#B22222", "#008080"
      )
      
      chartdata <- coord_pop_df %>%
        select(starts_with("V")) %>%
        as.matrix()
      
      num_pieces <- ncol(chartdata)
      dynamic_colors <- k_palette[1:num_pieces]
      
      leaflet_map <- leaflet(data = df_sf) %>%
        addTiles() %>%
        addMinicharts(
          lng = df_sf$LON,
          lat = df_sf$LAT,
          chartdata = chartdata,
          type = "pie",
          width = df_sf$scaled_size,
          colorPalette = dynamic_colors
        ) %>%
        addLabelOnlyMarkers(
          lng = df_sf$LON,
          lat = df_sf$LAT,
          label = lapply(df_sf$Pop_name, htmltools::HTML),
          labelOptions = labelOptions(noHide = TRUE, direction = "top", textOnly = TRUE, textsize = "12px")
        )
      
      temp_dir <- file.path(folder_path(), "plots")
      dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
      
      mapshot(
        leaflet_map,                # Tu objeto de mapa leaflet
        file = paste0(temp_dir,"/map_ADM.png"),  # Ruta donde guardar
        vwidth = 1650,   # Ancho del archivo PNG
        vheight = 1000   # Altura del archivo PNG
      )
      
      # return leaflet_map
      return(leaflet_map)
    }
    
    ################################################################################  
    output$map_plot <- renderLeaflet({
      map_with_save()
    })
    ################################################################################  
    
  }) # final input$admix
  
  observeEvent(input$lancaAdmixture, {
    req(input$fitxer, input$k)
    
  #  admixture_bin <- "/home/jfibla/miniconda3/bin/admixture"
    admixture_bin <- "www/admixture"
    
    fitxer <- paste0("data/", input$fitxer)
    k <- input$k
    log_path <- paste0("logs/", tools::file_path_sans_ext(input$fitxer), ".", k, ".log")
    
    # Executa la comanda i captura la sortida i l’error
    cmd <- paste(admixture_bin, fitxer, k)
    log_output <- tryCatch(
      {
        system(cmd, intern = TRUE)
      },
      error = function(e) {
        paste("Error d'execució:", e$message)
      }
    )
    
    # Desa-ho en un fitxer
    writeLines(log_output, log_path)
  })
  
  output$logAdmixture <- renderText({
    req(input$fitxer, input$k)
    log_path <- paste0("logs/", tools::file_path_sans_ext(input$fitxer), ".", input$k, ".log")
    if (file.exists(log_path)) {
      paste(readLines(log_path), collapse = "\n")
    } else {
      "Encara no s'ha generat cap log per aquest fitxer/K."
    }
  })
  ##############################################################################
  ############################## DOWNLOAD RESULTS STEP #########################
  ##############################################################################
  
  output$downloadPlotsZip <- downloadHandler(
    filename = function() {
      paste0("All_Plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Ruta temporal donde están los plots
      temp_dir <- file.path(folder_path(), "plots")
      
      # Verificamos que existan
      if (!dir.exists(temp_dir)) {
        showNotification("No hay plots para descargar.", type = "error")
        stop("Directorio de plots no encontrado.")
      }
      
      # Creamos archivo .zip
      zipfile <- tempfile(fileext = ".zip")
      zip(zipfile, files = list.files(temp_dir, full.names = TRUE), flags = "-r9Xj")
      
      # Copiamos el zip al archivo destino de descarga
      file.copy(zipfile, file)
    }
  )
  
  
  
  ################## to see temp_dir content at herohu logs ######################
  
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
  
  
}
# Reset All button functionality
print(sapply(sort(sapply(ls(), function(x) object.size(get(x))), decreasing = TRUE), function(size) round(size / 1024^2, 2)))

shinyApp(ui = ui, server = server)