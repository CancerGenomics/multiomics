shinyUI(

  fluidPage(

  useShinyjs(),
    
  titlePanel("multiOmics"),
  
  tabsetPanel(type = "pills", 

      tabPanel("mRNA-miRNA pipeline",    
			  
      sidebarLayout(
      
        sidebarPanel(
          #titlePanel("multiOmics"),
          fileInput("mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
          fileInput("mirnaFile",accept=c("text/csv"), label=h4("miRNA profile")),
          sliderInput("thresholdSlider", label=h4("Correlation coefficient"), 
				      min=0.3, max=1, value=0.7, step=0.05),
		      radioButtons("pearsons.method", label = h4("Correlation test"),
				       choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
				       selected = "pearson"),		
		        # directoryInput('outputDir', label = 'Select a directory for the output', value = '~'),
            # textInput('outputFileName', label = 'Output file name', value = 'inputStep2-matureMirnaXmrna.csv'),
          actionButton("runMRNAMiRNACorrelation", "Run pipeline")
        ),
        mainPanel(
          DT::dataTableOutput('result'),
          tags$div(id="downloadMrnaMirnaResultDiv",
              shinyjs::hidden(downloadButton("downloadMrnaMirnaResult", "Download csv"))),

          plotOutput('correlationPlot')
        )
      )
    ),
    tabPanel("mRNA-CNV pipeline",    
      sidebarLayout(
					
        sidebarPanel(
          fileInput("cnv.mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
		      fileInput("cnv.cnvFile",accept=c("text/csv"), label=h4("CNV profile")),
          sliderInput("cnv.thresholdSlider", label=h4("Correlation coefficient"), 
                      min=0.3, max=1, value=0.7, step=0.05),
          radioButtons("cnv.pearsons.method", label = h4("Correlation test"),
                       choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
                       selected = "pearson"),		
          actionButton("runMRNACNVCorrelation", "Run pipeline")
        ),
        mainPanel(
          DT::dataTableOutput('MRNACNVResult'),
          tags$div(id="downloadMrnaCNVResultDiv",
                   shinyjs::hidden(downloadButton("downloadMrnaCNVResult", "Download csv")))
          
        )
      )
    ),
    tabPanel("mRNA-methylation pipeline",    
             sidebarLayout(
               
               sidebarPanel(
                 fileInput("meth.mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
                 fileInput("meth.methFile",accept=c("text/csv"), label=h4("Methylation profile")),
                 sliderInput("meth.thresholdSlider", label=h4("Correlation coefficient"), 
                             min=0.3, max=1, value=0.7, step=0.05),
                 radioButtons("meth.pearsons.method", label = h4("Correlation test"),
                              choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
                              selected = "pearson"),		
                 actionButton("runMRNAMethylationCorrelation", "Run pipeline")
               ),
               mainPanel(
                 DT::dataTableOutput('MRNAMethylationResult'),
                 tags$div(id="downloadMrnaMethylationResultDiv",
                          shinyjs::hidden(downloadButton("downloadMrnaMethylationResult", "Download csv")))
                 
               )
             )
    ),
    tabPanel("TCGA Xena Hub",
      
      tags$div(class="row ",
                 tags$form(class="well",
                           tags$div(class="form-group shiny-input-container",
        #style = "background-color: #eeeeee; 
        #              border: 1px solid #dddddd;
        #border-radius: 5px;
        #margin: 20px" ,
        
      tags$br(),
      actionButton("connectToXenaHub", "Looks for TCGA cohorts"),
      tags$br(),tags$br(),tags$br(),
		  fluidPage(
       fluidRow(
          column(4,
            #bsModal("mrnaXenaSelector", "Xena connector", "", size = "large",
            selectInput("xenaCohorts","XenaHub available cohorts",NULL, size = 30, selectize = F)
          ),
          column(4,
            hidden(selectInput("xenaCohortDatasets","Selected cohort datasets",NULL, size = 30, selectize = F))
          )
        ),
        #downloadButton("downloadData", "Download selected dataset"),
        uiOutput("downloadLinkOutput")

        #onclick(id = "downloadData3", expr = "window.open(document.getElementById('link').value);")
        #shinyjs::runjs("window.open(document.getElementById('link').value);")
        
        #),
      ) 
    ))
    ))
  )
  )
)