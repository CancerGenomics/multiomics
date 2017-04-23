shinyUI(

  fluidPage(

  useShinyjs(),
    
  titlePanel("multiOmics"),
  
  tabsetPanel(type = "pills", 

      tabPanel("mRNA-miRNA pipeline",    
			  
      sidebarLayout(
      
        sidebarPanel(
          fileInput("mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
          fileInput("mirnaFile",accept=c("text/csv"), label=h4("miRNA profile")),
          fileInput("mirna.survivalFile", label=h4("Follow-up data")),
          hidden(selectInput("mirna.survival.column.name","Survival column name",NULL)),
          hidden(selectInput("mirna.event.column.name","Event column name",NULL)),
          sliderInput("thresholdSlider", label=h4("Correlation coefficient"), 
				      min=0.3, max=1, value=0.7, step=0.05),
          radioButtons("pearsons.method", label = h4("Correlation test"),
				       choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
				       selected = "pearson"),		
		      h4(id="asas2","multiMiR"),
	        checkboxInput("miRNA.runMultimir",p(id="multimirTooltipText","miRNA-target interaction")),
		      bsTooltip("multimirTooltipText",
		                "Recover miRNA-mRNA target interactions through the multiMiR reource that includes 11 validated/predicted miRNA–target databases (e.g.: miRecords, miRTar-Base, miRanda, etc.)", 
		                placement = "right",options = list(container = "body")),
		      tags$hr(),
          actionButton("runMRNAMiRNACorrelation", "Run pipeline")
        ),
        mainPanel(
          DT::dataTableOutput('result'),
          tags$div(id="downloadMrnaMirnaResultDiv",
              shinyjs::hidden(downloadButton("downloadMrnaMirnaResult", "Download csv"))),

          plotOutput('correlationPlot'),
          plotOutput('correlationSurvival')
        )
      )
    ),
    tabPanel("mRNA-CNV pipeline",    
      sidebarLayout(
					
        sidebarPanel(
          fileInput("cnv.mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
          fileInput("cnv.cnvFile",accept=c("text/csv"), label=h4("CNV profile")),
          fileInput("cnv.survivalFile", label=h4("Follow-up data")),
          hidden(selectInput("cnv.survival.column.name","Survival column name",NULL)),
          hidden(selectInput("cnv.event.column.name","Event column name",NULL)),
          sliderInput("cnv.thresholdSlider", label=h4("Correlation coefficient"), 
                      min=0, max=1, value=0.7, step=0.05),
          radioButtons("cnv.pearsons.method", label = h4("Correlation test"),
                       choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
                       selected = "pearson"),		
          tags$hr(),
          actionButton("runMRNACNVCorrelation", "Run pipeline")
        ),
        mainPanel(
          DT::dataTableOutput('MRNACNVResult'),
          tags$div(id="downloadMrnaCNVResultDiv",
                   shinyjs::hidden(downloadButton("downloadMrnaCNVResult", "Download csv"))),
          
          plotOutput('cnv.correlationPlot'),
          plotOutput('cnv.correlationSurvival')
        )
      )
    ),
    tabPanel("mRNA-methylation pipeline",    
             sidebarLayout(
               
               sidebarPanel(
                 fileInput("meth.mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
                 fileInput("meth.methFile",accept=c("text/csv"), label=h4("Methylation profile")),
                 # TODO habria que tomar los choices de las plataformas disponibles en getMethylationPlatformNames()
                 selectInput("meth.platform.select", label = h4("Platform"), choices = c("HumanMethylation450 BeadChip")),
                 sliderInput("meth.thresholdSlider", label=h4("Correlation coefficient"), 
                             min=0.3, max=1, value=0.7, step=0.05),
                 radioButtons("meth.pearsons.method", label = h4("Correlation test"),
                              choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
                              selected = "pearson"),		
                 tags$hr(),
                 actionButton("runMRNAMethylationCorrelation", "Run pipeline")
               ),
               mainPanel(
                 DT::dataTableOutput('MRNAMethResult'),
                 tags$div(id="downloadMrnaMethylationResultDiv",
                          shinyjs::hidden(downloadButton("downloadMrnaMethResult", "Download csv"))),

                 plotOutput('meth.correlationPlot'),
                 plotOutput('meth.correlationSurvival')
                 
                 
               )
             )
    ),
    tabPanel("TCGA Xena Hub",
      
      tags$div(class="row ",
                
        tags$form(class="well",
          tags$div(class="form-group shiny-input-container",

            fluidPage(
              fluidRow(
                column(5,
                   actionButton("connectToXenaHub", "Looks for TCGA cohorts")
                ),
                column(4,
                       uiOutput("downloadLinkOutput")
                )
              ),
              fluidRow(column(12,tags$br())),
              fluidRow(
                column(5,
                  selectInput("xenaCohorts","XenaHub available cohorts",NULL, size = 30, selectize = F)
                ),
                column(6,
                  hidden(selectInput("xenaCohortDatasets","Selected cohort datasets",NULL, size = 30,selectize = F))
                )
              )
            ) 
          )
        )
      )
    )
    )
  )
)