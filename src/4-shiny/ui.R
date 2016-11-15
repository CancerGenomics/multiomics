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
          tags$fieldset(class="form-group shiny-input-container input-group",
                        tags$legend(h4("clinical data")),
                        fileInput("mirna.survivalFile", label="File"),
                        hidden(selectInput("mirna.survival.column.name","Survival column name",NULL)),
                        hidden(selectInput("mirna.event.column.name","Event column name",NULL))
          ),
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
		      fileInput("cnv.survivalFile", label=h4("clinical data")),
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
                 fileInput("meth.survivalFile", label=h4("clinical data")),
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

            tags$br(),
            actionButton("connectToXenaHub", "Looks for TCGA cohorts"),
            tags$br(),tags$br(),tags$br(),
		        fluidPage(
              fluidRow(
                column(4,
                  selectInput("xenaCohorts","XenaHub available cohorts",NULL, size = 30, selectize = F)
                ),
                column(4,
                  hidden(selectInput("xenaCohortDatasets","Selected cohort datasets",NULL, size = 30, selectize = F))
                )
              ),
              uiOutput("downloadLinkOutput")
            ) 
          )
        )
      )
    )
    )
  )
)