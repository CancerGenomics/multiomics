shinyUI(
  
  #useShinyjs(),  
  
  #fluidPage(
  navbarPage("multiOmics",
    tabPanel("mRNA-miRNA pipeline",    
			  
      sidebarLayout(
      
        sidebarPanel(
          #titlePanel("multiOmics"),
          tags$fieldset(
          
            tags$legend("mRNA profile"),
            
            fileInput("mrnaFile", accept=c("text/csv"), label="select local file"),
          
            actionButton("mRNAXenaLookup", "Look up mRNA"),
	  	      bsModal("mrnaXenaSelector", "Xena connector", "", size = "large",
		          selectInput("mRNACohorts","XenaHub available cohorts",NULL),
		          selectInput("mRNACohortDatasets","Cohort datasets",NULL),
				      #textOutput("textnew"),
				      actionButton("mRNAXenaRDownload", "Download selected dataset"),
				      downloadButton('downloadData', 'Download'),
				      actionButton("mRNAXenaRCancel", "Cancel")
            ),
	  	      tags$br(), tags$br(), tags$br(), tags$br()

          ),
          tags$fieldset(
            
            tags$legend("miRNA profile"),

            fileInput("mirnaFile",accept=c("text/csv"), label="select local file"),
            sliderInput("thresholdSlider", label=h4("Correlation coefficient"), 
				      min=0, max=1, value=0.7, step=0.05),
		        radioButtons("pearsons.method", label = h4("Correlation test"),
				       choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
				       selected = "spearman"),		
		        # directoryInput('outputDir', label = 'Select a directory for the output', value = '~'),
            # textInput('outputFileName', label = 'Output file name', value = 'inputStep2-matureMirnaXmrna.csv'),
            actionButton("runMRNAMiRNACorrelation", "Run pipeline")
          )
        ),
        mainPanel(
          tags$fieldset("asas"),
          
          DT::dataTableOutput('result'),
          plotOutput('correlationPlot')
        )
      )
    ),
    tabPanel("mRNA-CNV pipeline",    
      sidebarLayout(
					
        sidebarPanel(
		  #titlePanel("multiOmics"),
          fileInput("cnv.mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
		  fileInput("cnv.cnvFile",accept=c("text/csv"), label=h4("CNV profile")),
          sliderInput("cnv.thresholdSlider", label=h4("Correlation coefficient"), 
                      min=0, max=1, value=0.7, step=0.05),
          radioButtons("cnv.pearsons.method", label = h4("Correlation test"),
                       choices = c("Pearson" = "pearson", "Spearman" = "spearman", "Kendall" = "kendall"), 
                       selected = "spearman"),		
          actionButton("runMRNACNVCorrelation", "Run pipeline")
        ),
        mainPanel(
          DT::dataTableOutput('MRNACNVResult')
        )
      )
    )
  )
)