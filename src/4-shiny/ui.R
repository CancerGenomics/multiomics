shinyUI(

  
  fluidPage(

  useShinyjs(),
  rclipboardSetup(),
    
  titlePanel("multiOmics"),
  
  # dialogo modal para mensajes generales          
  bsModal("generalMessageModal", title = "multiOmics message", 
          size="large", trigger = "",
          textOutput("generalMessageOutputText")
  ),
  
  tabsetPanel(type = "pills", 

      tabPanel("mRNA-miRNA pipeline",    
			  
      sidebarLayout(
      
        sidebarPanel(


          fileInput("mrnaFile", label=h4("mRNA profile",style="display: inline-block;",
                                                               actionLink(inputId="mrnaFileHelp", label="", 
                                                                          icon = icon("question-sign",lib = "glyphicon")),
                                                               bsModal("mrnaFileHelpModal", title = "mRNA file format", 
                                                                       trigger = "mrnaFileHelp",size="large",
                                                                       p("File must have this information structure, and columns must be delimited by the tab character"),
                                                                       img(src="mrna_format.png", width = "100%")
                                                                       )
                                         )
                    
                    ),
          tags$div(id="mirnaMrnaFileErrorMsgDiv",
                   shinyjs::hidden(
                     p(id="mirnaMrnaFileErrorMsg",
                       class="alert alert-danger",
                       "mRNA file has a wrong format, please verify the file and re-run the pipeline")
                                                               )
                    ),
          
          fileInput("mirnaFile", label=h4("miRNA profile",style="display: inline-block;",
                                                               actionLink(inputId="mirnaFileHelp", label="", 
                                                                          icon = icon("question-sign",lib = "glyphicon")),
                                                               bsModal("mirnaFileHelpModal", title = "miRNA file format", 
                                                                       trigger = "mirnaFileHelp",size="large",
                                                                       p("File must have this information structure, and columns must be delimited by the tab character"),
                                                                       img(src="mirna_format.png", width = "100%")
                                                                       )
                                                               )
                    ),
          tags$div(id="mirnaMirnaFileErrorMsgDiv",
                   shinyjs::hidden(
                     p(id="mirnaMirnaFileErrorMsg",
                       class="alert alert-danger",
                       "miRNA file has a wrong format, please verify the file and re-run the pipeline")
                   )
          ),
          fileInput("mirna.survivalFile", label=h4("Follow-up data",style="display: inline-block;",
                                                   actionLink(inputId="followUpFileHelp", label="", 
                                                              icon = icon("question-sign",lib = "glyphicon")),
                                                   bsModal("followUpFileHelpModal", title = "follow up data file format", 
                                                           trigger = "followUpFileHelp",size="large",
                                                           p("File must have this information structure, and columns must be delimited by the tab character"),
                                                           img(src="fllowup_format.png", width = "100%")
                                                           )
                                                   )
                    ),
          tags$div(id="mirnaSurvivalFileErrorMsgDiv",
                   shinyjs::hidden(
                     p(id="mirnaSurvivalFileErrorMsg",
                       class="alert alert-danger",
                       "Follow up file has a wrong format, please verify the file and re-run the pipeline")
                   )
          ),
          hidden(selectInput("mirna.survival.column.name","Survival column name",NULL)),
          hidden(selectInput("mirna.event.column.name","Event column name",NULL)),
          sliderInput("thresholdSlider", label=h4("Correlation coefficient"), 
				      min=0.3, max=1, value=0.7, step=0.05),
          radioButtons("mirna.correlation.type", label = h4("Correlation type"),
				       choices = c("Positive" = "positive", "Negative" = "negative", "Both" = "both"), 
				       selected = "negative"),		
		      h4(id="asas2","multiMiR"),
	        checkboxInput("miRNA.runMultimir",p(id="multimirTooltipText","miRNA-target interaction")),
		      bsTooltip("multimirTooltipText",
		                "Recover miRNA-mRNA target interactions through the multiMiR resource that includes 11 validated/predicted miRNA–target databases (e.g.: miRecords, miRTar-Base, miRanda, etc.)", 
		                placement = "right",options = list(container = "body")),
		      tags$hr(),
          actionButton("runMRNAMiRNACorrelation", "Run pipeline")
        ),
        mainPanel(
          DT::dataTableOutput('result'),
          tags$div(id="downloadMrnaMirnaResultDiv",
              shinyjs::hidden(downloadButton("downloadMrnaMirnaResult", "Download csv"))
          ),
          # UI ouputs for the copy-to-clipboard buttons
          shinyjs::hidden(uiOutput("mirnaClip")),
          
          plotOutput('correlationPlot'),
          plotOutput('correlationSurvival')
        )
      )
    ),
    tabPanel("mRNA-CNV pipeline",    
      sidebarLayout(
					
        sidebarPanel(
          fileInput("cnv.mrnaFile", label=h4("mRNA profile",style="display: inline-block;",
                                                                   actionLink(inputId="cnvMrnaFileHelp", label="", 
                                                                              icon = icon("question-sign",lib = "glyphicon")),
                                                                   bsModal("cnvMrnaFileHelpModal", title = "mRNA file format", 
                                                                           trigger = "cnvMrnaFileHelp",size="large",
                                                                           p("File must have this information structure, and columns must be delimited by the tab character"),
                                                                           img(src="mrna_format.png", width = "100%")
                                                                           )
                                                                   )
                    ),
          tags$div(id="cnvMrnaFileErrorMsgDiv",
                   shinyjs::hidden(
                     p(id="cnvMrnaFileErrorMsg",
                       class="alert alert-danger",
                       "mRNA file has a wrong format, please verify the file and re-run the pipeline")
                   )
          ),
          fileInput("cnv.cnvFile", label=h4("CNV profile",style="display: inline-block;",
                                                                 actionLink(inputId="cnvFileHelp", label="", 
                                                                            icon = icon("question-sign",lib = "glyphicon")),
                                                                 bsModal("cnvFileHelpModal", title = "cnv file format", 
                                                                         trigger = "cnvFileHelp",size="large",
                                                                         p("File must have this information structure, and columns must be delimited by the tab character"),
                                                                         img(src="cnv_format.png", width = "100%")
                                                                         )
                                                                 )
                    ),
          tags$div(id="cnvFileErrorMsgDiv",
                   shinyjs::hidden(
                     p(id="cnvFileErrorMsg",
                       class="alert alert-danger",
                       "CNV file has a wrong format, please verify the file and re-run the pipeline")
                   )
          ),
          fileInput("cnv.survivalFile", label=h4("Follow-up data",style="display: inline-block;",
                                                 actionLink(inputId="cnvFollowUpFileHelp", label="", 
                                                            icon = icon("question-sign",lib = "glyphicon")),
                                                 bsModal("cnvFollowUpFileHelpModal", title = "Follow-up file format", 
                                                         trigger = "cnvFollowUpFileHelp",size="large",
                                                         p("File must have this information structure, and columns must be delimited by the tab character"),
                                                         img(src="fllowup_format.png", width = "100%")
                                                         )
                                                 )
                    ),
          tags$div(id="cnvSurvivalFileErrorMsgDiv",
                   shinyjs::hidden(
                     p(id="cnvSurvivalFileErrorMsg",
                       class="alert alert-danger",
                       "Follow up file has a wrong format, please verify the file and re-run the pipeline")
                   )
          ),
          hidden(selectInput("cnv.survival.column.name","Survival column name",NULL)),
          hidden(selectInput("cnv.event.column.name","Event column name",NULL)),
          sliderInput("cnv.thresholdSlider", label=h4("Correlation coefficient"), 
                      min=0, max=1, value=0.7, step=0.05),
          radioButtons("cnv.correlation.type", label = h4("Correlation type"),
                       choices = c("Positive" = "positive", "Negative" = "negative", "Both" = "both"), 
                       selected = "negative"),				
          tags$hr(),
          actionButton("runMRNACNVCorrelation", "Run pipeline")
        ),
        mainPanel(
          DT::dataTableOutput('MRNACNVResult'),
          tags$div(id="downloadMrnaCNVResultDiv",
                   shinyjs::hidden(downloadButton("downloadMrnaCNVResult", "Download csv"))),
          
          shinyjs::hidden(uiOutput("cnvClip")),
          
          plotOutput('cnv.correlationPlot'),
          plotOutput('cnv.correlationSurvival')
        )
      )
    ),
    tabPanel("mRNA-methylation pipeline",    
             sidebarLayout(
               
               sidebarPanel(
                 fileInput("meth.mrnaFile", label=h4("mRNA profile",style="display: inline-block;",
                                                                           actionLink(inputId="methMrnaFileHelp", label="", 
                                                                                      icon = icon("question-sign",lib = "glyphicon")),
                                                                           bsModal("methMrnaFileHelpModal", title = "mRNA file format", 
                                                                                   trigger = "methMrnaFileHelp",size="large",
                                                                                   p("File must have this information structure, and columns must be delimited by the tab character"),
                                                                                   img(src="mrna_format.png", width = "100%")
                                                                                   )
                                                                           )
                           ),
                 tags$div(id="methMrnaFileErrorMsgDiv",
                          shinyjs::hidden(
                            p(id="methMrnaFileErrorMsg",
                              class="alert alert-danger",
                              "mRNA file has a wrong format, please verify the file and re-run the pipeline")
                          )
                 ),
                 fileInput("meth.methFile", label=h4("Methylation profile",style="display: inline-block;",
                                                                          actionLink(inputId="methFileHelp", label="", 
                                                                                     icon = icon("question-sign",lib = "glyphicon")),
                                                                          bsModal("methFileHelpModal", title = "Methylation file format", 
                                                                                  trigger = "methFileHelp",size="large",
                                                                                  p("File must have this information structure, and columns must be delimited by the tab character"),
                                                                                  img(src="meth_format.png", width = "100%")
                                                                                  )
                                                                          )
                           ),
                 tags$div(id="methFileErrorMsgDiv",
                          shinyjs::hidden(
                            p(id="methFileErrorMsg",
                              class="alert alert-danger",
                              "Methylation file has a wrong format, please verify the file and re-run the pipeline")
                          )
                 ),
                 # TODO habria que tomar los choices de las plataformas disponibles en getMethylationPlatformNames()
                 selectInput("meth.platform.select", label = h4("Platform"), choices = c("HumanMethylation450 BeadChip")),
                 sliderInput("meth.thresholdSlider", label=h4("Correlation coefficient"), 
                             min=0.3, max=1, value=0.7, step=0.05),
                 radioButtons("meth.correlation.type", label = h4("Correlation type"),
                              choices = c("Positive" = "positive", "Negative" = "negative", "Both" = "both"), 
                              selected = "positive"),	
                 tags$hr(),
                 actionButton("runMRNAMethylationCorrelation", "Run pipeline")
               ),
               mainPanel(
                 DT::dataTableOutput('MRNAMethResult'),
                 tags$div(id="downloadMrnaMethylationResultDiv",
                          shinyjs::hidden(downloadButton("downloadMrnaMethResult", "Download csv"))),

                 shinyjs::hidden(uiOutput("methClip")),
                 
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
                       hidden(selectInput("xenaCohortDatasetsFilter","Filter datasets",
                                          c("All","miRNA","RNA","CopyNumber","Methylation"),
                                          selected = "All")),
                       hidden(selectInput("xenaCohortDatasets","Selected cohort datasets",NULL, size = 26,selectize = F))
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