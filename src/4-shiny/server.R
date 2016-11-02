# max size is 500MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
  
  sharedValues <- reactiveValues(fromButton=F,correlations="",cnvMrnaCorrelations="")
  
  ###########################################################################
  ########################## MIRNA - MRNA PIPELINE TAB
  ###########################################################################
  
  observeEvent(input$runMRNAMiRNACorrelation, { 
    runMRNAMiRNACorrelation()
  })

  mrnaExpressionData <- reactive({
    print("mrnaExpressionData")
    readMrnaExpressionFile(input$mrnaFile$datapath)
  })
  
  mirnaExpressionData <- reactive({
    print("mirnaExpressionData")
    readMirnaExpressionFile(input$mirnaFile$datapath)
  })
  
  threshold <- reactive({
	input$thresholdSlider
  })

  pearsonsMethod <- reactive({
	input$pearsons.method
  })

 correlations <- reactive(quote({
    if(sharedValues$fromButton) {
	  sharedValues$correlations <- CalculateCorrelationsMirnaMrna(
			                       mrnaExpressionData(), mirnaExpressionData(),
		   						   output.path="~/", 
								   output.file.name = paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-outputFile.csv", sep = ""),
								   r.minimium = threshold(), inc.progress = T, 
								   pearsons.method = pearsonsMethod())
	} 
	return (sharedValues$correlations)
  }), quoted = T)
  
  runMRNAMiRNACorrelation <- function() { 
          withProgress(message = 'Please stand by...', 
          detail = "calculating correlation", 
          min=0, max=1, {
            
     if(!is.null(input$mrnaFile) && !is.null(input$mirnaFile)) {

       #Checks if both files has the same samples in the same order. If not, aborts the execution.
       print("Checking if both files has the same samples in the same order...")
       suppressWarnings(checkSamplesFormIRNArnaCorrelation(mrnaExpressionData(), mirnaExpressionData(), 1))
       print("Preparing...")
        
       sharedValues$fromButton <- T
       output$result <- DT::renderDataTable(correlations(), selection = 'single')
       #sharedValues$fromButton <- F
        
       if(nrow(correlations()) > 1) {
         print("Hay resultados")
       } else {
         print("NO hay resultados")
       }      
       sharedValues$fromButton <- F
    } else {
       print("No hay archivos cargados")
    }
            
      
  })
} 

output$correlationPlot <- renderPlot({
  if(!is.null(input$result_rows_selected)){
    selected.gene <- correlations()[input$result_rows_selected,1]
    selected.mirna <- correlations()[input$result_rows_selected,2]
    selected.gene.row <- which(mrnaExpressionData()==selected.gene)
    selected.mirna.row <- which(mirnaExpressionData()==selected.mirna)
    X <- as.numeric(as.vector(mrnaExpressionData()[selected.gene.row,2:ncol(mrnaExpressionData())]))
    Y <- as.numeric(as.vector(mirnaExpressionData()[selected.mirna.row,2:ncol(mirnaExpressionData())]))
    cor.test(X, Y)
    plot(X, Y, col='Black', pch=1) #col=Group
    line <- lm(Y ~ X)
    abline(line, col="blue")
  }
})


###########################################################################
########################## CNV - MRNA PIPELINE TAB
###########################################################################

  observeEvent(input$runMRNACNVCorrelation, { 
	runMRNACNVCorrelation()
   })

  cnvThreshold <- reactive({
    input$cnv.thresholdSlider
  })

  cnvMrnaExpressionData <- reactive({
	print("mrnaExpressionData")
	readMrnaExpressionFile(input$cnv.mrnaFile$datapath)
  })

  cnvExpressionData <- reactive({
	print("cnvExpressionData")
	readCNVFile(input$cnv.cnvFile$datapath)
  })

  cnvPearsonsMethod <- reactive({
    input$cnv.pearsons.method
  })
  
  cnvMrnaCorrelations <- reactive(quote({
    if(sharedValues$fromButton) {
		sharedValues$cnvMrnaCorrelations <- CnvXMrnas(cnvMrnaExpressionData(), cnvExpressionData(), output.path="~/", 
			                                        output.file.name=paste(input$cnv.mrnaFile$name,"-",input$cnv.cnvFile$name,"-outputFile.csv", sep = ""),
													r.minimium = cnvThreshold(), inc.progress = T, pearsons.method = cnvPearsonsMethod())
    }
    return (sharedValues$cnvMrnaCorrelations)
  }), quoted = T)


  runMRNACNVCorrelation <- function(){
	  withProgress(message = 'Please stand by...', 
			  detail = "calculating correlation", 
			  min=0, max=1, {
				  
				  if(!is.null(input$cnv.mrnaFile) && !is.null(input$cnv.cnvFile)) {
					  
					  print("Preparing...")
					  
					  sharedValues$fromButton <- T
					  output$MRNACNVResult <- DT::renderDataTable(cnvMrnaCorrelations(), selection = 'single')

					  if(nrow(cnvMrnaCorrelations()) > 1) {
						  print("Hay resultados")
					  } else {
						  print("NO hay resultados")
					  }      
					  sharedValues$fromButton <- F
				  } else {
					  print("No hay archivos cargados")
				  }
				  
    })
  }

  
  ###########################################################################
  ########################## XENA HUB CONNECTOR TAB
  ###########################################################################
  
  mrnaCohortsData <- reactive({
    withProgress(message = 'Connectig with Xena...', 
                 detail = "Getting cohorts", 
                 min=0, max=1, {    
                   print("Connecting to XenaHub...")
                   cohorts(XenaHub(hosts = "https://tcga.xenahubs.net"))
                   #cohorts(XenaHub())
                 })  
  })
  
  observeEvent(input$connectToXenaHub, { 
    updateSelectInput(session, "xenaCohorts","XenaHub available cohorts", choices = mrnaCohortsData())
    #toggleModal(session,"mrnaXenaSelector")
  })  
  
  observeEvent(input$xenaCohorts, {
    withProgress(message = 'Connectig with Xena...', 
                 detail = "Getting datasets", 
                 min=0, max=1, {    
                   if(input$xenaCohorts != "" ){
                     #print(input$xenaCohorts)
                     if(input$xenaCohorts != "(unassigned)"){
                       print(paste("Searching datesets for", input$xenaCohorts, sep=" "))
                       ds <- datasets(XenaHub(cohorts = input$xenaCohorts))
                       updateSelectInput(session, "xenaCohortDatasets","Cohort datasets", choices = ds)
                       shinyjs::show("xenaCohortDatasets")
                     } else {
                       ds <- c("")
                       updateSelectInput(session, "xenaCohortDatasets","Cohort datasets", choices = ds)
                       shinyjs::hide("xenaCohortDatasets")
                     }
                   }
                 })   
  })   
  
  observeEvent(input$xenaCohortDatasets, {  
    print("Update datasets")
    print(input$xenaCohortDatasets)
    temp <- getUrlFromTCGAXenaHub(input$xenaCohortDatasets)
    print(temp)
    output$url <- renderText(temp)
    updateTextInput(session, inputId = "link", value = temp)
    update
    #show("downloadLink")
  })
  
  output$downloadLinkOutput <- renderUI({
    if((!is.null(input$xenaCohortDatasets)) && (input$xenaCohortDatasets != "")){
      print(input$xenaCohortDatasets)
      tagList(
        tags$a(id="downloadLink",
               tags$i(class="fa fa-download"),
               "Download selected datasets",href=getUrlFromTCGAXenaHub(input$xenaCohortDatasets), target="_blank", 
               class="btn btn-default shiny-download-link"
               )
      )
    } else {
      tagList()
    }
  })  

#  downloadFile <- reactive({
#    withProgress(message = 'Connectig with Xena...', 
#                 detail = "Downloading file", 
#                 min=0, max=1, {    
#      if(input$xenaCohorts != "" && input$xenaCohortDatasets != ""){
#        print(paste("Download file for datesets for", input$xenaCohorts, input$xenaCohortDatasets, sep=" "))
#        url <- getUrlFromTCGAXenaHub(input$xenaCohortDatasets)
#        print(url)
#        download.file(url = url, destfile = tempfile())
#      }
#    })
#  })
  
#  output$downloadData <- downloadHandler(
#    filename = function() { 
#      paste(input$xenaCohortDatasets, ".csv", sep="")
#    },
#    content = function(file) {
#      if(input$xenaCohorts != "" && input$xenaCohortDatasets != ""){
#        print(file)
#        write.csv(downloadFile(), file)
#      }
#    }
#  )
  

})
