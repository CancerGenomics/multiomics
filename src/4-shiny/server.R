# max size is 500MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
  
  sharedValues <- reactiveValues(fromButton=F,correlations="",correlationsStep2="",cnvMrnaCorrelations="")
  
  ###########################################################################
  ########################## MIRNA - MRNA PIPELINE TAB
  ###########################################################################
  
  observeEvent(input$runMRNAMiRNACorrelation, { 
    if(!is.null(input$mrnaFile) && !is.null(input$mirnaFile)) {
      
      sharedValues$fromButton <- T
      sharedValues$correlationsStep2 <- matrix()
        
      runMRNAMiRNACorrelation()
      
      #chequear si quiere correr el step 2 y mandar a correr
      if(input$miRNA.runMultimir) {
        print("running with step 2")
        runMultimirAnalisys()
        colnames(sharedValues$correlationsStep2) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value","miRNA db","Predicted score","PubMed ID")
        matrix.to.render <- sharedValues$correlationsStep2
      } else {
        print("running only step 1")
        colnames(sharedValues$correlations) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value")
        matrix.to.render <- sharedValues$correlations
      }

      # must add clinical data to the matrix to render
      if(!is.null(input$mirna.survivalFile)) {
        number.of.clusters=1
        clinical.survival.column.name="OVERALL_SURVIVAL"
        clinical.event.column.name="overall_survival_indicator"
        progResult <- getPrognosticStatistic(mrnaExpressionData(), number.of.clusters, groupin.FUN=multiomics.cut2, 
                                             input$mirna.survivalFile$datapath, clinical.survival.column.name, 
                                             clinical.event.column.name, minimium.number.of.samples.in.a.group=10)

        # creating a matrix to bind to the actual result        
        tmp <- matrix(nrow = nrow(sharedValues$correlations), ncol = ncol(progResult)-1 )
        # add new colnames to show
        colnames(tmp) <- colnames(progResult)[2:ncol(progResult)]
        matrix.to.render <- cbind(sharedValues$correlations,tmp)
        
        # loop into the matrix to plot looking for the corresponding values per gen
        for (i in 1:nrow(matrix.to.render)) {
          # actual gen from row i
          gen <- matrix.to.render[i,1]
          # values corresponding to actual gen to add columns
          row <- which(progResult[,1] == gen)
          # adding columns for actual gen row
          matrix.to.render[i,(ncol(sharedValues$correlations) +1):ncol(matrix.to.render)] <- progResult[row,2:ncol(progResult)]
        }
        
      }      
      
      # render the matrix with corresponding values
      output$result <- DT::renderDataTable(matrix.to.render, selection = 'single')
	  
  	  if(!input$miRNA.runMultimir) {
	      if(nrow(matrix.to.render) > 0) {
		      shinyjs::show(id = "downloadMrnaMirnaResult")
	      } else {
  		    shinyjs::hide(id = "downloadMrnaMirnaResult")
	      }
	    }

      sharedValues$fromButton <- F      
      
    } else {
      print("No hay archivos cargados")
    }
    
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
            
       #Checks if both files has the same samples in the same order. If not, aborts the execution.
       print("Checking if both files has the same samples in the same order...")
       suppressWarnings(checkSamplesFormIRNArnaCorrelation(mrnaExpressionData(), mirnaExpressionData(), 1))
       print("Preparing...")
       correlations()
        
  })
  } 
  
  runMultimirAnalisys <- function() { 
	  withProgress(message = 'Please stand by...', 
			           detail = "Searching in miRNA databases...", 
			           min=0, max=1, {
				  
			output.file.name = paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-multiMiR-outputFile.csv", sep = "")    
			print(output.file.name)
      step2Res <- keepBestGeneXMirnaAccordingCorrelationAndAddMirnaDbInfo(sharedValues$correlations,output.path="~/",
                                                              output.file.name, predicted.cut.off=10
                                                              )
      print("Finish multimir analisys")
      print(nrow(step2Res))
      collapsed.output.file.name = paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-multiMiR-collapsed-outputFile.csv", sep = "")    
      collapsedResult <- ColapseMirnaXMrna(data.frame(step2Res), output.path = "~/", output.file = collapsed.output.file.name)
      print("Finish collapsing")
      
      sharedValues$correlationsStep2 <- collapsedResult
      
	})
  }   
  
  output$downloadMrnaMirnaResult <- downloadHandler(
    filename = function() { 
      if(input$miRNA.runMultimir) {
        paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-multiMiR-collapsed-outputFile.csv", sep = "")
      } else {
        paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-outputFile.csv", sep = "")
      }
    },
    content = function(file) {
      if(input$miRNA.runMultimir) {
        write.csv(sharedValues$correlationsStep2, file)
      } else {
        write.csv(sharedValues$correlations, file)
      }
    })  

  output$correlationPlot <- renderPlot({
    if(!is.null(input$result_rows_selected)){
      selected.gene <- correlations()[input$result_rows_selected,1]
      selected.mirna <- correlations()[input$result_rows_selected,2]
      selected.gene.row <- which(mrnaExpressionData()==selected.gene)
      selected.mirna.row <- which(mirnaExpressionData()==selected.mirna)
      X <- as.numeric(as.vector(mirnaExpressionData()[selected.mirna.row,2:ncol(mirnaExpressionData())]))
      Y <- as.numeric(as.vector(mrnaExpressionData()[selected.gene.row,2:ncol(mrnaExpressionData())]))
      cor.test(X, Y)
      plot(X, Y, xlab=selected.mirna, ylab=selected.gene, main='miRNA vs. mRNA correlation plot', col='Black', pch=21, frame.plot=TRUE) #col=Group
      line <- lm(Y ~ X)
      abline(line, col="blue")
    }
  })

  
  output$correlationSurvival <- renderPlot({
    ERROR.GROUPING="ERROR.GROUPING"
    ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING="ERROR.EXECUTING.SURVFIT.FOR.PLOTTING"
    
    if(!is.null(input$result_rows_selected)){

      selected.gene <- correlations()[input$result_rows_selected,1]
      selected.gene.row <- which(mrnaExpressionData()==selected.gene)
      expression.vector <- as.numeric(as.vector(mrnaExpressionData()[selected.gene.row,2:ncol(mrnaExpressionData())]))

      ######ENGANCHAR CON LA GUI####
      #time<-read.time()
      #event<-read.event()
      #number.of.clusters<-read.number.of.clusters()
      #grouping.FUN<-read.grouping.fun()
      #minimium.number.of.samples.in.a.group<-read.minimium.number.of.samples.in.a.group()
      
      #####ELIMINAR ESTE HARDCODEO UNA VEZ ENGANCHADO########
      minimium.number.of.samples.in.a.group<-1
      event<-c(0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,0,0,1,1,1,0,0,0,0)
      time<-c(0.00000,22.17659,181.25670,36.43532,105.62630,64.98563,123.53180,46.88296,133.88090,94.48871,121.65910,211.90970,31.29960,66.04920,31.43040,210.89120,27.48600,17.83984,44.84400,46.19302,153.19920,155.95890,181.61810,24.77207,169.52770,112.88710,65.83984,220.09030,32.49281,32.51520,95.14560,49.41360,66.66119,153.16630,158.09450)
      number.of.clusters<-2
      grouping.FUN<-multiomics.cut2

      
      tryCatch({
        #Grouping
        tryCatch(
          {	groups<-grouping.FUN(expression.vector, number.of.clusters)
          #Validation for checking if all groups have got well formed groups
          result.check.groups.are.well.formed<-checkGroupsAreWellFormed(groups, selected.gene, minimium.number.of.samples.in.a.group)	
          if (!result.check.groups.are.well.formed@OK) stop(result.check.groups.are.well.formed@message)
          
          },error=function(e){stop(formatErrorMessage(error.type=ERROR.GROUPING, error.detail=e$message))})
        
          n.groups <- length(unique(groups))
        
          #SurvFit for plotting
          tryCatch({
           surv.fit<-survfit(formula = Surv(time, event) ~ groups)
            #Los colores se asignan as?: el primer color del vector col se asigna al grupo m?s chico, el segundo color al segundo m?s chico y as? siguiendo.
           #Es decir, no importa el orden en el que aparezcan los elementos en el time.Considera solamente el groups. 
           plot(surv.fit,col=c("blue", "red"), xlab="Time", ylab="survival")
           title<-selected.gene
           drawLegends(groups, "Survival", title)
          },error=function(e){stop(formatErrorMessage(error.type=ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING, error.detail=e$message))})
      })
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

					  if(nrow(cnvMrnaCorrelations()) > 1) {
					    output$MRNACNVResult <- DT::renderDataTable(cnvMrnaCorrelations(), selection = 'single')
					    shinyjs::show(id = "downloadMrnaCNVResult")
						  print("Hay resultados")
					  } else {
						  print("NO hay resultados")
					    shinyjs::hide(id = "downloadMrnaCNVResult")
					  }      
					  sharedValues$fromButton <- F
				  } else {
					  print("No hay archivos cargados")
				  }
				  
    })
  }
  
  output$downloadMrnaCNVResult <- downloadHandler(
    filename = function() { 
      paste(input$cnv.cnvFile$name,"-",input$cnv.cnvFile$name,"-outputFile.csv", sep = "")
    },
    content = function(file) {
      write.csv(cnvMrnaCorrelations(), file)
    }) 
  
  ###########################################################################
  ########################## METHYLATION - MRNA PIPELINE TAB
  ###########################################################################
  
  observeEvent(input$runMRNAMethylationCorrelation, { 
    runMRNAMethylationCorrelation()
  })  
  
  runMRNAMethylationCorrelation <- function(){
    withProgress(message = 'Please stand by...', 
                 detail = "calculating correlation", 
                 min=0, max=1, {
                   
      print("Implementar!")
                   
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

})
