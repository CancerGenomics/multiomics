# max size is 500MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
  
  sharedValues <- reactiveValues(fromButton=F,correlations="",correlationsStep2="",
                                 mirna.matrix.to.render="",cnvMrnaCorrelations="",
                                 cnv.matrix.to.render="", meth.matrix.to.render="",
                                 methMrnaCorrelations="", xena.datasets="")
  
  
  geneListToClipboard <- function(gene.list){
    gene.list <- unique(gene.list)
    str <- ""
    for(i in 1:length(gene.list)) {
      str <- paste0(str,gene.list[i],"\n")
    }
    return(str)
  }
  
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
        ###MDB: 26/2/2018 - P.ADJUST
        #colnames(sharedValues$correlationsStep2) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value","miRNA db","Predicted score","PubMed ID")
        #colnames(sharedValues$correlationsStep2) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value", "adj-p-value", "id", "miRNA db","Predicted score","PubMed ID")
        colnames(sharedValues$correlationsStep2) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value", "adj-p-value", "miRNA db","Predicted score","PubMed ID")
        sharedValues$mirna.matrix.to.render <- sharedValues$correlationsStep2
      } else {
        print("running only step 1")
        
        ###MDB: 26/2/2018 - P.ADJUST
        #colnames(sharedValues$correlations) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value")
        colnames(sharedValues$correlations) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value", "adj-p-value")
        
        sharedValues$mirna.matrix.to.render <- sharedValues$correlations
      }
      
      # must add clinical data to the matrix to render
      if((!is.null(input$mirna.survivalFile) && (nrow(sharedValues$correlationsStep2)>0))) {
        number.of.clusters=1
        progResult <- getPrognosticStatistic(mrnaExpressionData(), number.of.clusters, groupin.FUN=multiomics.cut2, 
                                             input$mirna.survivalFile$datapath, input$mirna.survival.column.name , 
                                             input$mirna.event.column.name, minimium.number.of.samples.in.a.group=10)
        
        
        if(input$miRNA.runMultimir) {
          ###MDB: 26/2/2018 - P.ADJUST
          #colnames(sharedValues$correlationsStep2) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value","miRNA db","Predicted score","PubMed ID")
          colnames(sharedValues$correlationsStep2) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value", "adj-p-value", "miRNA db","Predicted score","PubMed ID" )
          
          # creating a matrix to bind to the actual result        
          tmp <- matrix(nrow = nrow(sharedValues$correlationsStep2), ncol = ncol(progResult)-1 )
          actual.result <- sharedValues$correlationsStep2
        } else {
          ###MDB: 26/2/2018 - P.ADJUST
          #colnames(sharedValues$correlations) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value")
          colnames(sharedValues$correlations) <- c("Gene","Mature miRNA","miRNA-mRNA correlation","p-value","adj-p-value")
          
          # creating a matrix to bind to the actual result        
          tmp <- matrix(nrow = nrow(sharedValues$correlations), ncol = ncol(progResult)-1 )
          actual.result <- sharedValues$correlations
        }
        
        # add new colnames to show
        colnames(tmp) <- colnames(progResult)[2:ncol(progResult)]
        sharedValues$mirna.matrix.to.render <- cbind(actual.result,tmp)
        
        # loop into the matrix to plot looking for the corresponding values per gen
        for (i in 1:nrow(sharedValues$mirna.matrix.to.render)) {
          # actual gen from row i
          gen <- sharedValues$mirna.matrix.to.render[i,1]
          # values corresponding to actual gen to add columns
          row <- which(progResult[,1] == gen)
          # adding columns for actual gen row
          sharedValues$mirna.matrix.to.render[i,(ncol(actual.result) +1):ncol(sharedValues$mirna.matrix.to.render)] <- progResult[row,2:ncol(progResult)]
        }
      }      
      
      # render the matrix with corresponding values
      output$result <- DT::renderDataTable(sharedValues$mirna.matrix.to.render, selection = 'single')
      
      if(!input$miRNA.runMultimir) {
        if(nrow(sharedValues$mirna.matrix.to.render) > 0) {
          shinyjs::show(id = "downloadMrnaMirnaResult")
          shinyjs::show(id = "mirnaClip")
          shinyjs::hide(id = "mirnaEmptyResultMessage")
        } else {
          shinyjs::hide(id = "downloadMrnaMirnaResult")
          shinyjs::hide(id = "mirnaClip")
          shinyjs::show(id = "mirnaEmptyResultMessage")
        }
      }
      
      sharedValues$fromButton <- F      
      
    } else {
      openGeneralInformationMessage("You must select at least mRNA and miRNA file")
      print("No hay archivos cargados para el pipeline de mirna")
    }
    
  })
  
  openGeneralInformationMessage <- function(message) {
    output$generalMessageOutputText <- renderText(message)
    shinyBS::toggleModal(session = session, modalId = "generalMessageModal", toggle = "open")
    
  }
  
  mrnaExpressionReadData <- reactive({
			  print("mrnaExpressionData")
			  return(readMrnaExpressionFile(input$mrnaFile$datapath))
		  })
  
  mirnaExpressionReadData <- reactive({
			  print("mirnaExpressionData")
			  readMirnaExpressionFile(input$mirnaFile$datapath)
		  })
  
  mrnaExpressionData <- reactive({
    return(mirnaMrnaIntersection()[[1]])
  })
  
  mirnaExpressionData <- reactive({
    return(mirnaMrnaIntersection()[[2]])
  })

  mirnaMrnaIntersection <- reactive({
	#Keep columns which are in both databases
	intersection<-keepSameColumns(mrnaExpressionReadData(),mirnaExpressionReadData())
	return(intersection)
  })
  
  threshold <- reactive({
    input$thresholdSlider
  })
  
  pearsonsMethod <- reactive({
    input$pearsons.method
  })
  
  correlations <- reactive(quote({
    if(sharedValues$fromButton) {

      keep.pos.cor <- (input$mirna.correlation.type == "positive") || (input$mirna.correlation.type == "both")
      keep.neg.cor <- (input$mirna.correlation.type == "negative") || (input$mirna.correlation.type == "both")
      
      sharedValues$correlations <- CalculateCorrelationsMirnaMrnaUsingWCGNA(
        mrnaExpressionData(), mirnaExpressionData(),
        output.path="~/", 
        output.file.name = paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-outputFile.csv", sep = ""),
        r.minimium = threshold(), inc.progress = T, keep.pos.cor = keep.pos.cor, keep.neg.cor = keep.neg.cor)
    } 
    return (sharedValues$correlations)
  }), quoted = T)
  
  
  runMRNAMiRNACorrelation <- function() { 
    withProgress(message = 'Please stand by...', 
                 detail = "calculating correlation", 
                 min=0, max=1, {
                   
                   #Checks if both files has the same samples in the same order. If not, aborts the execution.
                   #print("Checking if both files has the same samples in the same order...")
                   #suppressWarnings(checkSamplesFormIRNArnaCorrelation(mrnaExpressionData(), mirnaExpressionData(), 1))
                   print("Preparing...")
                   correlations()
                   
                   # Add clipboard buttons
                   output$mirnaClip <- renderUI({
                     rclipButton("mirnaCopyToClipboard", "Copy genes to clipboard", geneListToClipboard(correlations()[,1]) , icon("clipboard"))
                   })      
                   
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
                   if (nrow(step2Res)>0){
                     collapsed.output.file.name = paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-multiMiR-collapsed-outputFile.csv", sep = "")    
                     collapsedResult <- ColapseMirnaXMrna(data.frame(step2Res), output.path = "~/", output.file = collapsed.output.file.name)
                     print("Finish collapsing")
                   }
                   else
                   {
                     collapsedResult<-data.frame(matrix(nrow = 0, ncol = 8))
                     
                   }
                   sharedValues$correlationsStep2 <- collapsedResult
                   
                 })
  }  
  
  observeEvent(input$mirna.survivalFile,{
    clinical.data <- na.omit(read.table(input$mirna.survivalFile$datapath, nrows = 1, header=TRUE,fill=TRUE))
    updateClinicalDataSelects(session, clinical.data, "mirna.survival.column.name", "mirna.event.column.name")
  })
  
  updateClinicalDataSelects <- function(session, clinical.data, survivalSelectName, eventSelectName) {
    if("OVERALL_SURVIVAL" %in% colnames(clinical.data)) {
      updateSelectInput(session,survivalSelectName, choices = colnames(clinical.data), selected = "OVERALL_SURVIVAL")
    } else {
      updateSelectInput(session,survivalSelectName, choices = colnames(clinical.data))
    }
    shinyjs::show(survivalSelectName)
    
    if("overall_survival_indicator" %in% colnames(clinical.data)) {
      updateSelectInput(session,eventSelectName, choices = colnames(clinical.data), selected = "overall_survival_indicator")
    } else {
      updateSelectInput(session,eventSelectName, choices = colnames(clinical.data))
    }
    shinyjs::show(eventSelectName)   
  }
  
  validatemirnaExpressionData <- function(datapath) {
    result <- tryCatch({
      readMirnaExpressionFile(datapath)
      result <- T
    }, error = function(e) {
      result <- F
    })
    return(result)
  }  
  
  observeEvent(input$mirnaFile,{
    valid <- validatemirnaExpressionData(input$mirnaFile$datapath)
    if(!valid){
      shinyjs::show(id="mirnaMirnaFileErrorMsg")
      sharedValues$mirna.matrix.to.render <- NULL
      output$result <- DT::renderDataTable(sharedValues$mirna.matrix.to.render, selection = 'single')      
      shinyjs::hide(id="downloadMrnaMirnaResult")
      shinyjs::hide(id="mirnaClip")
      
    } else {
      shinyjs::hide(id="mirnaMirnaFileErrorMsg")
    }
    shiny::validate(
      need(valid,"miRNA file has a wrong format, please verify the file and retry")
    )
    
  })  
  
  validatemrnaExpressionData <- function(datapath) {
    result <- tryCatch({
      readMrnaExpressionFile(datapath)
      result <- T
    }, error = function(e) {
      result <- F
    })
    return(result)
  }
  
  observeEvent(input$mrnaFile,{
    valid <- validatemrnaExpressionData(input$mrnaFile$datapath)
    if(!valid){
      shinyjs::show(id="mirnaMrnaFileErrorMsg")
      sharedValues$mirna.matrix.to.render <- NULL
      output$result <- DT::renderDataTable(sharedValues$mirna.matrix.to.render, selection = 'single')      
      shinyjs::hide(id="downloadMrnaMirnaResult")
      shinyjs::hide(id="mirnaClip")
    } else {
      shinyjs::hide(id="mirnaMrnaFileErrorMsg")
    }
    shiny::validate(
      need(valid,"mRNA file has a wrong format, please verify the file and retry")
    )     
    
  })  
  
  output$downloadMrnaMirnaResult <- downloadHandler(
    filename = function() { 
      if(input$miRNA.runMultimir) {
        paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-multiMiR-collapsed-outputFile.csv", sep = "")
      } else {
        paste(input$mirnaFile$name,"-",input$mrnaFile$name,"-outputFile.csv", sep = "")
      }
    },
    content = function(file) {
      write.csv(sharedValues$mirna.matrix.to.render, file)
    })  
  
  output$correlationPlot <- renderPlot({
    if(!is.null(input$result_rows_selected)){
      selected.gene <- correlations()[input$result_rows_selected,1]
      selected.mirna <- correlations()[input$result_rows_selected,2]
      selected.gene.row <- which(mrnaExpressionData()==selected.gene)
      selected.mirna.row <- which(mirnaExpressionData()==selected.mirna)
      
      mrna.sample.names<- colnames(mrnaExpressionData())[2:ncol(mrnaExpressionData())]
      mirna.sample.names<- colnames(mirnaExpressionData())[2:ncol(mirnaExpressionData())]
      cols.in.common<-intersect(mrna.sample.names, mirna.sample.names)
      
      #X <- as.numeric(as.vector(mirnaExpressionData()[selected.mirna.row,2:ncol(mirnaExpressionData())]))
      #Y <- as.numeric(as.vector(mrnaExpressionData()[selected.gene.row,2:ncol(mrnaExpressionData())]))
      
      X <- as.numeric(as.vector(mirnaExpressionData()[selected.mirna.row,cols.in.common]))
      Y <- as.numeric(as.vector(mrnaExpressionData()[selected.gene.row,cols.in.common]))
      
      cor.test(X, Y)
      plot(X, Y, xlab=selected.mirna, ylab=selected.gene, main='miRNA vs. mRNA correlation plot', col='Black', pch=21, frame.plot=TRUE) #col=Group
      line <- lm(Y ~ X)
      abline(line, col="blue")
    }
  })
  
  
  output$correlationSurvival <- renderPlot({
    ERROR.GROUPING="ERROR.GROUPING"
    ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING="ERROR.EXECUTING.SURVFIT.FOR.PLOTTING"
    
    if(!is.null(input$mirna.survivalFile) && !is.null(input$result_rows_selected)){
      
      selected.gene <- correlations()[input$result_rows_selected,1]
      selected.gene.row <- which(mrnaExpressionData()==selected.gene)
      expression.vector <- as.numeric(as.vector(mrnaExpressionData()[selected.gene.row,2:ncol(mrnaExpressionData())]))
      
      
      mirna.survival.matrix <- read.table(input$mirna.survivalFile$datapath, header = T)
      mirna.survival.matrix<-mirna.survival.matrix[mirna.survival.matrix$X_SAMPLE_ID %in% colnames(mirnaExpressionData()),]
      survival.name.col <- which(colnames(mirna.survival.matrix)==input$mirna.survival.column.name)
      survival.event.col <- which(colnames(mirna.survival.matrix)==input$mirna.event.column.name)
      
      ######ENGANCHAR CON LA GUI Y ELIMINAR HARCODEO ####
      #number.of.clusters<-read.number.of.clusters()
      number.of.clusters<-2
      #grouping.FUN<-read.grouping.fun()
      grouping.FUN<-multiomics.cut2
      #minimium.number.of.samples.in.a.group<-read.minimium.number.of.samples.in.a.group()
      minimium.number.of.samples.in.a.group<-1
      
      time<-as.vector(mirna.survival.matrix[,survival.name.col])
      event<-as.vector(mirna.survival.matrix[,survival.event.col])
     
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
  
  cnvMrnaExpressionReadData <- reactive({
    print("cnvMrnaExpressionReadData")
    readMrnaExpressionFile(input$cnv.mrnaFile$datapath)
  })
  
  cnvExpressionReadData <- reactive({
    print("cnvExpressionReadData")
    readCNVFile(input$cnv.cnvFile$datapath)
  })
  
  cnvPearsonsMethod <- reactive({
    input$cnv.pearsons.method
  })

  cnvMrnaExpressionData <- reactive({
	return(cnvMrnaIntersection()[[1]])
  })

  cnvExpressionData <- reactive({
	return(cnvMrnaIntersection()[[2]])
  })

  cnvMrnaIntersection <- reactive({
	#Keep columns which are in both databases
	intersection<-keepSameColumns(cnvMrnaExpressionReadData(),cnvExpressionReadData())
	return(intersection)
  })
  
  cnvMrnaCorrelations <- reactive(quote({
    if(sharedValues$fromButton) {
      
      keep.pos.cor <- (input$cnv.correlation.type == "positive") || (input$cnv.correlation.type == "both")
      keep.neg.cor <- (input$cnv.correlation.type == "negative") || (input$cnv.correlation.type == "both")
      
      sharedValues$cnvMrnaCorrelations <- CnvXMrnas(cnvMrnaExpressionData(), cnvExpressionData(), output.path="~/", 
                                                         output.file.name=paste(input$cnv.mrnaFile$name,"-",input$cnv.cnvFile$name,"-outputFile.csv", sep = ""),
                                                         r.minimium = cnvThreshold(), inc.progress = T,keep.pos.cor=keep.pos.cor,keep.neg.cor=keep.neg.cor)
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
                     cnvMrnaCorrelations()
                     sharedValues$cnv.matrix.to.render <- sharedValues$cnvMrnaCorrelations
                     
                     if(!is.null(input$cnv.survivalFile) && (nrow(sharedValues$cnv.matrix.to.render) > 0)) {
                       number.of.clusters=1
                       print(cnvMrnaExpressionData())
                       progResult <- getPrognosticStatistic(cnvMrnaExpressionData(), number.of.clusters, groupin.FUN=multiomics.cut2, 
                                                            input$cnv.survivalFile$datapath, input$cnv.survival.column.name , 
                                                            input$cnv.event.column.name, minimium.number.of.samples.in.a.group=10)
                       colnames(sharedValues$cnvMrnaCorrelations) <- c("Gene","Location","CNV-mRNA correlation","p-value", "p_value_fdr_adjusted", "ID")
                       # creating a matrix to bind to the actual result        
                       tmp <- matrix(nrow = nrow(sharedValues$cnvMrnaCorrelations), ncol = ncol(progResult)-1 )
                       actual.result <- sharedValues$cnvMrnaCorrelations
                       
                       # add new colnames to show
                       colnames(tmp) <- colnames(progResult)[2:ncol(progResult)]
                       sharedValues$cnv.matrix.to.render <- cbind(actual.result,tmp)
                       
                       # loop into the matrix to plot looking for the corresponding values per gen
                       for (i in 1:nrow(sharedValues$cnv.matrix.to.render)) {
                         # actual gen from row i
                         gen <- sharedValues$cnv.matrix.to.render[i,1]
                         # values corresponding to actual gen to add columns
                         row <- which(progResult[,1] == gen)
                         if(!isEmpty(row)){
                           # adding columns for actual gen row
                           sharedValues$cnv.matrix.to.render[i,(ncol(actual.result) +1):ncol(sharedValues$cnv.matrix.to.render)] <- progResult[row,2:ncol(progResult)]
                         }
                       }					    
                     }
                     
                     if(nrow(sharedValues$cnv.matrix.to.render) > 0) {
                       output$MRNACNVResult <- DT::renderDataTable(sharedValues$cnv.matrix.to.render, selection = 'single')
                       shinyjs::show(id = "downloadMrnaCNVResult")
                       # Add clipboard buttons
                       output$cnvClip <- renderUI({
                         rclipButton("cnvCopyToClipboard", "Copy genes to clipboard", geneListToClipboard(cnvMrnaCorrelations()[,1]) , icon("clipboard"))
                       })  					  
                       
                       shinyjs::show(id = "cnvClip")
                       shinyjs::hide(id = "cnvEmptyResultMessage")
                       print("Hay resultados")
                     } else {
                       print("NO hay resultados")
                       shinyjs::hide(id = "downloadMrnaCNVResult")
                       shinyjs::hide(id = "cnvClip")
                       shinyjs::show(id = "cnvEmptyResultMessage")
                     }      
                     sharedValues$fromButton <- F
                   } else {
                     openGeneralInformationMessage("You must select at least CNV and miRNA file")
                     print("No hay archivos cargados para el pipeline de CNV")
                   }
                   
                 })
  }
  
  output$downloadMrnaCNVResult <- downloadHandler(
    filename = function() { 
      paste(input$cnv.cnvFile$name,"-",input$cnv.mrnaFile$name,"-outputFile.csv", sep = "")
    },
    content = function(file) {
      write.csv(sharedValues$cnv.matrix.to.render, file)
    }) 
  
  observeEvent(input$cnv.survivalFile,{
    clinical.data <- na.omit(read.table(input$cnv.survivalFile$datapath, nrows = 1, header=TRUE,fill=TRUE))
    updateClinicalDataSelects(session, clinical.data, "cnv.survival.column.name", "cnv.event.column.name")
    
  })  
  
  validateCnvExpressionData <- function(datapath) {
    result <- tryCatch({
      readCNVFile(datapath)
      result <- T
    }, error = function(e) {
      result <- F
    })
    return(result)
  }  
  
  observeEvent(input$cnv.cnvFile,{
    valid <- validateCnvExpressionData(input$cnv.cnvFile$datapath)
    if(!valid){
      shinyjs::show(id="cnvFileErrorMsg")
      sharedValues$cnv.matrix.to.render <- NULL
      output$MRNACNVResult <- DT::renderDataTable(sharedValues$cnv.matrix.to.render, selection = 'single')      
      shinyjs::hide(id="downloadMrnaCNVResult")
      shinyjs::hide(id="cnvClip")
    } else {
      shinyjs::hide(id="cnvFileErrorMsg")
    }
    shiny::validate(
      need(valid,"CNV file has a wrong format, please verify the file and retry")
    )
    
  })    
  
  observeEvent(input$cnv.mrnaFile,{
    valid <- validatemrnaExpressionData(input$cnv.mrnaFile$datapath)
    if(!valid){
      shinyjs::show(id="cnvMrnaFileErrorMsg")
      sharedValues$cnv.matrix.to.render <- NULL
      output$MRNACNVResult <- DT::renderDataTable(sharedValues$cnv.matrix.to.render, selection = 'single')      
      shinyjs::hide(id="downloadMrnaCNVResult")
      shinyjs::hide(id="cnvClip")
    } else {
      shinyjs::hide(id="cnvMrnaFileErrorMsg")
    }
    shiny::validate(
      need(valid,"mRNA file has a wrong format, please verify the file and retry")
    )     
    
  })   
  
  output$cnv.correlationPlot <- renderPlot({
    if(!is.null(input$MRNACNVResult_rows_selected)){
      selected.gene <- cnvMrnaCorrelations()[input$MRNACNVResult_rows_selected,1]
      selected.gene.row <- which(cnvMrnaExpressionData()==selected.gene)
      selected.cnv.row <- which(cnvExpressionData()==selected.gene)
      X <- as.numeric(as.vector(cnvExpressionData()[selected.cnv.row,2:ncol(cnvExpressionData())]))
      Y <- as.numeric(as.vector(cnvMrnaExpressionData()[selected.gene.row,2:ncol(cnvMrnaExpressionData())]))
      cor.test(X, Y)
      xlab <- paste0(selected.gene," Mrna")
      ylab <- paste0(selected.gene," CNV")
      plot(X, Y, xlab=xlab, ylab=ylab, main='CNV vs. mRNA correlation plot', col='Black', pch=21, frame.plot=TRUE) #col=Group
      line <- lm(Y ~ X)
      abline(line, col="blue")
    }
  })
  
  
  output$cnv.correlationSurvival <- renderPlot({
    ERROR.GROUPING="ERROR.GROUPING"
    ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING="ERROR.EXECUTING.SURVFIT.FOR.PLOTTING"
    
    if(!is.null(input$cnv.survivalFile) && !is.null(input$MRNACNVResult_rows_selected)){
      
      selected.gene <- cnvMrnaCorrelations()[input$MRNACNVResult_rows_selected,1]
      selected.gene.row <- which(cnvExpressionData()==selected.gene)
      expression.vector <- as.numeric(as.vector(cnvExpressionData()[selected.gene.row,2:ncol(cnvExpressionData())]))
      
      
      cnv.survival.matrix <- read.table(input$cnv.survivalFile$datapath, header = T, stringsAsFactors = F)
      cnv.survival.matrix<-cnv.survival.matrix[cnv.survival.matrix$X_SAMPLE_ID %in% colnames(cnvExpressionData()),]
      survival.name.col <- which(colnames(cnv.survival.matrix)==input$cnv.survival.column.name)
      survival.event.col <- which(colnames(cnv.survival.matrix)==input$cnv.event.column.name)
      
      ######ENGANCHAR CON LA GUI Y ELIMINAR HARCODEO ####
      number.of.clusters<-2
      grouping.FUN<-multiomics.cut2
      minimium.number.of.samples.in.a.group<-1
      
      time<-as.vector(cnv.survival.matrix[,survival.name.col])
      event<-as.vector(cnv.survival.matrix[,survival.event.col])
      
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
        },error=function(e){
            print(e)
            stop(formatErrorMessage(error.type=ERROR.EXECUTING.SURV.FIT.FOR.PLOTTING, error.detail=e$message))
          })
      })
    }
  })    
  
  
  ###########################################################################
  ########################## METHYLATION - MRNA PIPELINE TAB
  ###########################################################################
  
  
  methMrnaCorrelations <- reactive(quote({
    if(sharedValues$fromButton) {
     
      keep.pos.cor <- (input$meth.correlation.type == "positive") || (input$meth.correlation.type == "both")
      keep.neg.cor <- (input$meth.correlation.type == "negative") || (input$meth.correlation.type == "both")
      
      sharedValues$methMrnaCorrelations <- methXMrnasWCGNA(methMrnaExpressionData(), methExpressionData(), methPlatform() , output.path="~/",
                                                           output.file.name=paste(input$meth.mrnaFile$name,"-",input$meth.methFile$name,"-outputFile.csv", sep = ""),           
                                                           r.minimium = methThreshold(), 
                                                           inc.progress = T,keep.pos.cor=keep.pos.cor,
                                                           keep.neg.cor=keep.neg.cor)
      
    }
    return (sharedValues$methMrnaCorrelations)
  }), quoted = T)  
  
  methPlatform <- reactive({
    return(getMethylationPlatformTable(input$meth.platform.select))
  })  
  
  methThreshold <- reactive({
    input$meth.thresholdSlider
  })
  
  methMrnaExpressionReadData <- reactive({
    readMrnaExpressionFile(input$meth.mrnaFile$datapath)
  })
  
  methExpressionReadData <- reactive({
    readMethylationFile(input$meth.methFile$datapath)
  })

  methMrnaExpressionData <- reactive({
			return(methMrnaIntersection()[[1]])
		})

  methExpressionData <- reactive({
			return(methMrnaIntersection()[[2]])
		})

  methMrnaIntersection <- reactive({
			#Keep columns which are in both databases
			intersection<-keepSameColumns(methMrnaExpressionReadData(),methExpressionReadData())
			return(intersection)
		})
  
  methPearsonsMethod <- reactive({
    input$meth.pearsons.method
  })
  
  
  observeEvent(input$runMRNAMethylationCorrelation, { 
    runMRNAMethylationCorrelation()
  })  
  
  
  runMRNAMethylationCorrelation <- function(){
    withProgress(message = 'Please stand by...', 
                 detail = "calculating correlation", 
                 min=0, max=1, {
                   
                   if(!is.null(input$meth.mrnaFile) && !is.null(input$meth.methFile)) {
                     
                     print("Preparing...")
                     
                     sharedValues$fromButton <- T
                     methMrnaCorrelations()
                     sharedValues$meth.matrix.to.render <- sharedValues$methMrnaCorrelations
                     
                     if(nrow(sharedValues$meth.matrix.to.render) > 0) {
                       output$MRNAMethResult <- DT::renderDataTable(sharedValues$meth.matrix.to.render, selection = 'single')
                       shinyjs::show(id = "downloadMrnaMethResult")
                       shinyjs::hide(id = "methEmptyResultMessage")
                       # Add clipboard buttons
                       output$methClip <- renderUI({
                         rclipButton("methCopyToClipboard", "Copy genes to clipboard", geneListToClipboard(methMrnaCorrelations()[,1]) , icon("clipboard"))
                       })   
                       shinyjs::show(id = "methClip")
                       print("Hay resultados")
                     } else {
                       print("NO hay resultados")
                       shinyjs::hide(id = "downloadMrnaMethResult")
                       shinyjs::hide(id = "methClip")
                       shinyjs::show(id = "methEmptyResultMessage")
                     }
                     sharedValues$fromButton <- F
                   } else {
                     openGeneralInformationMessage("You must select at least Methylation and miRNA file")
                     print("No hay archivos cargados para el pipeline de meth")
                   }
                   
                 })    
  }
  
  
  validateMethExpressionData <- function(datapath) {
    result <- tryCatch({
      readMethylationFile(datapath)
      result <- T
    }, error = function(e) {
      result <- F
    })
    return(result)
  }  
  
  observeEvent(input$meth.methFile,{
    valid <- validateCnvExpressionData(input$meth.methFile$datapath)
    if(!valid){
      shinyjs::show(id="methFileErrorMsg")
      sharedValues$meth.matrix.to.render <- NULL
      output$MRNAMethResult <- DT::renderDataTable(sharedValues$meth.matrix.to.render, selection = 'single')      
      shinyjs::hide(id="downloadMrnaMethResult")
      shinyjs::hide(id="methClip")
    } else {
      shinyjs::hide(id="methFileErrorMsg")
    }
    shiny::validate(
      need(valid,"Methylation file has a wrong format, please verify the file and retry")
    )
    
  })    
  
  observeEvent(input$meth.mrnaFile,{
    valid <- validatemrnaExpressionData(input$meth.mrnaFile$datapath)
    if(!valid){
      shinyjs::show(id="methMrnaFileErrorMsg")
      sharedValues$meth.matrix.to.render <- NULL
      output$MRNAMethResult <- DT::renderDataTable(sharedValues$meth.matrix.to.render, selection = 'single')      
      shinyjs::hide(id="downloadMrnaMethResult")
      shinyjs::hide(id="methClip")
    } else {
      shinyjs::hide(id="methMrnaFileErrorMsg")
    }
    shiny::validate(
      need(valid,"mRNA file has a wrong format, please verify the file and retry")
    )     
    
  })   
  
  
  
  output$downloadMrnaMethResult <- downloadHandler(
    filename = function() { 
      paste(input$meth.methFile$name,"-",input$meth.mrnaFile$name,"-outputFile.csv", sep = "")
    },
    content = function(file) {
      write.csv(sharedValues$meth.matrix.to.render, file)
    })   
  
  output$meth.correlationPlot <- renderPlot({
    if(!is.null(input$MRNAMethResult_rows_selected)){
      selected.gene <- methMrnaCorrelations()[input$MRNAMethResult_rows_selected,1]
      selected.meth <- methMrnaCorrelations()[input$MRNAMethResult_rows_selected,3]
      selected.gene.row <- which(methMrnaExpressionData()==selected.gene)
      selected.meth.row <- which(methExpressionData()==selected.meth)
      X <- as.numeric(as.vector(methExpressionData()[selected.meth.row,2:ncol(methExpressionData())]))
      Y <- as.numeric(as.vector(methMrnaExpressionData()[selected.gene.row,2:ncol(methMrnaExpressionData())]))
      cor.test(X, Y)
      plot(X, Y, xlab=selected.gene, ylab=selected.gene, main='Methylation vs. mRNA correlation plot', col='Black', pch=21, frame.plot=TRUE) #col=Group
      line <- lm(Y ~ X)
      abline(line, col="blue")
    }
  })
  
  
  ###########################################################################
  ########################## XENA HUB CONNECTOR TAB
  ###########################################################################
  
  mrnaCohortsData <- reactive({
    tryCatch(
      withProgress(message = 'Connectig with Xena...', detail = "Getting cohorts", min=0, max=1, {    
        print("Connecting to XenaHub...")
        cohorts(XenaHub(hosts = "https://tcga.xenahubs.net"))
      })  
      , error=function(e){
        tryCatch(
          {
            library(curl) 
            readLines(curl("https://tcga.xenahubs.net"))
            withProgress(message = 'Connectig with Xena...', detail = "Getting cohorts", min=0, max=1, {    
              print("Connecting to XenaHub...")
              cohorts(XenaHub(hosts = "https://tcga.xenahubs.net"))
            })
          }, error=function(e){print("Error connecting to https://tcga.xenahubs.net.")})
      })
  })
  
  
  
  observeEvent(input$connectToXenaHub, { 
    updateSelectInput(session, "xenaCohorts","XenaHub available cohorts", choices = mrnaCohortsData())
  })  
  
  observeEvent(input$xenaCohorts, {
    withProgress(message = 'Connectig with Xena...', 
                 detail = "Getting datasets", 
                 min=0, max=1, {    
                   if(input$xenaCohorts != "" ){
                     #print(input$xenaCohorts)
                     if(input$xenaCohorts != "(unassigned)"){
                       print(paste("Searching datesets for", input$xenaCohorts, sep=" "))
                       sharedValues$xena.datasets <- XenaR::datasets(XenaHub(hosts = "https://tcga.xenahubs.net", cohorts = input$xenaCohorts))
                       shinyjs::show("xenaCohortDatasets")
                       shinyjs::show("xenaCohortDatasetsFilter")
                     } else {
                       sharedValues$xena.datasets <- c("")
                       shinyjs::hide("xenaCohortDatasets")
                       shinyjs::hide("xenaCohortDatasetsFilter")
                     }
                     updateSelectInput(session, "xenaCohortDatasets","Cohort datasets", choices = sharedValues$xena.datasets)
                     updateSelectInput(session, "xenaCohortDatasetsFilter",selected = "All")
                   }
                 })   
  })   
  
  observeEvent(input$xenaCohortDatasetsFilter, {  
    
    if(!is.null(input$xenaCohorts)) {
      filtered <- sharedValues$xena.datasets
      print(input$xenaCohortDatasetsFilter)
      if((input$xenaCohortDatasetsFilter != "All")) {
        filtered <- filtered[grep(input$xenaCohortDatasetsFilter,filtered)]
        if(input$xenaCohortDatasetsFilter == "RNA") {
          mirna <- filtered[grep("miRNA",filtered)]
          filtered <- filtered[!filtered %in% mirna]
        }
      }
      
      updateSelectInput(session, "xenaCohortDatasets","Cohort datasets", choices = filtered)
    }
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
               "Download selected dataset",href=getUrlFromTCGAXenaHub(input$xenaCohortDatasets), target="_blank", 
               class="btn btn-default shiny-download-link"
        )
      )
    } else {
      tagList()
    }
  })  
  
  
})
