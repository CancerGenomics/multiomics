# max size is 500MB
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output, session) {
  
  sharedValues <- reactiveValues(fromButton=F,correlations="")
  
  observeEvent(input$runCalc, { 
    runCalc()
    })
  observeEvent(input$input$result_rows_selected, {plotCorrelation()})
  
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
      sharedValues$correlations <- CalculateCorrelationsMirnaMrna(mrnaExpressionData(), mirnaExpressionData(),
                                   output.path="~/", 
                                   output.file.name = paste(input$mrnaFile$name,"-outputFile.csv", sep = ""),
                                   r.minimium = threshold(), inc.progress = T, pearsons.method = pearsonsMethod())
    } 
    return (sharedValues$correlations)
  }), quoted = T)
  
  runCalc <- function() { 
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
        #print(nrow(mrnaExpressionData()[selected.gene.row,]))
        #print(ncol(mrnaExpressionData()[selected.gene.row,]))
        #print(mirnaExpressionData()[selected.mirna.row,])
        X <- as.numeric(as.vector(mrnaExpressionData()[selected.gene.row,2:ncol(mrnaExpressionData())]))
        Y <- as.numeric(as.vector(mirnaExpressionData()[selected.mirna.row,2:ncol(mirnaExpressionData())]))
        #X<-c(50,8,70,65)
        #Y<-c(5,4,8,9) 
        cor.test(X, Y)
        plot(X, Y, col='Black', pch=1) #col=Group
        line <- lm(Y ~ X)
        abline(line, col="blue")
      }
    })    


})
