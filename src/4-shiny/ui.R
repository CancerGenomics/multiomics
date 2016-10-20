shinyUI(
  
  fluidPage(
  
    sidebarLayout(
      
      sidebarPanel(
        titlePanel("mRna-miRna analysis"),
        fileInput("mrnaFile", accept=c("text/csv"), label=h4("mrna File")),
        fileInput("mirnaFile",accept=c("text/csv"), label=h4("mirna File")),
        sliderInput("thresholdSlider", label=h4("Threshold"), min=0, max=1, value=0.7, step=0.05),
        # directoryInput('outputDir', label = 'Select a directory for the output', value = '~'),
        # textInput('outputFileName', label = 'Output file name', value = 'inputStep2-matureMirnaXmrna.csv'),
        actionButton("runCalc", "Run Calculation")
      ),
      mainPanel(
        DT::dataTableOutput('result'),
        plotOutput('correlationPlot')
      )
    )
  )
)