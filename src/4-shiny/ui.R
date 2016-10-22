shinyUI(
  
  fluidPage(
  
    sidebarLayout(
      
      sidebarPanel(
        titlePanel("multiOmics"),
        fileInput("mrnaFile", accept=c("text/csv"), label=h4("mRNA profile")),
        fileInput("mirnaFile",accept=c("text/csv"), label=h4("miRNA profile")),
        sliderInput("thresholdSlider", label=h4("Threshold"), min=0, max=1, value=0.7, step=0.05),
        # directoryInput('outputDir', label = 'Select a directory for the output', value = '~'),
        # textInput('outputFileName', label = 'Output file name', value = 'inputStep2-matureMirnaXmrna.csv'),
        actionButton("runCalc", "Run pipeline")
      ),
      mainPanel(
        DT::dataTableOutput('result'),
        plotOutput('correlationPlot')
      )
    )
  )
)