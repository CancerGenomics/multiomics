# multiOmics: an R package to infer genomics and epigenomics mechanisms involved with cancer progression  

## How to execute multiOmics?

### Requirements

RStudio and basic RStudio knowledge is required to run multiOmics.     

### First, install required packages 

Before running multiOmics, run [dependencies_multiomics.R](/src/1-installation/dependencies_multiomics.R) script or run the folowwing instructions in order to install the required dependencies:

``` R
install.packages("XML")    
install.packages("RCurl")    
multimirFile <- "https://gitlab.com/cancergenomics/multiomics/raw/master/src/1-installation/dependencies/multiMiR_1.0.1.tar.gz"    
install.packages(multimirFile, repos=NULL, type="source")    
install.packages("GGally")
install.packages("ggplot2")    
install.packages("Hmisc")  
install.packages("shiny")  
install.packages("DT")    
install.packages("shinyBS")      
install.packages("shinyjs")  
install.packages("rclipboard")
install.packages("clipr")
install.packages("heatmap.plus")  
source("http://bioconductor.org/biocLite.R")  
biocLite("ctc")  
biocLite("impute")  
biocLite("survcomp")  
biocLite("org.Hs.eg.db")  
biocLite("devtools")  
biocLite("mtmorgan/xenar")  
```

### Second, just run multiOmics!

You don't need to clone this git repository, RStudio will do it for you automatically. 
multiOmics provides a visual shiny application.   
In order to run this application from RStudio you can execute (after installing required packages):  

``` R
library("shiny")  
runGitHub("multiomics", "cancergenomics", subdir = "src/4-shiny/", ref = "multiomics-0.0.2-beta.12")
```

### Running cases studies for each pipeline

mRNA, CNV and DNA methylation data matrix are provided under the [example](/examples/) folder to test each multiOmics pipeline.


## Running behind internet proxy

multiOmics indirectly uses RCurl ahd httr libraries for accessing internet. 
If you are behind an internet proxy, you must run the following scripts for configuring proxy connections correctly:  

``` R
set proxy configuration here  
myProxyHost <- "set proxy host here"   
myProxyPort <- set proxy port here  

### RCurl configuration   
opts <- list(proxy=myProxyHost, proxyport=myProxyPort)  
options(RCurlOptions = opts)  
library("RCurl")  
### test connection  
getURL("http://www.google.com")  
  
### httr configuration  
library("httr")  
set_config(use_proxy(url=myProxyHost,port=myProxyPort))  
### test connection  
GET("http://www.google.com")  
```