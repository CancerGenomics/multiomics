# multiOmics: an R package to infer genomics and epigenomics mechanisms involved with cancer progression  

## required packages installation 

Before running multiOmics, run this R scripts in order to install the required dependencies:

install.packages("XML")    
install.packages("RCurl")    
multimirFile <- "https://gitlab.com/cancergenomics/multiomics/raw/master/src/1-installation/dependencies/multiMiR_1.0.1.tar.gz"    
install.packages(multimirFile, repos=NULL, type="source")    
install.packages("GGally")    
install.packages("Hmisc")  
install.packages("shiny")  
install.packages("DT")    
install.packages("shinyBS")      
install.packages("shinyjs")  
install.packages("heatmap.plus")  
source("http://bioconductor.org/biocLite.R")  
biocLite("ctc")  
biocLite("impute")  
biocLite("survcomp")  
biocLite("org.Hs.eg.db")  
biocLite("devtools")  
biocLite("mtmorgan/xenar")  

## running behind proxy

multiOmics indirectly uses RCurl ahd httr libraries for accessing internet. 
You must run the following scripts for configuring proxy connections:  

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


## running multiOmics from github

multiOmics provides a visual shiny application.   
In order to run this application, from RStudio you can execute (after installing required packages):  

library("shiny")  
runGitHub("multiomics", "cancergenomics", subdir = "src/4-shiny/", ref = "multiomics-0.0.2-beta.10")

## Running cases studies for each pipeline

mRNA, CNV and DNA methylation data matrix are provided under the [example](/examples/) folder to test each multiOmics pipeline.
