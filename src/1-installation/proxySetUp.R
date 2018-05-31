#seteo de proxy para RCurl (multimir)
opts <- list(proxy="proxy-url", proxyport=8080)
options(RCurlOptions = opts)
library("RCurl")
getURL("http://www.google.com.ar")


#seteo para httr (xenar )
library("httr")
set_config(use_proxy(url="proxy-url",port=8080))
GET("http://www.google.com.ar")
