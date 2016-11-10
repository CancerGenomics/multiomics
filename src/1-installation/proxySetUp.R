#seteo de proxy para RCurl (multimir)
opts <- list(proxy="vip-proxy-cli-ka", proxyport=8080)
options(RCurlOptions = opts)
library("RCurl")
getURL("http://www.google.com.ar")


#seteo para httr (xenar )
library("httr")
set_config(use_proxy(url="vip-proxy-cli-ka",port=8080))
GET("http://www.google.com.ar")
