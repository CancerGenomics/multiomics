
# retorna una matriz con toda la info de las plataformas disponibles
getMethylationPlatforms <- function(path="") {
  if (path=="") path=paste0(getwd(),"/../../resources/methilation.platforms/illuminaMethyl450_hg19_GPL16304.txt")
  platforms <- matrix(nrow = 1,ncol = 2)
  platforms[1,] <- c("HumanMethylation450 BeadChip",path)
  return(platforms)
}

# retorna un arreglo con los nombres de las plataformas
getMethylationPlatformNames <- function() {
  return(getMethylationPlatforms()[,2])
}

# returna el contenido de la plataforma indicada
getMethylationPlatformTable <- function(meth.platform) {
  platforms <- getMethylationPlatforms()
  index <- which(platforms == meth.platform)
  return(read.table(platforms[index,2], header = TRUE, sep="\t"))
}

# returna el contenido de la plataforma indicada
getMethylationPlatformTableForPipeline <- function(meth.platform, path) {
  platforms <- getMethylationPlatforms(path)
  index <- which(platforms == meth.platform)
  return(read.table(platforms[index,2], header = TRUE, sep="\t"))
}