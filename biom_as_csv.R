biom.as.csv <- function(path){
  if (!"dplyr"%in%row.names(installed.packages())) {
    install.packages("dplyr")
  }
  if (!"biomformat"%in%row.names(installed.packages())){
    source("https://bioconductor.org/biocLite.R")
    biocLite("biomformat")
  }
  require(dplyr)
  require(biomformat)
  
  biom <- read_biom(path)
  
  comm <- biom$data
  
  den <- function(x){
    y <- x$id
    return(y)
  }
  
  taxon <- function(x){
    y <- paste0(unlist(x$metadata), collapse="; ")
    return(y)
  }
  
  tax1 <- lapply(biom$rows, FUN=den)
  tax2 <- lapply(biom$rows, FUN=taxon)
  
  tax <- data.frame(den.otu=unlist(tax1), taxonomy=unlist(tax2))
  
  comm <- as.data.frame(comm) 
  
  attach(comm)
  comm <- comm[order(row.names(comm)),]
  detach(comm)
  
  comm <- t(comm) %>% as.data.frame()
  
  row.names(comm) <- tax1
  
  write.csv(comm, "./feature_table.csv")
  write.csv(tax, "taxonomy.csv")
}