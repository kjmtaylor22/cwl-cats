## The following functions were written by Kara JM Taylor for use with 
## 16S rRNA metabarcoding data, as output in table format by Qiime 1.8-1.9.
## It was written to produce the data presented in Ngunjiri et al. (2019). 
## "Longitudinal dynamics and interaction of respiratory and gut microbiota in 
## commercial chicken layers across all stages of the farm sequence."

## To perform properly, the functions must be used with Qiime formatted data,
## according to the pipeline given in CWL85_pipeline.R: OTUs must have an 
## associated taxonomy, either from the GreenGenes database, or from the Silva 
## release 132 classifier for Qiime. Other table formats and taxonomic classifications
## from other versions of these databases may not function as expected. 

## Please read the README.txt file prior to use. Descriptions of the functions
## and their options are given in README.txt. An example pipeline is given in 
## the CWL85_pipeline.R file also in this directory.

## These functions may need to be modified to perform correctly with other data.
## Any questions, feel free to contact Kara at taylor.1895@osu.edu

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





bact.tax <- function(taxonomy, database=NULL){
  if (!"dplyr"%in%row.names(installed.packages())) {
    install.packages("dplyr")
  }
  require(dplyr)
  
  splt <- function(x){
    y <- strsplit(x, split="; ", fixed=TRUE)
  }
  
  trim <- function(x){ #where x is a character string
    if ( is.null(database)==TRUE ){
      stop("Please specify database: SILVA (`silva`) or GreenGenes (`green`)")
    }
    if ( !(database%in%c("silva","green")) ){
      stop("Please specify database: SILVA (`silva`) or GreenGenes (`green`)")
    }
    if (database=="silva"){
      y <- substring(x, first=6)
    }
    if (database=="green"){
      y <- substring(x, first=4)
    } 
    return(y)
  }
  
  silv <- function(x){
    require(stringr)
    
    for (w in c("^un", "^gut", "^meta", " bacterium", "taxa", "taxon", "[-_()]", "group", "sensu stricto")){
      w1 <- grep(w, x)
      if (w %in% c("^un", "^gut", "^meta", " bacterium", "taxa", "taxon")){
        x[w1] <- ""
      }
      if (w=="[-_()]" & length(w1)>0){
        w2 <- grep("^[A-Z]", x[7])
        if (length(w2)==1){
          if (w1==7){
            x[7] <- ""
          } else {
            x[w1] <- strsplit(x[7], split=" ")[[1]][1]
          }
        } 
        if (length(w1)>=1 & length(w2)==0){
          a1 <- grep("Escherichia", x[6])
          a2 <- grep("Shigella", x[6])
          if (length(a1)==0 | length(a2)==0){
            x[w1] <- ""
          }
        }
      }
      
      if (w=="group" & length(w1)>0){
        for (i in w1){
          w2 <- strsplit(x[i], split=" ")
          w3 <- grep("[A-Z]", w2[[1]])
          if (length(w3) > 1){
            x[i] <- ""
          }
          if (length(w3)==1){
            x[i] <- w2[[1]][1]
          }
        }
      }
      if (w=="sensu stricto" & length(w1)>0){
        for (i in w1){
          x[i] <- strsplit(x[i], split=" ")[[1]][1]
        }
      }
    }
    
    if (length(x)==7 & x[7]!=""){
      w2 <- strsplit(x[7], split=" ")
      if (length(w2[[1]])==2){
        x[7] <- w2[[1]][2]
      } 
      if (length(w2[[1]])>2){
        x[7] <- ""
      }
      if (x[7]=="sp."){
        x[7] <- ""
      }
    }
    
    if (length(x)==7 & x[6]==""){
      x[7] <- ""
    }
    
    
    s <- str_count(x, "[0-9]")
    if (length(s) > 0){
      for (i in 1:length(s)){
        if (s[i]==1){
          x[i] <- strsplit(x[i], split=" ")[[1]][1]
        } else {next}
      }
    }
    
    if (length(x) > 5){ 
      w1 <- str_count(x[6], "[A-Z]")
      if (w1==2 & str_count(x[6], "[ ]")==1){
        w3 <- strsplit(x[6], split=" ")[[1]]
        x[6] <- w3[1]
        x[7] <- w3[2] #tolower(w3[2])
      }
    }
    return(x)
  }
  
  add <- function(x){
    if (length(x)==7){
      y <- x
    } else {
      y <- append(x, rep("", 7-length(x)))
    }
  }
  
  org <- function(x){
    w <- grep("[1-9]", x)
    if (length(w)==0){
      y <- max(which(x!=""))
      if (y==7){
        z <- paste(x[6], x[7], sep=" ")
      } else {
        z <- paste(x[y])
      }
    } else {
      #z <- paste(x[min(w)-1], paste(x[w], collapse=" "), sep=" ")
      z <- x[min(w)-1]
    }
    pull <- grep("\\[", z)
    if (length(pull) > 0){
      ext <- unlist(strsplit(z, split=""))
      p1 <- grep("\\[", ext)
      p2 <- grep("\\]", ext)
      ext <- ext[-c(p1, p2)]
      z <- paste(ext,collapse="")
    }
    return(z)
  }
  
  tidy <- function(x){
    pull <- grep("\\[", x, value=TRUE)
    get <- grep("\\[", x)
    if (length(pull)==0){
      return(x)
    } else {
      for (i in length(pull)){
        ext <- unlist(strsplit(pull[i], split=""))
        p1 <- grep("\\[", ext)
        p2 <- grep("\\]", ext)
        ext <- ext[-c(p1, p2)]
        z <- paste(ext, collapse="")
        x[get] <- z              
      }
      return(x)
    }
  }
  
  single <- function(x){
    y <- strsplit(as.character(x), split=" ")
    y1 <- lapply(y, FUN=length)
    y2 <- which(y1==1)
  }
  
  separated <- sapply(as.character(taxonomy[,2]), FUN=splt)
  trimmed <- lapply(separated, FUN=trim)
  
  if (database=="silva"){
    trimmed <- lapply(trimmed, FUN=silv)
  }
  
  evened <- sapply(trimmed, FUN=add)
  
  named <- apply(evened, MARGIN=2, FUN=org)
  names(named) <- NULL
  named <- unlist(named)
  
  out <- data.frame(den.otu=taxonomy[,1], taxonomy=named)
  out <- mutate(out, tag=paste("otu", row.names(out), sep=""))
  
  y <- rep(NA, dim(out)[1])
  y[single(out[,2])] <- paste(out[single(out[,2]),2], " (", 
                              out[single(out[,2]),3], ")", sep="")
  y[-single(out[,2])] <- paste(out[-single(out[,2]),2])
  
  
  out <- data.frame(out, otu.name=y)
  
  z <- transmute(out, tag.name=paste0(otu.name, " (", tag, ")")) 
  z$tag.name[grep("[(]", y)] <- y[grep("[(]", y)]
  
  out <- data.frame(out, tag.name=z$tag.name)
  
  
  evened <- apply(evened, MARGIN=2, FUN=tidy)
  classified <- as.data.frame(t(evened), row.names=1:dim(evened)[2])
  names(classified) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  out <- data.frame(out, classified)
  
  
  #if genus-species not available, use smallest level possible (but not strain codes; 
  #use regular expressions to pull out any names given as alphanumeric codes)
  #and then follow that name with either "sp" or the strain code
  tax <- out
  save(tax, file="./taxonomy.RD")
  
  return(out)
}





pcoa.plotting <- function(dist, meta, group, colors="rainbow", method="", axes=c(1,2,3), fixed=NULL){
  
  
  pkgs <- c("dplyr", "ape", "ggplot2", "gridExtra", "RColorBrewer", "pairwiseAdonis")
  if (any(!pkgs%in%row.names(installed.packages()))) {
    if (!"pairwiseAdonis"%in%row.names(installed.packages())){
      devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
    }
    install.packages(pkgs[which(!pkgs%in%row.names(installed.packages())==T)])
  }
  
  require(dplyr)
  require(ape)
  require(ggplot2)
  require(gridExtra)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/pcoa/")==FALSE){
    dir.create("./coreplots/pcoa/")
  }
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    names(meta)[1] <- "sample"
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
  }
  
  
  g <- which(names(meta)==group)
  
  if (class(dist)=="dist"){
    twPCoA <- pcoa(dist, rn=labels(dist))
  }
  if (class(dist)%in%c("matrix", "data.frame")){
    twPCoA <- pcoa(dist, rn=row.names(dist))
  }  
  
  if (any(!row.names(twPCoA)%in%meta[,which(names(twPCoA)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }  
  
  eig <- twPCoA$values$Eigenvalues[1:length(twPCoA$vectors[1,])] #only positive eigenvalues
  
  rescale <- eig/sum(eig)
  varexp <- round(rescale/sum(rescale)*100, 1)
  twPCOA <- data.frame(twPCoA$vectors)
  new_names <- rep("", ncol(twPCOA))
  for(i in 1:ncol(twPCOA)){
    new_names[i] <- paste("PCo",i, sep="")
  }
  names(twPCOA) <- new_names
  
  twPCOA <- left_join(data.frame(twPCOA, sample=row.names(twPCOA)), 
                      meta[,c(1,g)])
  
  
  
  if (colors=="rainbow"){
    palette(c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet"))
    if (length(unique(as.factor(meta[,g]))) > 9){
      palette(colorRampPalette(palette())(length(unique(as.factor(twPCOA[,g2])))))
    }
  }
  if (colors=="blind"){
    palette(c("#0072B2", "#D55E00", "#CC79A7", "#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442"))
    if (length(unique(twPCOA[,dim(twPCOA)[2]])) > 8){
      palette(colorRampPalette(palette())(length(unique(as.factor(twPCOA[,g2])))))
    } 
  }
  if (colors=="Dark2"){
    palette(RColorBrewer::brewer.pal(8, "Dark2"))
    if (length(unique(twPCOA[,dim(twPCOA)[2]])) > 8){
      palette(colorRampPalette(palette())(length(unique(as.factor(twPCOA[,g2])))))
    }  
  }
  if (colors=="Set1"){
    palette(RColorBrewer::brewer.pal(8, "Set1"))
    if (length(unique(twPCOA[,dim(twPCOA)[2]])) > 8){
      palette(colorRampPalette(palette())(length(unique(as.factor(twPCOA[,g2])))))
    }
  }
  if (colors=="reds"){
    palette(RColorBrewer::brewer.pal(10, "RdGy"))
    palette(colorRampPalette(palette())(length(unique(twPCOA[,dim(twPCOA)[2]]))))
  }
  
  if (is.null(fixed)){
    yl <- expression(guides())
    xl <- expression(guides())
  } else {
    yl <- expression(ylim(fixed[[2]][1], fixed[[2]][2]))
    xl <- expression(xlim(fixed[[1]][1], fixed[[1]][2]))
  }
  
  pch <- c(0,1,2,3,5,6,15,16,17,18)
  shps <- expression(scale_shape())
  
  if (class(with(meta, eval(parse(text=group))))=="factor"){
    grad <- expression(scale_fill_manual(values=palette()[1:length(unique(twPCOA[,dim(twPCOA)[2]]))]))
    if (length(unique(twPCOA[,dim(twPCOA)[2]])) > 10){nr=8} else {nr=2}
    
    guids <- expression(guides(fill=guide_legend(nrow=nr, title=group, 
                                                 override.aes = list(size=5,
                                                                     shape=21,
                                                                     fill=palette()[1:length(unique(twPCOA[,dim(twPCOA)[2]]))],
                                                                     color="black"))))
  }
  if (class(with(meta, eval(parse(text=group))))=="numeric"){
    grad <- expression(scale_fill_gradientn(colors=c("blue", "green", "red"),
                                            guide=guide_colorbar(title=group,
                                                                 label.theme=element_text(angle=315))))
    
    guids <- expression(labs(caption=""))
  }
  
  x = twPCOA[,axes[1]]
  y = twPCOA[,axes[2]]
  z = twPCOA[,axes[3]]
  
  
  points1 <- expression(geom_point(size = 4, aes(x, y, fill = eval(parse(text=group))), shape=21))
  points2 <- expression(geom_point(size = 4, aes(x, z, fill = eval(parse(text=group))), shape=21))
  points3 <- expression(geom_point(size = 4, aes(y, z, fill = eval(parse(text=group))), shape=21))
  
  
  p12 <- ggplot(twPCOA) +
    eval(expr=points1) +
    geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
    geom_vline(xintercept=0, color="grey80", linetype="dotdash") +
    labs(title=paste("Beta diversity in", deparse(substitute(dist))), 
         subtitle=paste("Distance method: ", method)) +
    xlab(names(twPCOA)[axes[1]]) +
    ylab(names(twPCOA)[axes[2]]) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(), 
          panel.background = element_rect(fill = "white", color = "black", size = 1),
          plot.title = element_text(size = 18, face="bold", hjust = 0.5),
          axis.title = element_text(size = 15, face="bold"),
          axis.text = element_text(color = "black", size = 12),
          legend.key = element_rect(fill="white"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12, face="bold"),
          legend.position="bottom",
          legend.justification = "center") +
    eval(expr=grad) +
    eval(expr=shps) +
    eval(expr=guids) +
    labs(caption=paste("Variance explained: ", "\n",
                       paste(names(twPCOA)[axes[1]],"(", as.character(varexp[axes[1]]), "%)"), "\n",
                       paste(names(twPCOA)[axes[2]],"(", as.character(varexp[axes[2]]), "%)"), sep=""))
  
  p13 <- ggplot(twPCOA) +
    eval(expr=points2) +
    geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
    geom_vline(xintercept=0, color="grey80", linetype="dotdash") +
    labs(title=paste("Beta diversity in", deparse(substitute(dist))), 
         subtitle=paste("Distance method: ", method)) +
    xlab(names(twPCOA)[axes[1]]) +
    ylab(names(twPCOA)[axes[3]]) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(), panel.background = element_rect(fill = "white", color = "black", size = 1),
          plot.title = element_text(size = 18, face="bold", hjust = 0.5),
          axis.title = element_text(size = 15, face="bold"),
          axis.text = element_text(color = "black", size = 12),
          legend.key = element_rect(fill="white"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12, face="bold"),
          legend.position="bottom",
          legend.justification = "center") +
    eval(expr=grad) +
    eval(expr=shps) +
    eval(expr=guids) +
    labs(caption=paste("Variance explained: ", "\n",
                       paste(names(twPCOA)[axes[1]],"(", as.character(varexp[axes[1]]), "%)"), "\n",
                       paste(names(twPCOA)[axes[3]],"(", as.character(varexp[axes[3]]), "%)"), sep=""))
  
  p23 <- ggplot(twPCOA) +
    eval(expr=points3) +
    geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
    geom_vline(xintercept=0, color="grey80", linetype="dotdash") +
    labs(title=paste("Beta diversity in", deparse(substitute(dist))), 
         subtitle=paste("Distance method: ", method)) +
    xlab(names(twPCOA)[axes[2]]) +
    ylab(names(twPCOA)[axes[3]]) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(), panel.background = element_rect(fill = "white", color = "black", size = 1),
          plot.title = element_text(size = 18, face="bold", hjust = 0.5),
          axis.title = element_text(size = 15, face="bold"),
          axis.text = element_text(color = "black", size = 12),
          legend.key = element_rect(fill="white"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12, face="bold"),
          legend.position="bottom",
          legend.justification = "center") +
    eval(expr=grad) +
    eval(expr=shps) +
    eval(expr=guids) +
    labs(caption=paste("Variance explained: ", "\n",
                       paste(names(twPCOA)[axes[2]],"(", as.character(varexp[axes[2]]), "%)"), "\n",
                       paste(names(twPCOA)[axes[3]],"(", as.character(varexp[axes[3]]), "%)"), sep=""))
  
  tiff(paste("./coreplots/pcoa/PCOA_", as.character(axes)[1], as.character(axes)[2], 
             as.character(axes)[3], paste(deparse(substitute(dist)), group, sep="_"), 
             ".tiff", sep=""),
       width=9000, height=4000, res=600)
  grid.arrange(p12, p13, p23, newpage=FALSE, nrow=1, ncol=3)
  dev.off()
  
  
  sink(paste("./coreplots/pcoa/PCOA_", as.character(axes)[1], as.character(axes)[2], 
             as.character(axes)[3], paste(deparse(substitute(dist)), group, sep="_"), 
             ".txt", sep=""))
  print(pairwiseAdonis::pairwise.adonis(dist, twPCOA[,dim(twPCOA)[2]]))
  sink()
  
}






core.id <- function(comm, meta, tax, group, subgroup, margin){ 
  # where x is a sample-by-OTU table only, rows as samples, columns as otus
  # Row.names should be the sample IDs.
  # Any metadata should be input as 'factors', 
  # and make sure that the ID column in the metadata is labeled "SampleID"
  if (!"dplyr"%in%row.names(installed.packages())) {install.packages("dplyr")}
  require(dplyr)
  
  den <- colnames(comm) 
  
  new <- data.frame(sample=row.names(comm), comm)
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
  }
  
  if (any(!row.names(comm)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  data <- left_join(new, meta)
  
  rm(new)
  colnames(data) <- c("ID", names(data)[-1])
  
  core.mb <- function(x){
    z1 <- length(x[which(x>0)])
    if (length(margin)==1){
      if (z1/length(x) >= margin){
        return(paste("+"))
      }
    } else {
      if (z1/length(x) >= margin[i]){
        return(paste("+"))
      }
    }
  }
  
  
  out <- as.list(unique(as.character(with(data, eval(parse(text=group))))))
  names(out) <- unique(as.character(with(data, eval(parse(text=group)))))
  
  g <- which(names(data)==group)
  sg <- which(names(data)==subgroup)
  
  for (i in 1:length(out)){
    sub1 <- data[which(data[,g]==out[[i]]),]
    grp <- unique(as.character(sub1[,sg]))
    out1 <- as.list(unique(as.character(sub1[,sg])))
    
    for (j in 1:length(grp)){
      if (as.character(grp[j])%in%unique(as.character(sub1[,sg]))){
        cat(paste(out[[i]]," : ", grp[j], "\n", sep=""))
        sub2 <- sub1[which(sub1[,sg]==grp[j]),
                     which(names(sub1)%in%den)]
        z <- apply(sub2, MARGIN=2, FUN=core.mb) %>% unlist()
        out1[[j]] <- append(out1[[j]], names(z))
        out1[[j]] <- out1[[j]][-1]
        names(out1)[[j]] <- as.character(grp[j])
        if (length(out1[[j]])==0){
          out1[[j]] <- NA
        }
      }
    }
    
    out[[i]] <- append(out[[i]], out1)
    out[[i]] <- out[[i]][-1]
  }
  return(out)
}






core.stack <- function(data, list, tax=NULL, factors, group, subgroup, hi.tax="tag", threshold=NULL, margin, 
                       date, core.only=FALSE, fixed=FALSE, landscape=TRUE, view.all=FALSE){
  
  cran <- c("dplyr", "RColorBrewer", "ggplot2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  require(dplyr)
  require(ggplot2)
  require(RColorBrewer)
  
  colors <- c("#330033", "#F0E442", "#0072B2", "#FF0033", "#FF9933", "#999999", "#56B4E9", "#FFCC00", "#00FF00", "#FF66FF", "#666666")
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/profiles/")==FALSE){
    dir.create("./coreplots/profiles/")
  } 
  
  fmx <- as.data.frame(t(data)) %>%
    summarize_all(funs(sum=sum)) %>%
    max()
  
  ra <- function(x){ #where x is a vector
    y <- sum(x)
    z <- sapply(x, FUN=function(x){z1 <- x/y})
    return(z)
  }
  avg <- function(x){
    z <- mean(x, na.rm=TRUE)
    return(z)
  }
  gcd <- function(x,y) {
    r <- x%%y
    return(ifelse(r, gcd(y, r), y))
  }
  
  id <- grep("sample", names(factors), ignore.case = T)
  if (length(id)==1){
    names(factors)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(factors)[1], ") as sample ID. "))
    names(factors)[1] <- "sample"
  }
  
  if (any(!row.names(data)%in%factors[,which(names(factors)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  g <- which(names(factors)==group)
  sg <- which(names(factors)==subgroup)
  s <- which(names(factors)=="sample")
  features <- colnames(data)
  
  
  new <- data.frame(ID=row.names(data), data)
  data <- right_join(factors[,c(s,g,sg)], new, by=c("sample"="ID"))
  data <- mutate(data, supergroup=paste(eval(parse(text=group)), 
                                        eval(parse(text=subgroup)), sep="."))
  
  gsg <- unique(data[,c(2,3,which(names(data)=="supergroup"))])
  
  cat("Calculating as relative abundance (%) \n")
  
  data <- data[,-c(1:3)] %>%
    group_by(supergroup) %>%
    summarize_all(funs(sum=sum))
  
  supergroup <- data[,1]
  
  
  data <- apply(data[,-1], MARGIN=1, FUN=ra) %>% t() %>% data.frame()
  
  
  names(data) <- features
  data <- cbind(supergroup, data)
  
  if (is.null(threshold)==TRUE){
    t.av <- function(x){
      y <- which(x>0)
      z <- mean(x[y], na.rm=TRUE)
      return(z)
    }
    new <- apply(data, MARGIN=1, FUN=t.av)
    split <- strsplit(names(new), split=NULL)
    splt <- function(x){
      y <- grep("[A-Z]", x)
      if (y[1]==y[length(y)]){
        z <- x[y[1]]
      } else {
        if (x[y[1]]=="B"){
          z <- x[y[1]]
        } else {
          z <- paste(x[y[1]], x[y[length(y)]], sep="")
        }
      }
      return(z)
    }
    fromgroup <- lapply(split, splt) %>% as.data.frame() %>% t()
    new <- data.frame(new, fromgroup=fromgroup)
    new <- summarize(group_by(new, fromgroup), mean=mean(new))
    thresh <- new$mean
    rm(new, fromgroup, split)
  }
  
  data <- right_join(gsg, data, by=c("supergroup"="supergroup"))
  
  g <- which(names(data)==group)
  sg <- which(names(data)==subgroup)
  
  df <- data.frame()
  for (i in 1:length(list)){
    if (exists("thresh")){
      threshold <- thresh[i]
      threshold <- round(threshold, digits=4)
    }
    
    cat(paste("Calculating core microbiome for ", names(list)[i], "\n", sep=""))
    group.core <- data.frame()
    
    if (all(is.na(unlist(list[[i]])))){next}
    
    
    
    core.pool <- data.frame(count=rep(1, length(unlist(list[i]))),
                            otu=unlist(list[i]))
    core <- names(which(summary(as.factor(core.pool$otu))==length(which(is.na(list[[i]])==FALSE))))
    
    
    if (length(core)>0 & length(which(is.na(list[[i]])==FALSE))>=1){
      sub <- data[which(data[,g]==names(list[i])),
                  which(names(data)%in%c(subgroup, core))]
      
      core <- names(sub)[which(names(sub)%in%core)]
      
      if (is.null(tax)==FALSE){
        core <- data.frame(v1=c(1:length(core)), v2=core)
        core <- left_join(core, tax, by=c("v2"="tag"))
        
        if (hi.tax%in%names(tax)){
          get <- which(names(core)==hi.tax)
          core <- as.character(unlist(core[,get]))
        } else {core <- as.character(core$otu.name)}
        
        grouped <- data.frame(features=core, t(sub[,-1]))
        
        grouped <- grouped %>%
          group_by(features) %>%
          summarize_all(funs(sum=sum))
        
        names(grouped) <- c("features", unlist(as.character(sub[,1])))
        
        core <- as.character(unlist(grouped$features))
        grouped <- grouped[,-1]
        grouped <- t(grouped) %>% as.data.frame()
        names(grouped) <- core
        
        grouped <- cbind(v1=as.character(sub[,1]), grouped)
        
        names(grouped) <- c(subgroup, names(grouped)[-1])
        row.names(grouped) <- grouped[,1]
      }
      
      if (is.null(tax)==TRUE){
        grouped <- sub
      }
      
      
      grouped <- reshape(grouped, varying=core,
                         v.names="rel.abund",
                         timevar="otu",
                         times=core,
                         direction="long")
      
      grouped <- data.frame(grouped, groupvar=rep(names(list)[i], dim(grouped)[1]))
      names(grouped) <- c(subgroup, names(grouped)[-1])
      
      grouped$otu <- as.factor(grouped$otu)
      grouped$otu <- factor(grouped$otu, levels = rev(levels(grouped$otu)))
      pull <- grouped$otu[which(grouped$rel.abund > threshold)]
      pull <- pull[which(duplicated(pull)==FALSE)]
      grouped <- grouped[which(grouped$otu%in%pull),]
      group.core <- cbind(grouped, core=rep("core", dim(grouped)[1]))
    }
    
    if (core.only==FALSE){
      for (j in 2:length(list[[i]])){
        if (is.na(list[[i]][j])==FALSE){
          pull <- which(unlist(list[[i]][j:length(list[[i]])])%in%unlist(list[[i]][1:j-1])==FALSE)
          pull <- unlist(list[[i]][j:length(list[[i]])], use.names = FALSE)[pull]
          pool <- data.frame(count=rep(1, length(pull)), otu=pull)
          
          otus <- names(which(summary(pool$otu)==
                                length(which(is.na(list[[i]][j:length(list[[i]])])==FALSE))))
          
          if (length(otus)>0 & length(which(is.na(list[[i]])==FALSE))>=1){
            sub <- data[which(data[,g]==names(list[i])),
                        which(names(data)%in%c(subgroup, otus))]
            
            if (is.null(tax)==FALSE){
              otus <- data.frame(v1=c(1:length(otus)), v2=otus)
              otus <- left_join(otus, tax, by=c("v2"="tag"))
              
              if (hi.tax%in%names(tax)){
                get <- which(names(otus)==hi.tax)
                otus <- as.character(unlist(otus[,get]))
              } else {otus <- as.character(otus$otu.name)}
              
              grouped <- data.frame(features=otus, t(sub[,-1]))
              
              grouped <- grouped %>%
                group_by(features) %>%
                summarize_all(funs(sum=sum))
              
              otus <- as.character(unlist(grouped$features))
              grouped <- grouped[,-1]
              grouped <- t(grouped) %>% as.data.frame()
              names(grouped) <- otus
              
              
              grouped <- cbind(v1=as.character(sub[,1]), grouped)
              
              names(grouped) <- c(subgroup, names(grouped)[-1])
              row.names(grouped) <- grouped[,1]
            }
            if (is.null(tax)==TRUE){
              grouped <- sub
            }
            
            grouped <- reshape(grouped, varying=otus,
                               v.names="rel.abund",
                               timevar="otu",
                               times=otus,
                               direction="long")
            
            grouped <- data.frame(grouped, groupvar=rep(names(list)[i], dim(grouped)[1]))
            names(grouped) <- c(subgroup, names(grouped)[-1])
            grouped$otu <- as.factor(grouped$otu)
            grouped$otu <- factor(grouped$otu, levels = rev(levels(grouped$otu)))
            pull <- grouped$otu[which(grouped$rel.abund > threshold)]
            
            
            
            if(length(pull)>0){
              pull <- pull[which(duplicated(pull)==FALSE)]
              grouped <- grouped[which(grouped$otu%in%pull),]
              grouped <- cbind(grouped, core=paste(names(list[[i]][j]), "-->", names(list[[i]][length(list[[i]])]), sep=""))
              if (dim(grouped)[1]==0){next}
            }
            if (length(pull)==0){
              grouped <- group.core[1,]
              grouped$rel.abund[1] <- 0
              grouped$core <- paste(names(list[[i]][j]), "-->", names(list[[i]][length(list[[i]])]), sep="")
              if (dim(grouped)[1]==0){next}
            }
            group.core <- rbind(group.core, grouped)
          } else {
            grouped <- group.core[1,]
            grouped$rel.abund[1] <- 0
            grouped$core <- paste(names(list[[i]][j]), "-->", names(list[[i]][length(list[[i]])]), sep="")
            if (dim(grouped)[1]==0){next}
            group.core <- rbind(group.core, grouped)
          }
        } else {
          grouped <- group.core[1,]
          grouped$rel.abund[1] <- 0
          grouped$core <- paste(names(list[[i]][j]), "-->", names(list[[i]][length(list[[i]])]), sep="")
          if (dim(grouped)[1]==0){next}
          group.core <- rbind(group.core, grouped)
        }
      }
    }
    df <- rbind(df, group.core)
    
    write.csv(group.core, paste("./coreplots/profiles/", names(list)[i], date, "_", 
                                threshold, "_", as.character(margin), ".csv", sep=""))
    
    
    
    face <- c(floor(length(unique(group.core$core))/3), ceiling(length(unique(group.core$core))/3))
    r <- expression(gcd(h, w))
    kh <- 1
    
    if (face[2]>face[1]){
      if (length(unique(group.core$otu)) > 25){
        cols <- floor(length(unique(group.core$otu))/25)
      } else {cols <- 1}
      leg <- list(c(0.833, 0.22), "center")
    } else {
      cols <- 4
      leg <- list("bottom", "center")
    }
    
    if (dim(group.core)[1]==0){next}
    group.core <- group.core[which(group.core$core%in%names(summary(group.core$core)[which(summary(group.core$core)!=1)])),]
    group.core$core <- factor(group.core$core)
    
    h <- 800*length(unique(group.core$core))+1200
    w <- 1000*(ceiling(length(levels(group.core[,1]))/2))+800
    
    fexp <- expression(facet_grid(core~., drop=TRUE))
    fexp2 <- expression(facet_grid(groupvar~core, drop=FALSE))
    
    leg <- list("right", "right")
    
    if (core.only==FALSE){
      if (length(unique(group.core$otu)) > 40){
        kh <- floor(length(unique(group.core$otu))/40)
      } 
      if (length(unique(group.core$otu))/length(unique(group.core$core)) > 4){
        kh <- length(unique(group.core$otu))/4
      }
    }
    cols <- 1
    
    
    
    
    if (fixed==TRUE){
      fix.ylim <- expression(ylim(0,100))
    } else {
      fix.ylim <- expression(ylim(0, NA))
    }
    
    fs <- w/gcd(h, w)
    fs <- round(c(fs+(fs/2), fs+(fs/3.2), fs+(fs/4), fs+(fs/8), fs+(fs/16), fs*2))
    
    group.core$rel.abund <- group.core$rel.abund*100
    
    lev <- match(levels(with(meta, eval(parse(text=subgroup)))),
                 unique(with(group.core, eval(parse(text=subgroup)))))
    
    group.core[,which(names(group.core)==subgroup)] <-
      factor(group.core[,which(names(group.core)==subgroup)],
             levels=levels(with(group.core, eval(parse(text=subgroup))))[
               lev[which(is.na(lev)==FALSE)]])
    
    if (view.all==T){
      fexp <- expression(facet_grid(groupvar~., drop=TRUE))
      fexp2 <- expression(facet_grid(groupvar~., drop=TRUE))
      h <- 2000
      w <- 1000*(ceiling(length(levels(group.core[,1]))/2))+800
    }
    
    f <- ggplot(group.core, aes(eval(parse(text=subgroup)), rel.abund, fill=otu)) +
      eval(expr=fexp) +
      geom_bar(stat="identity") +
      eval(expr=fix.ylim) +
      theme_minimal() +
      scale_fill_manual(values=
                          colorRampPalette(colors)(length(unique(as.factor(group.core$otu))))) +
      theme(panel.border=element_rect(linetype="solid", fill=NA),
            legend.position=leg[[1]],
            legend.justification=leg[[2]],
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = fs[3], face="bold"),
            plot.title = element_text(size = fs[1], face="bold"),
            axis.title = element_text(size = fs[2], face="bold"),
            axis.text.x = element_text(color = "black", size = fs[6]),
            axis.text.y = element_text(color = "black", size = fs[2], face="bold"),
            legend.title = element_text(size = fs[3], face="bold"),
            legend.text = element_text(size = fs[3]),
            plot.caption = element_text(size=fs[5]),
            legend.background=element_rect(fill="white", colour="black"))
    
    jpeg(file=paste("./coreplots/profiles/", names(list[i]), date, "_", threshold, "_", as.character(margin), ".jpeg", sep=""),
         width=w, height=h, res=eval(expr=r), units="px")
    print(f + labs(title=paste("Progression of Core Microbiome Composition in", names(list)[i])) + 
            xlab(subgroup) + ylab("Relative Abundance (%)") +
            guides(fill=guide_legend(keywidth = fs[2]/10, keyheight=fs[2]/(10*kh), ncol=cols)) +
            labs(caption=paste("Occurence margin: ", as.character(margin*100), "%", "\n",
                               "Abundance threshold: ", as.character(threshold*100), "%")))
    dev.off()
  }
  write.csv(df, paste("./coreplots/profiles/facet", date, "_",
                      threshold, "_", as.character(margin), ".csv", sep=""))
  if ("BG" %in% unique(df$groupvar)){
    df <- subset(df, subset=groupvar!="BG")
    df$core <- factor(df$core)
  }
  
  sg <- which(names(df)==subgroup)
  pull <- unique(data.frame(df$groupvar, df[,sg]))
  pull <- names(which(summary(pull[,1])==length(unique(df[,sg]))))
  df <- df[which(df$groupvar%in%pull),]
  df$groupvar <- factor(df$groupvar)
  df$rel.abund <- df$rel.abund*100
  if (landscape==FALSE){
    w <- 4000
    h <- 8000
    legs <- expression(theme(legend.position="right",
                             legend.justification="right"))
    gui <- expression(guides(fill=guide_legend(keywidth = 1, keyheight=0.75, ncol=1)))
    
    
  } else {
    w <- 8000
    h <- 4000
    legs <- expression(theme(legend.position="bottom",
                             legend.justification="center"))
    gui <- expression(guides(fill=guide_legend(keywidth = 1, keyheight=0.75, ncol=7)))
    
  }
  
  df[,which(names(df)==subgroup)] <-
    factor(df[,which(names(df)==subgroup)],
           levels=levels(with(df, eval(parse(text=subgroup))))[
             lev[which(is.na(lev)==FALSE)]])
  
  bff <- ggplot(df, aes(eval(parse(text=subgroup)), rel.abund, fill=otu)) +
    geom_bar(stat="identity") +
    eval(expr=fix.ylim) +
    eval(expr=fexp2) +
    theme_minimal() +
    scale_fill_manual(values=
                        colorRampPalette(colors)(length(unique(as.factor(df$otu))))) +
    theme(panel.border=element_rect(linetype="solid", fill=NA),
          panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour="grey90", linetype="36"),
          panel.grid.minor.y = element_blank(),
          axis.text.x=element_text(angle=330),
          legend.background=element_rect(fill="white", colour="black")) +
    eval(expr=legs)
  jpeg(file=paste("./coreplots/profiles/facet", date, "_", threshold, "_", as.character(margin), ".jpeg", sep=""),
       width=w, res=600, units="px", height=h) #height=h*length(list))
  print(bff + labs(title=paste("Progression of Core Microbiome Composition")) + 
          xlab(subgroup) + ylab("Relative Abundance (%)") +
          eval(expr=gui) +
          labs(caption=paste("Occurence margin: ", as.character(margin*100), "%", "\n",
                             "Abundance threshold: ", as.character(threshold*100), "%")))
  dev.off()
  return(df)
}






dendro.heatmap <- function(comm, tax, meta, path="", group, subgroup, core.list=NULL, 
                           hi.tax="tag", file="", classification=NULL, trans="log", stop="no"){
  #### Series of nested functions to produce grouped heatmap
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
  }
  
  if (any(!row.names(comm)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/profiles/")==FALSE){
    dir.create("./coreplots/profiles/")
  } 
  
  cran <- c("dplyr", "egg", "ape", "vegan", "ggplot2", "reshape2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  biocon <- c("ggtree", "ggdendro", "dendextend")
  if (any(!biocon%in%row.names(installed.packages()))){
    source("https://bioconductor.org/biocLite.R")
    for (i in biocon){
      biocLite(i, suppressUpdates=T)
    }
  }
  
  require(dplyr)
  require(egg)
  require(ape)
  require(vegan)
  require(ggplot2)
  require(ggtree)
  require(ggdendro)
  require(dendextend)
  
  sub.tax <- function(comm, tax, meta, group, subgroup,
                      classification=NULL, core.list=NULL, rank="log"){
    
    if (is.null(core.list)==FALSE){
      if (class(core.list)=="numeric"){
        otus <- names(sort(colMeans(comm), decreasing = T))[1:core.list]
      } else {
        test <- grep("^otu", core.list$otu)
        test2 <- grep("otu", core.list$otu)
        
        if (length(test) < length(core.list$otu) & length(test2) > 0){
          if (length(test2)==length(core.list$otu)){
            if (dim(tax)[2]>12){
              otus <- tax$tag[which(tax$tag.name%in%as.character(unique(core.list$otu)))]
            } else {otus <- tax$tag[which(tax$tag.name%in%as.character(unique(core.list$otu)))]}
          } else {
            otus <- tax$tag[which(tax$otu.name%in%as.character(unique(core.list$otu)))]
          }
          otus <- otus[which(duplicated(otus)==FALSE)]
          otus <- otus[which(otus %in% names(comm))]
        } 
        if (length(test)==length(test2)){
          otus <- as.character(unique(core.list$otu))
          otus <- otus[which(otus %in% names(comm))]
        } 
        if (length(test2)==0 & is.null(classification)==TRUE){
          for (i in 1:dim(tax)[2]){
            t <- all(core.list$otu %in% tax[,i])
            if (t==TRUE){break}
          }
          otus <- tax$tag[which(tax[,i] %in% as.character(unique(core.list$otu)))]
          otus <- otus[which(otus %in% names(comm))]
        }
        if (length(test2)==0 & is.null(classification)==FALSE){
          cat(" `core.list` providing taxonomic level; \n Please set `classification` to NULL")
          return(NULL)
        }
        comm <- comm[,which(names(comm) %in% otus)]
      }
    } else {otus <- colnames(comm)}
    
    data <- comm
    
    if (is.null(classification)==FALSE){
      data <- t(data)
      data <- data.frame(feature=row.names(data), data)
      
      get <- which(names(tax)==classification)
      data <- right_join(tax[,c(3,get)], data, by=c("tag"="feature"))
      
      data <- data[,-1] %>%
        group_by_(classification) %>%
        summarize_all(funs(sum=sum)) %>%
        as.data.frame()
      
      cols <- unlist(data[,1])
      
      data <- data[,-1]
      data <- t(data) %>% as.data.frame()
      row.names(data) <- row.names(comm)
      names(data) <- cols
      
      data <- data[,which(names(data)!="")]
      
      otus <- tax$tag[match(names(data), tax[,get])]
      names(data) <- otus
    }
    
    sub <- data
    
    if (trans=="none"){
      rank1 <- sub
    }
    if (trans=="rank"){
      rank1 <- rank(unlist(sub))
      try <- matrix(rank1, nrow=dim(sub)[1], ncol=length(otus), 
                    dimnames=list(row.names(sub), otus))
      rank1 <- as.data.frame(try)
    }
    
    if (trans=="log"){
      rank1 <- log10(sub+1)
    }
    
    
    sub <- data.frame(ID=row.names(rank1), rank1)
    sub <- right_join(meta, sub, by=c("sample"="ID"))
    sub <- mutate(sub, !!"group.subgroup" := 
                    paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep="."))
    sub <- subset(sub, select=c("group.subgroup", otus))
    
    row.names(sub) <- row.names(comm)
    top <- sub
    return(top)
    
  }
  
  
  
  
  ## function for creating cluster dendrogram for the groups and subgroups
  group.dendro <- function(data, meta, group, subgroup){
    
    
    g <- which(names(meta)==group)
    sg  <- which(names(meta)==subgroup)
    
    columns <- colnames(data)
    rows <- row.names(data)
    
    new <- data.frame(ID=row.names(data), data)
    
    data <- left_join(new, meta[,c(1, g, sg)], by=c("ID"="sample"))
    row.names(data) <- row.names(new)
    data <- data[,-1]
    rm(new)
    
    data <- mutate(data, group.subgroup=paste(with(data, eval(parse(text=group))), 
                                              with(data, eval(parse(text=subgroup))), sep="."))
    
    otus <- names(data)[grep("^otu", names(data))]
    n <- length(otus)
    
    data <- data[,c(which(names(data)%in%c("group.subgroup", otus)))]
    row.names(data) <- rows
    
    g.sg <- names(data)[-grep("^otu", names(data))]
    ngsg <- length(unique(with(data, eval(parse(text=g.sg)))))
    
    
    
    avgd <- data %>%
      group_by_(g.sg) %>% 
      summarize_all(funs(sum=sum)) %>%
      as.data.frame()
    row.names(avgd) <- avgd[,1]
    avgd <- subset(avgd, select=grep("^otu", names(avgd)))
    names(avgd) <- otus
    
    bdiv <- vegdist(avgd, "bray")
    
    clust <- hclust(d=bdiv, method="average")
    
    dend <- as.dendrogram(clust) %>%
      dendextend::rotate(sort(as.character(unique(data$group.subgroup))))
    
    ddend <- dendro_data(dend, type="rectangle")
    
    order <- ddend$labels$label
    zoom <- ngsg*(0.0595-0.0001*ngsg)
    
    d <- ggplot(segment(ddend)) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=2) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,0,0, "pt")) +
      coord_cartesian(xlim=c(zoom, ngsg+1-zoom),
                      ylim=c(max(ddend$segments$yend)/10, 
                             max(ddend$segments$yend)))
    d[["tip.labels"]] <- order
    
    return(d)
  }
  
  
  ## function for creating phylogenetic cladogram for top taxa
  otu.phylo <- function(top, tax, path, hi.tax){#tax= output from bact.tax
    
    g.sg <- names(top)[-grep("^otu", names(top))]
    otus <- names(top)[grep("^otu", names(top))]
    ngsg <- length(unique(with(top, eval(parse(text=g.sg)))))
    
    tree <- ape::read.tree(path)
    
    tips <- as.character(tax$den.otu[-match(otus, tax$tag)])
    
    prune <- match(tax$den.otu, tree$tip.label)
    
    if (length(prune)>0){
      if(any(is.na(prune))){
        prune <- prune[which(!is.na(prune))]
      }
      tips <- c(tips, tree$tip.label[-prune])
    }
    
    drop.tree <- ape::drop.tip(tree, tip=tips)
    
    labs <- data.frame(tips=drop.tree$tip.label)
    labs <- left_join(labs, tax, by=c("tips"="den.otu"))
    
    if (hi.tax%in%names(tax)){
      get <- which(names(labs)==hi.tax)
      drop.tree$tip.label <- as.character(unlist(labs[,get]))
    } else {
      drop.tree$tip.label <- as.character(labs$otu.name)
    }
    
    
    p <- ggtree::ggtree(drop.tree, branch.length = "none", ladderize=FALSE) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,0,0, "pt"))
    
    pull <- data.frame(node=p$data$node, label=p$data$label)
    dup <- which(duplicated(pull$label[1:length(otus)])==TRUE)
    new <- drop.tip(drop.tree, tip=pull$node[dup])
    
    
    n <- length(new$tip.label)
    
    if (n<150){
      zoom <- n*(0.0595-0.0001*n)
    }
    if (n>150){
      zoom <- n*(0.0595-0.000045*n)
    }
    
    cp <- ggtree::ggtree(new, branch.length = "none", size=2) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,0,0, "pt")) +
      coord_cartesian(ylim=c(zoom, n+1-zoom))
    
    new <- cp$data[1:n, c(4,5,7)]
    
    out <- new[order(new$y),]
    
    
    cp[["new.labels"]] <- out$label
    
    return(cp)
  }
  
  ## function for creating heatmap using the top n otus
  top.heat <- function(top, tax, phylo, dendro, hi.tax, file){ 
    
    g.sg <- names(top)[-grep("^otu", names(top))]
    otus <- names(top)[grep("^otu", names(top))]
    
    group <- top %>%
      group_by_(g.sg) %>%
      summarize_all(funs(mean=mean)) %>%
      as.data.frame()
    
    row.names(group) <- group[,1]
    group <- group[,-1]
    
    names(group) <- otus
    group <- group[match(dendro$tip.labels, row.names(group)),]
    
    group2 <- reshape(group, varying=otus,
                      v.names="abund",
                      timevar="otu",
                      times=otus,
                      direction="long",
                      ids=row.names(group))
    
    group2 <- left_join(group2, tax, by=c("otu"="tag"))
    
    if (hi.tax%in%names(tax)){
      get <- which(names(group2)==hi.tax)
      group3 <- group2[,c(2, 3, get)] %>%
        group_by_("id", hi.tax) %>%
        summarize_at("abund", funs(abund=max)) %>%
        as.data.frame()
      names(group3)[which(names(group3)==hi.tax)] <- "otu.name"
    } else {
      group3 <- group2[,c(2, 3, 6)] %>%
        group_by(id, otu.name) %>%
        summarize_at("abund", funs(abund=max)) %>%
        as.data.frame()
      
    }
    
    
    group3$id <- factor(group3$id,
                        levels=dendro$tip.labels)
    
    group3$otu.name <- factor(group3$otu.name,
                              levels=phylo$new.labels)
    write.csv(group3, paste0("./coreplots/profiles/heatmapdata_", file, ".csv"))
    
    h <- ggplot(group3, aes(id, otu.name)) +
      geom_raster(aes(fill=abund)) +
      scale_fill_gradientn(colors=rev(c("magenta", "red", "orange", "yellow", "green", "turquoise", "blue", "navy")), 
                           guide="colorbar") +
      theme(axis.ticks=element_line(color="white"), 
            axis.text.x=element_text(angle=90, size=60, hjust=1, vjust=0.5, face="bold"),
            axis.text.y=element_text(size=60, vjust=0.5),
            plot.margin = margin(1,1,5,1, "pt"),
            plot.caption=element_text(size=40),
            legend.title=element_text(size=45),
            legend.text=element_text(size=40),
            legend.key.size=unit(5, "lines")) +
      scale_y_discrete(position="right") +
      labs(caption=paste("Features: n = ", as.character(length(levels(group3$otu.name)))))
    
    return(h)
  }
  
  top <- sub.tax(comm, tax, meta, group, subgroup, classification, core.list, rank)
  if (stop=="top"){return(top)}
  den <- group.dendro(comm, meta, group, subgroup)
  if (stop=="den"){return(den)}
  phy <- otu.phylo(top, tax, path, hi.tax) #use hi.tax if you want core OR all taxa at level other than in `sub`
  if (stop=="phy"){return(phy)}
  heat <- top.heat(top, tax, phy, den, hi.tax, file) #use hi.tax if you want core OR all taxa
  
  e <- ggplot() + 
    geom_blank() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0, "pt"))
  
  n <- length(phy$new.labels)
  
  if (n>=50){
    ht <- 56*n
  } 
  if (n<50 & n>=25){
    ht <- 3000
  }
  if (n<25){
    ht <- 2000
  }
  
  hts <- c(0.3-0.00075*n, 3)
  wds <- c(0.5+0.003*n, 3)
  
  tiff(paste0("./coreplots/profiles/phyloheatmap_", file, ".tiff"), width=3600, height=ht)
  ggarrange(e, den, phy, heat, heights=hts, widths=c(.5,3))
  dev.off()
  
}





mbiom.venn <- function(mat, meta, group, tax=NULL, file="./", tax.grp="tag") {
  
  cran <- c("dplyr", "VennDiagram", "RColorBrewer")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  require(dplyr)
  require(VennDiagram)
  require(RColorBrewer)
  
  palette(c(RColorBrewer::brewer.pal(8, "Set1"), "#000000"))
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
  }
  
  if (any(!row.names(mat)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/venn/")==FALSE){
    dir.create("./coreplots/venn/")
  } 
  
  g <- which(names(meta)==group)
  
  mat2 <- data.frame(sample=row.names(mat), mat)
  mat <- right_join(meta[,c(1,g)], mat2)
  
  
  iff <- function(x){ 
    if (x > 0){
      x <- 1
    } else {
      x <- 0
    }
  }
  
  top <- function(x){
    y <- sum(x)
    return(y)
  }
  
  cumul <- function(x){
    y <- seq(1:length(x))
    z1 <- which(x==0)
    y[z1] <- 0
    x <- y
  }
  
  rem0 <- function(x){
    if (any(x==0)){
      x <- x[-which(x==0)]
    } else {
      x <- x
    }
  }
  
  
  df <- mat[,-1] %>%
    group_by_(group) %>%
    summarize_all(funs(mean=mean)) %>%
    as.data.frame()
  row.names(df) <- with(df, eval(parse(text=group)))
  
  test <- apply(df[,-1], MARGIN=1, FUN=cumul) %>% as.data.frame()
  
  renamex <- names(df[,-1]) 
  renamex <- unlist(strsplit(renamex, split="_"))[seq(1,2*length(renamex), 2)]
  renamex <- as.numeric(substring(renamex, first=4))
  
  for (f in 1:dim(test)[2]){
    test[which(test[,f]!=0),f] <- renamex[which(test[,f]!=0)]
  }
  
  sorted <- sort(apply(df[,-1], MARGIN=2, FUN=top), decreasing=T)
  
  matched <- match(names(sorted), names(df[,-1]))
  test <- sapply(as.list(test[matched,]), FUN=rem0)
  
  if (!is.null(tax)){
    if (tax.grp %in% names(tax)){
      names(tax)[which(names(tax)==tax.grp)] <- "taxgrp"
      subtax <- tax[which(tax$taxgrp!=""),]
      test3 <- list()
      for (i in 1:length(test)){
        test3[[names(test)[i]]] <- as.character(unique(subtax$taxgrp[test[[i]]]))
      }
      sink(paste("./coreplots/venn/",file, "euler_group_tax.txt", sep=""))
      print(test3)
      sink()
    }
    
    if (length(test3)<=4){
      venn.diagram(test3, paste("./coreplots/venn/",file, "euler_group_tax.tiff", sep=""),
                   category.names = names(test3), margin=0.12,
                   col=palette()[1:length(test3)], lwd=4,
                   fill=palette()[1:length(test3)], alpha=.125,
                   cat.dist=c(0.1, 0.1, 0.1, 0.1)[1:length(test3)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1.5,
                   fontface="bold", fontfamily="sans")
    } 
    if (length(test3)==5){
      venn.diagram(test3, paste("./coreplots/venn/",file, "euler_group_tax.tiff", sep=""),
                   category.names = names(test3), margin=0.12,
                   col=palette()[1:length(test3)], lwd=4,
                   fill=palette()[1:length(test3)], alpha=.125,
                   cat.dist=c(0.22, 0.22, 0.22, 0.22, 0.22)[1:length(test3)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1,
                   fontface="bold", fontfamily="sans")
    }
    if (length(test3)>5){
      return(test3)
    }
    
    return(test3)
  } else {
    
    sink(paste("./coreplots/venn/",file, "euler_group_relative.txt", sep=""))
    print(test)
    sink()
    
    if (length(test) < 5){
      venn.diagram(test, paste("./coreplots/venn/",file, "euler_group_relative.tiff", sep=""),
                   category.names = names(test), margin=0.12,
                   col=palette()[1:length(test)], lwd=4,
                   fill=palette()[1:length(test)], alpha=.125,
                   cat.dist=c(0.1, 0.1, 0.1, 0.1)[1:length(test)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1.5,
                   fontface="bold", fontfamily="sans")
    } else {
      venn.diagram(test, paste("./coreplots/venn/",file, "euler_group_genus.tiff", sep=""),
                   category.names = names(test), margin=0.12,
                   col=palette()[1:length(test)], lwd=4,
                   fill=palette()[1:length(test)], alpha=.125,
                   cat.dist=c(0.22, 0.22, 0.22, 0.22, 0.22)[1:length(test)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1,
                   fontface="bold", fontfamily="sans")
    }
    return(test)
  }
}






shared.taxa <- function(mbiom, file){
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/venn/")==FALSE){
    dir.create("./coreplots/venn/")
  } 
  
  if (!"xlsx"%in%row.names(installed.packages())) {
    install.packages("xlsx")
  }
  require(xlsx)
  
  a <- combn(names(mbiom), 2)
  
  if (length(mbiom)>3){
    b <- combn(names(mbiom), 3)
  }
  if (length(mbiom)>4){
    c <- combn(names(mbiom), 4)
  }
  
  out <- list()
  for (i in 1:length(mbiom)){
    x <- mbiom[[i]][which(!mbiom[[i]]%in%unlist(mbiom[-i], use.names=F))]
    
    out[[names(mbiom)[i]]] <- x
  }
  
  for (i in 1:dim(a)[2]){
    x <- mbiom[[a[1,i]]][which(mbiom[[a[1,i]]]%in%mbiom[[a[2,i]]] & 
                                 !mbiom[[a[1,i]]]%in%unlist(mbiom[which(!names(mbiom)%in%a[,i])]))]
    
    out[[paste(a[1,i],a[2,i], sep="-")]] <- x
  }
  
  if (length(mbiom)==3){
    all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]] & 
                              mbiom[[1]]%in%mbiom[[3]])]
  }
  
  if (length(mbiom)>3){
    for (i in 1:dim(b)[2]){
      x <- mbiom[[b[1,i]]][which(mbiom[[b[1,i]]]%in%mbiom[[b[2,i]]] & 
                                   mbiom[[b[1,i]]]%in%mbiom[[b[3,i]]] &
                                   !mbiom[[b[1,i]]]%in%unlist(mbiom[which(!names(mbiom)%in%b[,i])]))]
      
      out[[paste(b[1,i],b[2,i],b[3,i], sep="-")]] <- x
    }
    if (length(mbiom)==4){
      all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]] & 
                                mbiom[[1]]%in%mbiom[[3]] &
                                mbiom[[1]]%in%mbiom[[4]])]
    }
  }
  
  if (length(mbiom)==5){
    for (i in 1:dim(c)[2]){
      x <- mbiom[[c[1,i]]][which(mbiom[[c[1,i]]]%in%mbiom[[c[2,i]]] & 
                                   mbiom[[c[1,i]]]%in%mbiom[[c[3,i]]] &
                                   mbiom[[c[1,i]]]%in%mbiom[[c[4,i]]] &
                                   !mbiom[[c[1,i]]]%in%unlist(mbiom[which(!names(mbiom)%in%c[,i])]))]
      
      out[[paste(c[1,i],c[2,i],c[3,i],c[4,i], sep="-")]] <- x
    }
    
    all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]] & 
                              mbiom[[1]]%in%mbiom[[3]] &
                              mbiom[[1]]%in%mbiom[[4]] &
                              mbiom[[1]]%in%mbiom[[5]])]
  }
  
  
  out[["all"]] <- all
  
  if (file.exists(paste0("./coreplots/venn/",file, ".xlsx"))){
    file.remove(paste0("./coreplots/venn/",file, ".xlsx"))
  }
  for (i in 1:length(out)){
    write.xlsx(out[i], file=paste0("./coreplots/venn/",file, ".xlsx"), sheetName = names(out)[i], append=T)
  }
  
  return(out)
}






mbiom.bar <- function(shared, comm, select=NULL, tax, hi.tax, meta, group, subgroup, file="mbiom_bar", legend=F){
  
  cran <- c("dplyr", "ggplot2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  require(dplyr)
  require(ggplot2)
  
  colors <- c("#330033", "#F0E442", "#0072B2", "#FF0033", "#FF9933", "#999999", "#56B4E9", "#FFCC00", "#00FF00", "#FF66FF", "#666666")
  
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
  }
  
  if (any(!row.names(comm)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/venn/")==FALSE){
    dir.create("./coreplots/venn/")
  } 
  if(dir.exists(paste0("./coreplots/venn/", file))==FALSE){
    dir.create(paste0("./coreplots/venn/", file))
  } 
  
  splt <- function(x){
    y <- strsplit(x, split="-")
    y <- unlist(y)
  }
  
  ra <- function(x){ #where x is a vector
    y <- sum(x)
    z <- sapply(x, FUN=function(x){z1 <- x/y})
    return(z)
  }
  
  sel.list <- as.list(names(shared))
  sel.list <- lapply(sel.list, FUN=splt)
  
  if (is.null(select)){select <- names(comm)} 
  
  s <- which(names(meta)=="sample")
  g <- which(names(meta)==group)
  sg <- which(names(meta)==subgroup)
  
  meta <- mutate(meta, !!"group.subgroup" := 
                   paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep="."))
  g.sg <- which(names(meta)=="group.subgroup")
  
  ht <- which(names(tax)==hi.tax)
  
  for (j in 1:dim(tax)[2]){
    t <- all(!is.na(match(shared$all, tax[,j])))
    if (t==TRUE){break}
  }
  
  for (i in 1:length(sel.list)){
    
    selection <- sel.list[[i]]
    
    if (selection=="all"){selection <- unique(unlist(sel.list))}
    
    otus <- tax$tag[match(shared[[i]], tax[,j])]
    otus <- select[select%in%otus]
    
    if (length(otus)==0){next}
    
    
    data <- data.frame(sample=row.names(comm), comm)
    data <- right_join(meta[,c(1,g.sg)], data)
    
    data2 <- data[,-1] %>%
      group_by(group.subgroup) %>%
      summarize_all(funs(mean=mean)) %>%
      as.data.frame()
    
    new <- apply(data2[,-1], MARGIN=1, FUN=ra) %>% t()
    
    new <- data.frame(groupvar=data2[,1], new)
    names(new)[-1] <- names(comm)
    
    data2 <- reshape2::melt(new)
    
    
    data2 <- right_join(tax[,c(3, ht)], data2, by=c("tag"="variable"))
    data2 <- right_join(unique(data.frame(meta[,c(g.sg, g, sg)])), data2, by=c("group.subgroup"="groupvar"))
    
    ce.il <- data2[data2[,group]%in%selection & data2$tag%in%otus,]
    ce.il$value <- ce.il$value*100
    
    names(ce.il)[which(names(ce.il)==hi.tax)] <- "Taxonomy"
    names(ce.il)[which(names(ce.il)==group)] <- "Group"
    names(ce.il)[which(names(ce.il)==subgroup)] <- "Subgroup"
    
    colors <- c("#330033", "#F0E442", "#0072B2", "#FF0033", "#FF9933", "#999999", "#56B4E9", "#FFCC00", "#00FF00", "#FF66FF", "#666666")
    
    if (legend==T){
      leg <- expression(guides(fill=guide_legend(title="", ncol=1)))
    } else {leg <- expression(guides(fill=F))}
    
    f <- ggplot(ce.il, aes(Subgroup, value, fill=Taxonomy)) +
      facet_wrap(vars(Group), dir="v", scales="free_y", strip.position = "right") +
      geom_bar(stat="identity", width=.96) +
      theme_minimal() +
      scale_fill_manual(values=colorRampPalette(colors)(length(unique(as.factor(ce.il$Taxonomy))))) +
      theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
            legend.position="right",
            legend.justification="right",
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = 80, face="bold"),
            axis.title = element_text(size = 65, face="bold"),
            axis.text.x = element_text(color = "black", size = 60, angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = 55),
            legend.title = element_text(size = 50, face="bold"),
            legend.text = element_text(size = 50),
            panel.spacing.y = unit(2, "lines"),
            plot.margin = margin(3, 3, 3, 3, "lines")) +
      ylab(NULL)+
      xlab(NULL)+
      eval(expr=leg)
    
    w <- 
      
      if (length(selection)==4){
        jpeg(paste0("./coreplots/venn/", file, "/", selection[1],"-", selection[2], 
                    "-", selection[3],"-", selection[4],"_sharedOTUs.jpg"),
             width=2300, height=1250)
      }
    if (length(selection)==3){
      jpeg(paste0("./coreplots/venn/", file, "/",selection[1],"-", selection[2], 
                  "-", selection[3],"_sharedOTUs.jpg"),
           width=1300, height=1650)
    }
    if (length(selection)==2){
      jpeg(paste0("./coreplots/venn/", file, "/", selection[1],"-", selection[2], "_sharedOTUs.jpg"),
           width=1300, height=1150)
    }
    if (length(selection)==1){
      jpeg(paste0("./coreplots/venn/", file, "/", selection[1],"_sharedOTUs.jpg"),
           width=1300, height=700)
    }
    print(f)
    dev.off()
  }
  
}