## The following functions were written by Kara JM Taylor for use with 
## 16S rRNA metabarcoding data, as output in table format by QIIME2-2019.1.
## It was written to produce the data presented in Taylor et al. (2020). 
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
    #source("https://bioconductor.org/biocLite.R")
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")}
      #biocLite(i, suppressUpdates=T)
      BiocManager::install("biomformat")
  }
  library(dplyr)
  library(biomformat)
  
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
  write.csv(tax, "taxonomy.csv", row.names=F)
}





bact.tax <- function(taxonomy, database=NULL){
  if (!"dplyr"%in%row.names(installed.packages())) {
    install.packages("dplyr")
  }
  library(dplyr)
  
  splt <- function(x){
    y <- strsplit(x, split=";", fixed=TRUE)
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
    library(stringr)
    
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
  #return(evened)
  
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
  
  out$phylum <- as.character(out$phylum)
  out$phylum[out$phylum=="Deinococcus"] <- "Deinococcus-Thermus"
  out$phylum <- as.factor(out$phylum)
  
  #if genus-species not available, use smallest level possible (but not strain codes; 
  #use regular expressions to pull out any names given as alphanumeric codes)
  #and then follow that name with either "sp" or the strain code
  tax <- out
  save(tax, file="./taxonomy.RD")
  
  return(out)
}





pcoa.plotting <- function(dist, meta, group, colors="rainbow", method="", axes=c(1,2,3), fixed=NULL){
  
  
  pkgs <- c("dplyr", "ape", "ggplot2", "gridExtra", "RColorBrewer", "pairwiseAdonis", "ggpubr")
  if (any(!pkgs%in%row.names(installed.packages()))) {
    if (!"pairwiseAdonis"%in%row.names(installed.packages())){
      devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
    }
    install.packages(pkgs[which(!pkgs%in%row.names(installed.packages())==T)])
  }
  
  library(dplyr)
  library(ape)
  library(ggplot2)
  library(gridExtra)
  library(ggpubr)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/pcoa/")==FALSE){
    dir.create("./coreplots/pcoa/")
  }
  if(dir.exists(paste0("./coreplots/pcoa/", deparse(substitute(dist))))==FALSE){
    dir.create(paste0("./coreplots/pcoa/", deparse(substitute(dist))))
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
    new_names[i] <- paste("PC",i, sep="")
  }
  names(twPCOA) <- new_names
  
  twPCOA <- left_join(data.frame(twPCOA, sample=row.names(twPCOA)), 
                      meta[,c(id,g)])
  
  g2 <- which(names(twPCOA)==group)
  
  if (colors=="rainbow"){
    palette(c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet"))
    if (length(unique(as.factor(meta[,g]))) > 9){
      palette(colorRampPalette(palette())(length(unique(as.factor(twPCOA[,g2])))))
    }
  }
  if (colors=="default"){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    palette(gg_color_hue(length(unique(as.factor(meta[,g])))))
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
  if (colors=="paired"){
    gg_color_hue2 <- function(n) {
      hues = seq(240, 580, length = n+1)
      hcl(h = hues, l = 50, c = 120)[1:n]}
    gg_color_hue3 <- function(n) {
      hues = seq(240, 580, length = n+1)
      hcl(h = hues, l = 90, c = 120)[1:n]}
    palette(c(rbind(gg_color_hue2(length(unique(twPCOA[,dim(twPCOA)[2]]))/2),
                    gg_color_hue3(length(unique(twPCOA[,dim(twPCOA)[2]]))/2))))
    palette(colorRampPalette(palette())(length(unique(twPCOA[,dim(twPCOA)[2]]))))
  }
  #return(palette())
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
  
  pleg <- as_ggplot(get_legend(p12))
  
  
  
  tiff(paste("./coreplots/pcoa/", deparse(substitute(dist)), "/PCOA_", as.character(axes)[1], as.character(axes)[2], 
             as.character(axes)[3], paste(deparse(substitute(dist)), group, sep="_"), 
             ".tiff", sep=""),
       width=9000, height=4000, res=600)
  pp <- arrangeGrob(p12+guides(fill=F), p13+guides(fill=F), p23+guides(fill=F), pleg, 
                    layout_matrix=rbind(c(1,2,3), c(1,2,3), c(4,4,4)))
  #grid.arrange(p12, p13, p23, newpage=FALSE, nrow=1, ncol=3)
  grid.arrange(pp)
  dev.off()
  
  
  sink(paste("./coreplots/pcoa/", deparse(substitute(dist)), "/PCOA_", as.character(axes)[1], as.character(axes)[2], 
             as.character(axes)[3], paste(deparse(substitute(dist)), group, sep="_"), 
             ".txt", sep=""))
  print(vegan::adonis(dist ~ twPCOA[,g2]))
  print(pairwiseAdonis::pairwise.adonis(dist, twPCOA[,g2]))
  sink()
  
}






core.id <- function(comm, meta, tax, group, subgroup, margin){ 
  # where x is a sample-by-OTU table only, rows as samples, columns as otus
  # Row.names should be the sample IDs.
  # Any metadata should be input as 'factors', 
  # and make sure that the ID column in the metadata is labeled "SampleID"
  if (!"dplyr"%in%row.names(installed.packages())) {install.packages("dplyr")}
  library(dplyr)
  
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
                       date, core.only=FALSE, fixed=FALSE, landscape=TRUE, view.all=FALSE, facet=FALSE){
  
  cran <- c("dplyr", "RColorBrewer", "ggplot2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/profiles/")==FALSE){
    dir.create("./coreplots/profiles/")
  } 
  if(dir.exists(paste0("./coreplots/profiles/", group, "_", date, "/"))==FALSE){
    dir.create(paste0("./coreplots/profiles/", group, "_", date, "/"))
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
          
          otus <- names(which(summary(pool$otu)== ##Note: this line means that the function ignores if some age groups don't meet the occurence margin parameter
                                length(which(is.na(list[[i]][j:length(list[[i]])])==FALSE))))
          #cat(paste0(names(list)[i], ": ", names(list[[i]][j]), "-->", names(list[[i]][length(list[[i]])])), otus, sep="\n")
          if ("(Other)"%in%otus){otus <- otus[-length(otus)]}
          if (length(otus)>0 & length(which(is.na(list[[i]])==FALSE))>=1){
            sub <- data[which(data[,g]==names(list[i])),
                        which(names(data)%in%c(subgroup, otus))]
            
            if (is.null(tax)==FALSE){
              otus <- otus[match(names(sub)[-1], otus)]
              otus <- data.frame(v1=c(1:length(otus)), v2=otus)
              otus <- left_join(otus, tax, by=c("v2"="tag"))
              
              
              if (hi.tax%in%names(tax)){
                get <- which(names(otus)==hi.tax)
                otus <- as.character(unlist(otus[,get]))
              } else {otus <- as.character(otus$otu.name)}
              
              #if (length(otus)!=dim(t(sub[,-1]))[1]){return(list(features=otus, comm=t(sub[,-1])))}
              #if (length(otus)!=dim(t(sub[,-1]))[1]){return(summary(pool$otu))}
              
              #otus <- otus[!is.na(otus)] ##This prevents the summary(pool$otu) from trying to return "(Other)" as an otu by accident
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
            #if (i==5 & j==4){return(pool)}
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
    
  }
  
  colors <- c("#330033", "#F0E442", "#0072B2", "#FF0033", "#FF9933", "#999999", "#56B4E9", "#FFCC00", "#00FF00", "#FF66FF", "#666666")
  colors <- colorRampPalette(colors)(length(unique(df$otu)))
  names(colors) <- sort(as.character(unique(df$otu)))
  
  
  for (i in unique(df$groupvar)){
    
    group.core <- subset(df, groupvar==i)
    
    group.core[,subgroup] <- factor(group.core[,subgroup], levels=unique(factors[,subgroup]))
    
    
    face <- c(floor(length(unique(group.core$core))/3), ceiling(length(unique(group.core$core))/3))
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
    
    fexp <- expression(facet_grid(core~., drop=TRUE))
    fexp2 <- expression(facet_grid(groupvar~core, drop=FALSE))
    
    leg <- list("right", "right")
    
    if (core.only==FALSE){
      if (length(unique(group.core$otu)) > 40){
        kh <- floor(length(unique(group.core$otu))/40)
      } 
      if (length(unique(group.core$otu))/length(unique(group.core$core)) > 8){
        kh <- length(unique(group.core$otu))/8
      }
    }
    cols <- 1
    
    
    
    
    if (fixed==TRUE){
      fix.ylim <- expression(ylim(0,100))
    } 
    if (fixed==FALSE){
      fix.ylim <- expression(ylim(0, NA))
    }
    if (fixed=="free"){
      fix.ylim <- expression(facet_wrap(vars(core), dir="v", scales="free_y", strip.position = "right", drop=T))
    }
    
    
    
    group.core$rel.abund <- group.core$rel.abund*100
    
    #lev <- match(levels(with(meta, eval(parse(text=subgroup)))),
    #             unique(with(group.core, eval(parse(text=subgroup)))))
    #lev <- lev[which(!is.na(lev))]
    
    #group.core[,which(names(group.core)==subgroup)] <-
    #  factor(group.core[,which(names(group.core)==subgroup)],
    #         levels=levels(with(group.core, eval(parse(text=subgroup))))[
    #           lev[which(is.na(lev)==FALSE)]])
    
    w <- 500*length(unique(group.core[,1]))+1000
    r <- expression(300)
    fs <- w/300
    fs <- round(c(fs+(fs/2), fs+(fs/3.2), fs+(fs/4), fs+(fs/8), fs+(fs/16), fs*2))
    
    if (view.all==T){
      fexp <- expression(facet_grid(groupvar~., drop=TRUE))
      fexp2 <- expression(facet_grid(groupvar~., drop=TRUE))
      h <- 2300
    } else {
      h <- 800*length(unique(group.core$core))+1200
    }
    
    group.core$otu <- factor(group.core$otu, levels=unique(group.core$otu)) #makes sure the otus are in the right order
    group.core$core <- factor(group.core$core, levels=unique(group.core$core))
    
    f <- ggplot(group.core, aes(eval(parse(text=subgroup)), rel.abund, fill=otu)) +
      eval(expr=fexp) +
      geom_bar(stat="identity") +
      eval(expr=fix.ylim) +
      theme_minimal() +
      scale_fill_manual(values=colors) +
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
            axis.text.x = element_text(color = "black", size = fs[2], angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = fs[2], face="bold"),
            legend.title = element_text(size = fs[3], face="bold"),
            legend.text = element_text(size = fs[3]),
            plot.caption = element_text(size=fs[5]),
            legend.background=element_rect(fill="white", colour="black"))
    
    jpeg(file=paste("./coreplots/profiles/", group, "_", date, "/", i, date, "_", threshold, "_", as.character(margin), ".jpeg", sep=""),
         width=w, height=h, units="px", res=eval(expr=r) )
    print(f + labs(title=paste("Progression of Core Microbiome Composition in", i)) + 
            xlab(subgroup) + ylab("Relative Abundance (%)") +
            guides(fill=guide_legend(keywidth = fs[2]/10, keyheight=fs[2]/(10*kh), ncol=cols)) +
            labs(caption=paste("Occurence margin: ", as.character(margin*100), "%", "\n",
                               "Abundance threshold: ", as.character(threshold*100), "%")))
    dev.off()
  }
  write.csv(df, paste("./coreplots/profiles/", group, "_", date, "/", date, "_",
                      threshold, "_", as.character(margin), ".csv", sep=""))
  
  if (facet==TRUE){
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
    
    #df[,which(names(df)==subgroup)] <-
    #  factor(df[,which(names(df)==subgroup)],
    #         levels=levels(with(df, eval(parse(text=subgroup))))[
    #           lev[which(is.na(lev)==FALSE)]])
    
    bff <- ggplot(df, aes(eval(parse(text=subgroup)), rel.abund, fill=otu)) +
      geom_bar(stat="identity") +
      eval(expr=fix.ylim) +
      eval(expr=fexp2) +
      theme_minimal() +
      scale_fill_manual(values=colors) +
      theme(panel.border=element_rect(linetype="solid", fill=NA),
            legend.position=leg[[1]],
            legend.justification=leg[[2]],
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = fs[6], face="bold"),
            plot.title = element_text(size = fs[1], face="bold"),
            axis.title = element_text(size = fs[2], face="bold"),
            axis.text.x = element_text(color = "black", size = fs[6]),#, angle=330, hjust=0, vjust=0.5),
            axis.text.y = element_text(color = "black", size = fs[6], face="bold"),
            legend.title = element_text(size = fs[3], face="bold"),
            legend.text = element_text(size = fs[1]),
            plot.caption = element_text(size=fs[5]),
            legend.background=element_rect(fill="white", colour="black")) +
      eval(expr=legs)
    jpeg(file=paste("./coreplots/profiles/", group, "_", date, "/facet", date, "_", threshold, "_", as.character(margin), ".jpeg", sep=""),
         width=w, res=600, units="px", height=h) #height=h*length(list))
    print(bff + labs(title=paste("Progression of Core Microbiome Composition")) + 
            xlab(subgroup) + ylab("Relative Abundance (%)") +
            eval(expr=gui) +
            labs(caption=paste("Occurence margin: ", as.character(margin*100), "%", "\n",
                               "Abundance threshold: ", as.character(threshold*100), "%")))
    dev.off()
  }
  return(df)
}






dendro.heatmap <- function(comm, tax, meta, path="", group, subgroup, core.list=NULL, 
                           hi.tax="tag", file="", classification=NULL, trans="log", stop="no", wd=NULL, ht=NULL){
  #### Series of nested functions to produce grouped heatmap
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
    id <- 1
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
  if(dir.exists(paste0("./coreplots/profiles/", group, "_", file, "/"))==FALSE){
    dir.create(paste0("./coreplots/profiles/", group, "_", file, "/"))
  } 
  
  cran <- c("dplyr", "egg", "ape", "vegan", "ggplot2", "reshape2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  biocon <- c("ggtree", "ggdendro", "dendextend")
  if (any(!biocon%in%row.names(installed.packages()))){
    #source("https://bioconductor.org/biocLite.R")
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")}
    for (i in biocon){
      #biocLite(i, suppressUpdates=T)
      BiocManager::install(i)
    }
  }
  
  library(dplyr)
  library(egg)
  library(ape)
  library(vegan)
  library(ggplot2)
  library(ggtree)
  library(ggdendro)
  library(dendextend)
  
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
    
    data <- left_join(new, meta[,c(id, g, sg)], by=c("ID"="sample"))
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
    
    if (stop=="data"){return(group3)}
    
    return(h)
  }
  
  top <- sub.tax(comm, tax, meta, group, subgroup, classification, core.list, rank)
  if (stop=="top"){return(top)}
  den <- group.dendro(comm, meta, group, subgroup)
  if (stop=="den"){return(den)}
  phy <- otu.phylo(top, tax, path, hi.tax) #use hi.tax if you want core OR all taxa at level other than in `sub`
  if (stop=="phy"){return(phy)}
  heat <- top.heat(top, tax, phy, den, hi.tax, file) #use hi.tax if you want core OR all taxa
  if (stop=="heat" | stop=="data"){return(heat)}
  
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
  n2 <- length(den$tip.labels)
  
  if (is.null(wd)){
    wd <- 1000+65*n2
  }
  if (is.null(ht)){
    if (n>=50){
      ht <- 56*n
    } 
    if (n<50 & n>=25){
      ht <- 800+65*n
    }
    if (n<25){
      ht <- 2000
    }
  }
  
  hts <- c(0.4-0.003*n, 3)
  wds <- c(0.5-0.003*n2, 3)
  
  
  tiff(paste0("./coreplots/profiles/", group, "_", file, "/phyloheatmap_", file, ".tiff"), width=wd, height=ht)
  egg::ggarrange(e, den, phy, heat, heights=hts, widths=wds)
  dev.off()
  
}





mbiom.venn <- function(mat, meta, group, tax=NULL, file="./", tax.grp="tag") {
  
  cran <- c("dplyr", "VennDiagram", "RColorBrewer")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  library(dplyr)
  library(VennDiagram)
  library(RColorBrewer)
  
  palette(c(RColorBrewer::brewer.pal(8, "Set1"), "#000000"))
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
    id <- 1
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
  if(dir.exists(paste0("./coreplots/venn/", file))==FALSE){
    dir.create(paste0("./coreplots/venn/", file))
  } 
  
  g <- which(names(meta)==group)
  
  mat2 <- data.frame(sample=row.names(mat), mat)
  mat <- right_join(meta[,c(id,g)], mat2)
  
  
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
  #return(test)
  sorted <- sort(apply(df[,-1], MARGIN=2, FUN=top), decreasing=T)
  #return(sorted)
  matched <- match(names(sorted), names(df[,-1]))
  #return(matched)
  #return(as.list(test[matched,]))
  test <- sapply(as.list(test[matched,]), FUN=rem0)
  #return(test)
  if (!is.null(tax)){
    if (tax.grp %in% names(tax)){
      names(tax)[which(names(tax)==tax.grp)] <- "taxgrp"
      subtax <- tax[which(tax$taxgrp!=""),]
      test3 <- list()
      for (i in 1:length(test)){
        test3[[names(test)[i]]] <- as.character(unique(subtax$taxgrp[test[[i]]]))
      }
      sink(paste("./coreplots/venn/",file, "/euler_group_tax.txt", sep=""))
      print(test3)
      sink()
    }
    #return(test3)
    if (length(test3)<=4){
      venn.diagram(test3, paste("./coreplots/venn/",file, "/euler_group_tax.tiff", sep=""),
                   category.names = names(test3), margin=0.12,
                   col=palette()[1:length(test3)], lwd=4,
                   fill=palette()[1:length(test3)], alpha=.125,
                   cat.dist=c(0.1, 0.1, 0.1, 0.1)[1:length(test3)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1.5,
                   fontface="bold", fontfamily="sans")
    } 
    if (length(test3)==5){
      venn.diagram(test3, paste("./coreplots/venn/",file, "/euler_group_tax.tiff", sep=""),
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
    
    sink(paste("./coreplots/venn/",file, "/euler_group_relative.txt", sep=""))
    print(test)
    sink()
    
    if (length(test) < 5){
      venn.diagram(test, paste("./coreplots/venn/",file, "/euler_group_relative.tiff", sep=""),
                   category.names = names(test), margin=0.12,
                   col=palette()[1:length(test)], lwd=4,
                   fill=palette()[1:length(test)], alpha=.125,
                   cat.dist=c(0.1, 0.1, 0.1, 0.1)[1:length(test)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1.5,
                   fontface="bold", fontfamily="sans")
    } else {
      venn.diagram(test, paste("./coreplots/venn/",file, "/euler_group_genus.tiff", sep=""),
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






shared.taxa <- function(mbiom, file, excel=T){
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/venn/")==FALSE){
    dir.create("./coreplots/venn/")
  } 
  if(dir.exists(paste0("./coreplots/venn/", file))==FALSE){
    dir.create(paste0("./coreplots/venn/", file))
  }
  
  if (!"xlsx"%in%row.names(installed.packages())) {
    install.packages("xlsx")
  }
  library(xlsx)
  
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
    
    if (length(mbiom)==2){
      all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]])]
    }
  }
  
  for (i in 1:dim(a)[2]){
    if (length(mbiom)!=2){
      x <- mbiom[[a[1,i]]][which(mbiom[[a[1,i]]]%in%mbiom[[a[2,i]]] & 
                                   !mbiom[[a[1,i]]]%in%unlist(mbiom[which(!names(mbiom)%in%a[,i])]))]
      
      out[[paste(a[1,i],a[2,i], sep=".")]] <- x
    }
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
      
      out[[paste(b[1,i],b[2,i],b[3,i], sep=".")]] <- x
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
      
      out[[paste(c[1,i],c[2,i],c[3,i],c[4,i], sep=".")]] <- x
    }
    
    all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]] & 
                              mbiom[[1]]%in%mbiom[[3]] &
                              mbiom[[1]]%in%mbiom[[4]] &
                              mbiom[[1]]%in%mbiom[[5]])]
  }
  
  
  out[["all"]] <- all
  
  if (excel==T){
    if (file.exists(paste0("./coreplots/venn/",file, "/", file, ".xlsx"))){
      file.remove(paste0("./coreplots/venn/",file, "/", file, ".xlsx"))
    }
    for (i in 1:length(out)){
      if (length(out[[i]])>0){
        write.xlsx(out[i], file=paste0("./coreplots/venn/",file, "/", file, ".xlsx"), sheetName = names(out)[i], append=T)
      }
    }
  }
  return(out)
}






mbiom.bar <- function(shared, comm, select=NULL, tax, hi.tax, meta, group, subgroup, file="mbiom_bar", legend=F, compare.otus=F, set.y=F){
  
  cran <- c("dplyr", "ggplot2", "ggpubr")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
    id <- 1
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
    y <- strsplit(x, split="\\.")
    y <- unlist(y)
  }
  
  ra <- function(x){ #where x is a vector
    y <- sum(x)
    z <- sapply(x, FUN=function(x){z1 <- x/y})
    return(z)
  }
  
  sel.list <- as.list(names(shared))
  
  if (length(sel.list)>4){
    sel.list <- lapply(sel.list, FUN=splt)
  }
  
  if (is.null(select)){select <- names(comm)} 
  
  s <- which(names(meta)=="sample")
  g <- which(names(meta)==group)
  sg <- which(names(meta)==subgroup)
  
  meta <- mutate(meta, !!"group.subgroup" := 
                   paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep=":"))
  g.sg <- which(names(meta)=="group.subgroup")
  
  ht <- which(names(tax)==hi.tax)
  
  for (j in 1:dim(tax)[2]){
    t <- all(!is.na(match(unlist(shared), tax[,j])))
    if (t==TRUE){break}
  }
  
  df <- data.frame()
  for (i in 1:length(sel.list)){
    
    selection <- sel.list[[i]]
    
    if ("all"%in%selection){
      selection <- unique(unlist(sel.list))
      selection <- selection[-length(selection)]
    }
    
    otus <- tax$tag[match(shared[[i]], tax[,j])]
    otus <- select[select%in%otus]
    
    if (length(otus)==0){next}
    
    
    data <- data.frame(sample=row.names(comm), comm)
    data <- right_join(meta[,c(id,g.sg)], data)

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
    
    ce.il <- data.frame(ce.il, sec.list=paste(selection, collapse="."))
    
    df <- rbind(df, ce.il)
    
  }
  
  write.csv(df, paste0("./coreplots/venn/", file, "/dataframe.csv"))
  

  if (compare.otus==T){
    g <- ggplot(df, aes(Subgroup, value, fill=Group)) +
      facet_wrap(vars(Taxonomy), ncol=4, dir="h", scales="free_y", strip.position = "top") +
      geom_bar(stat="summary", fun.y="sum", width=.96, position="dodge") +
      theme_minimal() +
      theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
            legend.position="bottom",
            legend.justification="center",
            legend.direction="horizontal",
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = 45, face="bold"),
            axis.title.y = element_text(size = 40, face="bold"),
            axis.text.x = element_text(color = "black", size = 40, angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = 40),
            legend.title = element_text(size = 40, face="bold"),
            legend.text = element_text(size = 30),
            panel.spacing.y = unit(2, "lines"),
            plot.margin = margin(3, 3, 3, 3, "lines")) +
      ylab("Relative Abundance (%)")+
      xlab(NULL)
    h <- ceiling(length(unique(df$Taxonomy))/4)
    jpeg(paste0("./coreplots/venn/", file, "/", group,"_sharedOTUs.jpg"),
         width=2300, height=420*h+280)
    print(g)
    dev.off()
  }
  
  #return(ce.il)
  colors <- c("#330033", "#F0E442", "#0072B2", "#FF0033", "#FF9933", "#999999", "#56B4E9", "#FFCC00", "#00FF00", "#FF66FF", "#666666")
  colors <- colorRampPalette(colors)(length(unique(as.factor(df$Taxonomy))))
  names(colors) <- as.character(unique(df$Taxonomy))
  
  if (compare.otus==T){
    k <- ggplot(df, aes(Subgroup, value, fill=Taxonomy)) +
      facet_wrap(vars(Group), ncol=2, dir="v", scales="free_y", strip.position = "right") +
      geom_bar(stat="summary", fun.y="sum", width=.96, position="stack") +
      theme_minimal() +
      scale_fill_manual(values=colors) +
      theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
            legend.position="right",
            legend.justification="center",
            #legend.direction="horizontal",
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = 60, face="bold"),
            axis.title.y = element_text(size = 40, face="bold"),
            axis.text.x = element_text(color = "black", size = 40, angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = 40),
            legend.title = element_text(size = 40, face="bold"),
            legend.text = element_text(size = 30),
            panel.spacing.y = unit(2, "lines"),
            plot.margin = margin(3, 3, 3, 3, "lines"))+
      ylab("Relative Abundance (%)")+
      xlab(NULL)
    h <- ceiling(length(unique(df$Group))/2)
    jpeg(paste0("./coreplots/venn/", file, "/", group,"_NicheSpace.jpg"),
         width=2300, height=460*h)
    print(k + guides(fill=F))
    dev.off()
    
    pleg <- as_ggplot(get_legend(k))
    if (length(unique(as.factor(df$Taxonomy))) > 13){
      wid <- length(unique(as.factor(df$Taxonomy)))/14*900
    } else {wid <- 900}
    jpeg(paste0("./coreplots/venn/", file, "/", group,"_NicheSpace_legend.jpg"),
         width=wid, height=1150)
    print(pleg)
    dev.off()
  }
  
  if (set.y==T){
    fix.ylim <- expression(ylim(0,100))
  } else {fix.ylim <- expression(ylim(0,NA))}
  
  for (i in unique(df$sec.list)){
    
    sub.sel <- subset(df, subset=sec.list==i)
    
    
    f <- ggplot(sub.sel, aes(Subgroup, value, fill=Taxonomy)) +
      facet_wrap(vars(Group), dir="v", scales="free_y", strip.position = "right") +
      geom_bar(stat="identity", width=.96) +
      theme_minimal() +
      scale_fill_manual(values=colors) +
      eval(expr=fix.ylim) +
      theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
            legend.position="right",
            legend.justification="center",
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
      xlab(NULL)
    #eval(expr=leg)
    w <- length(unique(df$Subgroup))
    
    if (legend==T){
      #leg <- expression(guides(fill=guide_legend(title="", ncol=1)))
      pleg <- as_ggplot(get_legend(f))
      if (length(unique(as.factor(sub.sel$Taxonomy))) > 13){
        wid <- length(unique(as.factor(sub.sel$Taxonomy)))/14*900
      } else {wid <- 900}
    } #else {leg <- expression(guides(fill=F))}
    
    count <- unlist(strsplit(i, split="\\."))
    
    if (length(count)==5){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()
      f <- f + facet_wrap(vars(Group), dir="v", scales="fixed", strip.position = "right")
      
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=2300, height=1650)
    }
    if (length(count)==4){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()
      f <- f + facet_wrap(vars(Group), dir="v", scales="fixed", strip.position = "right")
      
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=2300, height=1250)
    }
    if (length(count)==3){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()
      
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=1300, height=1650)
    }
    if (length(count)==2){
      jpeg(paste0("./coreplots/venn/", file, "/", i, "_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()
      
      jpeg(paste0("./coreplots/venn/", file, "/", i, "_sharedOTUs.jpg"),
           width=1300, height=1150)
    }
    if (length(count)==1){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()
      
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=1300, height=700)
    }
    print(f+guides(fill=F))
    dev.off()
    
    
  }
  graphics.off()
  return(df)
}


alpha.boxplots <- function(comm, meta, select=NULL, group, legend=TRUE, fixed=NULL, type="box"){
  
  pkgs <- c("vegan", "dplyr","fossil")
  if (any(!pkgs%in%row.names(installed.packages()))) {
    install.packages(pkgs[which(!pkgs%in%row.names(installed.packages())==T)])
  }
  library(dplyr)
  library(vegan)
  library(fossil)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/alpha/")==FALSE){
    dir.create("./coreplots/alpha/")
  }
  if(dir.exists(paste0("./coreplots/alpha/", deparse(substitute(comm))))==FALSE){
    dir.create(paste0("./coreplots/alpha/", deparse(substitute(comm))))
  }
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    names(meta)[1] <- "sample"
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    id <- 1
  }
  
  
  g <- which(names(meta)%in%group)
  
  
  if (is.null(select)==TRUE){
    sep.comm <- comm
  } else {
    sep.comm <- comm[grep(select, row.names(comm)),]
  }
  
  if (any(!row.names(sep.comm)%in%meta$sample)){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  shan <- diversity(sep.comm)
  rich <- specnumber(sep.comm)
  even <- shan/log(rich)
  #chao <- apply(sep.comm, MARGIN=1, FUN=chao1)
  
  shan <- data.frame(sample=names(shan), shan=shan, rich=rich, even=even)#, chao1=chao)
  
  
  
  shan <- left_join(shan, meta[,c(id,g)])
  #shan[,group] <- factor(shan[,group], levels = unique(shan[,group]))
  
  names(shan) <- c("sample", "shan", "rich", "even", "group.id")
  
  
  shan$group.id <- factor(shan$group.id,
                          levels=sort(unique(shan$group.id)))
  shan <- shan[order(shan$group.id),]
  
  clr <- c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet")
  clr <- colorRampPalette(clr)(length(unique(shan$group.id)))
  names(clr) <- as.character(unique(shan$group.id))
  palette(clr)
  
  or.par <- par( no.readonly = TRUE )
  
  if (type=="box"){
    jpeg(paste("./coreplots/alpha/",deparse(substitute(comm)),"/ALPHA_boxplot_", group, select, ".jpg", sep=""), 6000, 2300, res=800)
    par(mfrow=c(1,3), cex.axis=1.1, font.axis=2, las=3, xaxt="n")
    
    out <- list()
    fixedx <- 0
    for (j in c("shan", "rich", "even")){#, "chao1")){
      if (is.null(fixed)){
        boxplot(eval(parse(text=j)) ~ group.id, data=shan, boxfill=palette(), 
                boxwex=.8, main=select, ylab=j)
      } else {
        fixedx <- fixedx+1
        boxplot(eval(parse(text=j)) ~ group.id, data=shan, boxfill=palette(), 
                boxwex=.8, main=select, ylab=j, ylim=fixed[[fixedx]])
      }
      axis(1, at=c(1:length(levels(shan$group.id))), labels=FALSE)
      text(c(1:length(levels(shan$group.id))), par("usr")[3], font=2, 
           labels=c(levels(shan$group.id)), srt=45, adj=c(1, 1.5), offset=2.5)
      if (legend==TRUE){
        if (j=="shan"){
          legend(c(0,-0.5), horiz = TRUE, fill = palette(), legend = levels(shan$group.id), xpd=T, bty = 'n', cex=0.5)
        }
      }
      
      atw <- aov(eval(parse(text=j)) ~ group.id, data=shan)
      tuktw <- TukeyHSD(atw)
      source("D:/functions/tukey_select.R")
      tuksel <- tukey.select(tuktw$`group.id`, 1)
      
      test1 <- pairwise.t.test(with(shan, eval(parse(text=j))), shan$group.id, "bonferroni", pool.sd=T)
      test <- pairwise.wilcox.test(with(shan, eval(parse(text=j))), shan$group.id, "bonferroni")
      out[[j]] <- list("ANOVA"=summary(atw), "M-W U test"=test, "T-test"=test1, "TukeyHSD"=tuksel)
    }
    dev.off()
    
    sink(paste("./coreplots/alpha/",deparse(substitute(comm)),"/ALPHA_boxplot-stats_", group, select, ".txt", sep=""))
    print(out)
    sink() 
  }
  
  if (type=="point"){
    jpeg(paste("./coreplots/alpha/",deparse(substitute(comm)),"/ALPHA_scatter_", group, select, ".jpg", sep=""), 6000, 2100, res=800)
    par(mfrow=c(1,3), cex.axis=1.1, mar=c(4,4,3,2), font.axis=2, las=3, xaxt="n")
    
    out <- list()
    fixedx <- 0
    for (j in c("shan", "rich", "even")){#, "chao1")){
      if (is.null(fixed)){
        yax <- shan[,j]
        plot(eval(parse(text=j)) ~ jitter(as.numeric(as.character(group.id)),1), data=shan, bg=group.id, main=select, ylab=j, xlab=" ", pch=22, cex=1.5)
      } else {
        fixedx <- fixedx+1
        yax <- fixed[[fixedx]]
        plot(eval(parse(text=j)) ~ jitter(as.numeric(as.character(group.id)),1), data=shan, bg=group.id, main=select, ylab=j, xlab=" ", ylim=yax, pch=22, cex=1.5)
      }
      
      #axis(1, at=c(1:length(levels(shan$group.id))), labels=FALSE)
      #text(c(1:length(levels(shan$group.id))), par("usr")[3], font=2, 
      #     labels=c(levels(shan$group.id)), srt=45, adj=c(1, 1.5), xpd=T, offset=2.5)
      
      m <- lm(eval(parse(text=j)) ~ as.numeric(as.character(group.id)), data=shan)
      sum <- summary(m)
      abline(m$coefficients[[1]], m$coefficients[[2]])
      
      equation <- as.expression(bquote(y == .(signif(m$coefficients[[2]], 3)) ~ x + .(signif(m$coefficients[[1]], 3)) ~~~~~~~~~ R^2 == .(signif(sum$r.squared, 3))))
      
      legend(1, (max(yax)-(max(yax)-min(yax))/10), box.lty=0, bg=alpha("white", 0.75), y.intersp=0.05, cex=0.5, legend=equation)
             
             
      if (legend==TRUE){
        if (j=="rich"){
          xax <- unique(as.numeric(as.character(shan$group.id)))
          
          legend(-(max(xax)-min(xax))/4, (min(yax)-(max(yax)-min(yax))/6), 
                 x.intersp=0.2, horiz = TRUE, fill = palette(), legend = levels(shan$group.id), xpd=T, bty = 'n', pt.cex=1)
        }
      }
      
    }
    dev.off()
  }
  
  return(shan)
}


diss.boxplots <- function(dist, meta, select, group, legend=TRUE, fixed=NULL){
  palette(c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet"))
  
  if (!"dplyr"%in%row.names(installed.packages())) {install.packages("dplyr")}
  library(dplyr)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/dissimilarity/")==FALSE){
    dir.create("./coreplots/dissimilarity/")
  }
  
  if(dir.exists(paste0("./coreplots/dissimilarity/", deparse(substitute(dist))))==FALSE){
    dir.create(paste0("./coreplots/dissimilarity/", deparse(substitute(dist))))
  }
  if(dir.exists(paste0("./coreplots/dissimilarity/", deparse(substitute(dist)), "/", group))==FALSE){
    dir.create(paste0("./coreplots/dissimilarity/", deparse(substitute(dist)), "/", group))
  }
  
  g <- which(names(meta)%in%group)
  if (class(dist)=="dist"){
    sep.dist <- data.frame(sample=labels(dist), as.matrix(dist))
    row.names(sep.dist) <- labels(dist)
  }
  if (class(dist)%in%c("data.frame", "matrix")){
    sep.dist <- data.frame(sample=row.names(dist), as.matrix(dist))
    row.names(sep.dist) <- row.names(dist)
  }
  
  if (is.null(select)==FALSE){
    TW <- sep.dist[grep(select, sep.dist$sample),
                   grep(select, names(sep.dist))]
    tw <- reshape(TW,
                  varying=grep(select, names(TW), value=TRUE),
                  #varying=grep(paste("^", select, sep=""), names(TW), value=TRUE),
                  v.names="distance",
                  timevar="sample",
                  times=paste(grep(select, names(TW), value=TRUE)),
                  ids = grep(select, row.names(TW), value=TRUE),
                  #times=paste(grep(paste("^", select, sep=""), names(TW), value=TRUE)),
                  #ids = grep(paste("^", select, sep=""), row.names(TW), value=TRUE),
                  direction="long")
  } else {
    tw <- reshape(sep.dist, varying=row.names(sep.dist),
                  v.names="distance",
                  timevar="sample",
                  times=row.names(sep.dist),
                  ids = row.names(sep.dist),
                  direction="long")
  }
  
  tw <- tw[which(duplicated(tw$distance)==FALSE),]
  tw <- tw[-which(tw$distance==0),]
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    names(meta)[1] <- "sample"
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    id <- 1
  }
  if (any(!row.names(sep.dist)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  
  tw <- left_join(tw, meta[,c(id,g)])
  tw <- left_join(tw, meta[,c(id,g)], by=c("id"="sample"))
  names(tw) <- c("sample", "distance", "id", "group.sample", "group.id")
  
  tw2 <- tw[which(tw$group.sample==tw$group.id),]
  tw2$group.id <- droplevels(as.factor(tw2$group.id))#factor(tw2$group.id, levels=unique(tw2$group.id))
  
  tw2 <- tw2[order(tw2$group.id),]
  
  clr <- c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet")
  clr <- colorRampPalette(clr)(length(unique(tw2$group.id)))
  names(clr) <- as.character(unique(tw2$group.id))
  palette(clr)
  
  jpeg(paste("./coreplots/dissimilarity/", deparse(substitute(dist)), "/", group,"/DIST_boxplot_", group, select, ".jpg", sep=""), 2300, 2650, res=600)
  par(cex.axis=1.1, font.axis=2, las=3, xaxt="n")
  if (is.null(fixed)){
    boxplot(distance ~ group.id, data=tw2, boxfill=palette(), boxwex=.8, 
         main=select, ylab="dissimilarity")
  } else {
    boxplot(distance ~ group.id, data=tw2, boxfill=palette(), 
            boxwex=.8, main=select, ylab="dissimilarity", ylim=fixed)
  }
  axis(1, at=c(1:length(levels(tw2$group.id))), labels=FALSE)
  text(c(1:length(levels(tw2$group.id))), par("usr")[3], font=2, 
       labels=c(levels(tw2$group.id)), srt=45, adj=c(1, 1.5), xpd=T, offset=2.5)
  if (legend==TRUE){
    legend('top', horiz = TRUE, fill = palette(), legend = levels(tw2$group.id), bty = 'n', cex=0.5)
  }
  dev.off()
  
  atw <- aov(distance ~ group.id, data=tw2)
  tuktw <- TukeyHSD(atw)
  source("D:/functions/tukey_select.R")
  tuksel <- tukey.select(tuktw$`group.id`, 1)
  
  test1 <- pairwise.t.test(tw2$distance, tw2$group.id, "bonferroni", pool.sd=T)
  test <- pairwise.wilcox.test(tw2$distance, tw2$group.id, "bonferroni")
  
  sink(paste("./coreplots/dissimilarity/",  deparse(substitute(dist)), "/", group,"/DIST_boxplot-stats_", group, select, ".txt", sep=""))
  print(list("ANOVA"=summary(atw), "MW U test"=test, "t-test"=test1, "tukey"=tuksel))
  sink()  
  
  return(tw2)
}


