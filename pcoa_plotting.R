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