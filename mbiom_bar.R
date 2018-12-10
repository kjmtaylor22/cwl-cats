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