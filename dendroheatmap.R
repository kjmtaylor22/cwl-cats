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