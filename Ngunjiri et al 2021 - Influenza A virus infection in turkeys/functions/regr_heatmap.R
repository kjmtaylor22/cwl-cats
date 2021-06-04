regr.heatmap <- function(comm, tax, meta, select.gen, cyt.list, group, stdv, name, reg, w, h, swap=F, set.lim=F, ref=NULL){
  
  library(dplyr)
  library(ggplot2)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/cytokines/")==FALSE){
    dir.create("./coreplots/cytokines/")
  }
  
  ra <- function(x){ #where x is a vector
    y <- sum(x, na.rm=T)
    z <- sapply(x, FUN=function(x){z1 <- x/y})
    return(z)
  }
  
  rg <- function(x){
    m <- lm(log2(x+1) ~ factor[which(!is.na(factor))], na.action = "na.omit")
    return(coefficients(m)[2])
    #if (coefficients(m)[2]<0){
    #  y <- -summary(m)[[8]]
    #} else {y <- summary(m)[[8]]}
    #return(y)
  }
  
  pv <- function(x){
    m <- lm(log2(x +1) ~ factor[which(!is.na(factor))], na.action = "na.omit")
    return(summary(m)[[4]][8])
  }
  
  meta[,group] <- droplevels(meta[,group])
  
  g <- levels(meta[,group])
  
  #comm <- apply(comm, MARGIN=1, FUN=ra) %>% t() %>% as.data.frame()
  
  for (i in dim(tax)[2]){
    if(all(select.gen%in%tax[,i])){find <- names(tax)[i];break}else{next}
  }
  
  otu <- tax$tag[tax[,find]%in%select.gen] %>% .[.%in%names(comm)]
  
  select.gen <- gsub("-", "\\.", select.gen)
  
  data <- data.frame(tag=otu, t(comm[,otu])) %>%
    left_join(tax[,c("tag", find)]) %>%
    .[,-1] %>%
    group_by_(find) %>%
    summarize_all(funs(sum)) %>%
    as.data.frame() %>%
    `rownames<-` (as.character(.[,find])) %>%
    .[,-1] %>% t() %>% as.data.frame() %>%
    +1 %>% log2() %>%
    data.frame(SampleID=row.names(.), .) %>% 
    left_join(meta[,c("SampleID", group, cyt.list)])

  out1 <- data.frame()
  out2 <- data.frame()
  
  if (!is.null(ref)){g <- g[-grep(ref, g)]}
  
  for (i in g){
    
    sub <- data[data[,group]==i, select.gen]
    
    if (!is.null(ref)){sub <- rbind(sub, data[data[,group]==ref, select.gen])}
    
    for (j in cyt.list){
      
      factor <- data[data[,group]==i, j]
      
      if (!is.null(ref)){factor <- c(factor, data[data[,group]==ref, j])}
      
      if (!is.null(stdv)){factor[factor > stdv*sd(data[,j], na.rm=T)] <- NA}
      
      if (reg=="slope"){
      
        R1 <- summarize_all(as_tibble(sub[which(!is.na(factor)),]), funs(rg)) %>%
          as.data.frame() %>% 
          data.frame(., group=i, factor=j)
        
  
        out1 <- rbind(out1, R1)
        
        R2 <- summarize_all(as_tibble(sub[which(!is.na(factor)),]), funs(pv)) %>%
          as.data.frame() %>%
          data.frame(., group=i, factor=j)
        
        out2 <- rbind(out2, R2)
      
      }
      
      if (reg=="cor"){
        
        R2 <- psych::corr.test(sub[which(!is.na(factor)),], factor[which(!is.na(factor))], method="spearman", adjust="none")
        
        out1 <- rbind(out1, data.frame(t(R2$r), group=i, factor=j))
        out2 <- rbind(out2, data.frame(t(R2$p), group=i, factor=j))
        
      }
      
    }
    
  }

  R <- reshape2::melt(out1, variable.name=find)
  p <- reshape2::melt(out2, value.name="signif", variable.name=find)
  
  #if (reg=="cor"){p$signif <- 1-((1-p$signif)/(length(cyt.list)*length(g)))}
  #if (reg=="slope"){p$signif <- p.adjust(p$signif, "bonferroni")}
  
  p$signif <- p.adjust(p$signif, "bonferroni")
  
  signif <- p$signif
  signif[signif<=0.05] <- "**"
  signif[signif<0.1&signif>0.05] <- " " # "*"
  signif[signif>=0.1|is.nan(signif)] <- " "
  
  out <- data.frame(R, p=p$signif, signif=signif)
  
  out$value[is.nan(out$value)] <- 0
  
  tax[,find] <- gsub("-", "\\.", tax[,find])
  
  out <- data.frame(out, phylum=tax$phylum[match(out[,find], tax[,find])])
  
  out[,find] <- as.character(out[,find])
  
  out <- out[order(out$phylum, out[,find], decreasing=T),]

  out[,find] <- factor(out[,find], levels=unique(out[,find]))

  if (set.lim==F){
    lims <- c(-max(abs(out$value), na.rm=T), 
              max(abs(out$value), na.rm=T))
  } else {lims <- c(-set.lim, set.lim)}
  
  if (any(out$value>0.7, na.rm=T)){sig.col <- "white"} else {sig.col <- "black"}
  
  if (swap==T){
    
    hm <- ggplot(out, aes(group, eval(parse(text=find)))) +
      facet_wrap(vars(factor), nrow=1, strip.position = "top") +
      #geom_raster(aes(fill=value)) +
      geom_tile(fill="white", color="grey60") +
      geom_point(aes(size=abs(value), color=value), shape=19) +
      geom_text(aes(label=signif), color=sig.col) +
      scale_color_gradientn(colors=rev(c("red3", "red", "orange", "white", "turquoise", "blue", "navy")), 
                            guide=guide_colorbar(title=element_blank(), ticks.colour="black", frame.colour="black"),
                            limits=lims) +
      scale_size_continuous(range=c(0,6), limits=c(0,max(lims))) +
      theme(axis.ticks.y=element_blank(),
            axis.ticks.x=element_line(color="grey60"),
            axis.ticks.length.x = unit(0.25, "lines"),
            axis.text.x=element_text(angle=45, size=10, hjust=1, face="bold"),
            axis.text.y=element_text(size=10, vjust=0.5),
            strip.text=element_text(size=12, face="bold"),
            strip.background = element_rect(fill='white', color="white"),
            panel.border = element_rect(color="grey60", fill=NA),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10),
            legend.key.height =unit(0.15, "npc"),
            legend.key.width =unit(1, "lines")) +
      scale_y_discrete(position="left") +
      guides(size="none") +
      xlab(NULL) +
      ylab(NULL)
  
  } else {
    
    hm <- ggplot(out, aes(factor, eval(parse(text=find)))) +
      facet_wrap(vars(group), nrow=1, strip.position = "top") +
      #geom_raster(aes(fill=value)) +
      geom_tile(fill="white", color="grey60") +
      geom_point(aes(size=abs(value), color=value), shape=19) +
      geom_text(aes(label=signif), nudge_y=-0.12, color=sig.col) +
      scale_color_gradientn(colors=rev(c("red3", "tomato", "orange", "white", "turquoise", "royalblue", "navy")), 
                            guide=guide_colorbar(title=element_blank(), ticks.colour="black", frame.colour="black"),
                            limits=lims) +
      scale_size_continuous(range=c(0,6), limits=c(0,max(lims))) +
      theme(axis.ticks.y=element_blank(),
            axis.ticks.x=element_line(color="grey60"),
            axis.ticks.length.x = unit(0.25, "lines"),
            axis.text.x=element_text(angle=45, size=10, hjust=1, face="bold"),
            axis.text.y=element_text(size=10, vjust=0.5),
            strip.text=element_text(size=12, face="bold"),
            strip.background = element_rect(fill='white', color="white"),
            panel.border = element_rect(color="grey60", fill=NA),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10),
            legend.key.height =unit(0.15, "npc"),
            legend.key.width =unit(1, "lines")) +
      scale_y_discrete(position="left") +
      guides(size="none") +
      xlab(NULL) +
      ylab(NULL)
  }


  write.csv(out, paste0("./coreplots/cytokines/heatmapdata", "_", name, ".csv"))
  tiff(file=paste0("./coreplots/cytokines/heatmap", "_", name, ".tiff"), width=w, height=h, pointsize=10, res=300) 
  print(hm)
  dev.off()

  #return(out)
  return(data)
}