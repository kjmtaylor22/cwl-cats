facet.updown <- function(comm, tax, meta, group, subgroup, control.grp, folder=Sys.Date(), 
                         filename=Sys.time(), level="genus", set.lim=F, w, h){
  
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/facetupdown/")==FALSE){
    dir.create("./coreplots/facetupdown/")
  } 
  if(dir.exists(paste0("./coreplots/facetupdown/", folder, "/"))==FALSE){
    dir.create(paste0("./coreplots/facetupdown/", folder, "/"))
  } 
  
  cran <- c("dplyr", "ggplot2", "reshape2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  ## this data is log-transformed to adhere to the t-test
  tran <- comm %>% t() %>% as.data.frame() %>%
    data.frame(tag=row.names(.), .) %>%
    left_join(tax[,c("tag",level)]) %>%
    .[,-1] %>% group_by_(level) %>%
    summarize_all(funs(sum)) %>% as.data.frame() %>%
    #.[-which(.[,level]%in%c("", NA)),] %>%
    `rownames<-` (as.character(.[,level])) %>%
    .[,-1] %>% t() %>% as.data.frame() %>%
    data.frame(SampleID=row.names(.), .) %>%
    left_join(meta[,c("SampleID", group, subgroup)]) %>%
    reshape2::melt(.) %>%
    group_by_(group, subgroup, "variable") %>%
    mutate(setnum=1:length(SampleID)) %>%
    droplevels(.)
    #mutate(value=log2(value+1)) 
    names(tran)[2:4] <- c("group", "subgroup", "variable")
    tran2 <- reshape2::dcast(tran, setnum + group + variable ~ subgroup, value.var = "value") 

    
  
  #ggplot(tran) + geom_histogram(aes(value, fill=eval(parse(text=subgroup))), color="white", alpha=0.5) +
  #  facet_grid(variable ~ eval(parse(text=group)), scales = "free") + ggthemes::theme_few()
  #ad.p <- function(x){if (sum(x!=0)){y <- nortest::ad.test(x)$p.value}else{y<-1}; return(y)}
  
  #normtest <- tran %>% group_by_(group, subgroup, "variable") %>%
  #  summarize_at(vars(value), funs(ad.p))
    
  grps <- levels(tran$subgroup)

  for (i in grps[2:length(grps)]){
    tran2 <- tran2 %>%
      group_by(group, variable) %>% ## within these groups, the abundance data are log-normal
      mutate("pval.{i}" := wilcox.test(x=eval(parse(text=grps[1])), y=eval(parse(text=i)))$p.value)
  }
  
  tran3 <- tran2[,-1] %>%
    group_by(group, variable) %>%
    summarize_all(funs(mean), na.rm=T)

  final <- data.frame()
  for (i in grps[2:length(grps)]){
    ctrl <- grep(grps[1], names(tran3))
    trmt <- grep(i, names(tran3))

    temp <- tran3[,c(1,2,ctrl,trmt)] %>%
      mutate("{i}" := log2((eval(parse(text=i))+1)/(eval(parse(text=grps[1]))+1))) 
    temp[,ctrl] <- 0
    names(temp)[grep("pval", names(temp))] <- "p.value"
    
    temp2 <- reshape2::melt(temp, id.vars=c("group", "variable", "p.value"), 
                            variable.name=subgroup, value.name="Log2 Abundance")
    
    final <- rbind(final, temp2)
  }
  
  final$p.adj <- p.adjust(final$p.value, method = "holm") 
  final$signif[final$p.adj<0.05] <- "*"
  names(final)[2] <- level
  
  
  tax$genus <- gsub("-", "\\.", tax$genus)
  
  final <- data.frame(final, phylum=tax$phylum[match(final[,level], tax[,level])])
  
  final[,level] <- as.character(final[,level])
  
  final <- final[order(final$phylum, final[,level], decreasing=T),]
  
  final[,level] <- factor(final[,level], levels=rev(unique(final[,level])))
  
  final$signif[final[,subgroup]==grps[1]] <- NA
  
  final[,subgroup] <- factor(final[,subgroup], levels=grps)
  
  
  if (set.lim==F){
    lims <- c(-max(abs(final$Log2.Abundance)),
              max(abs(final$Log2.Abundance)))
  } else {lims <- c(-set.lim, set.lim)}
  
  p <- ggplot(final, aes(eval(parse(text=subgroup)), eval(parse(text=level)))) +
    geom_raster(aes(fill=Log2.Abundance)) +
    geom_text(aes(label=signif), size=30, nudge_y = -0.3) +
    facet_wrap(vars(group), dir="h", nrow=1, strip.position = "top") +
    scale_fill_gradientn(colors=c("blue", "white", "red"), 
                         guide="colorbar", limits=lims, breaks=sort(c(lims, 0)),
                         labels=as.character(sort(c(round(lims), 0)))) +
    theme_bw()+
    theme(axis.ticks=element_line(color="white"), 
          axis.text.x=element_text(angle=90, size=50, hjust=1, vjust=0.5, face="bold"),
          axis.text.y=element_text(size=50, vjust=0.5),
          axis.title=element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          plot.caption=element_text(size=40),
          panel.border = element_rect(color="black"),
          strip.text = element_text(size=60, face="bold"),
          strip.background = element_rect(color = "white", fill="white"),
          legend.title=element_text(size=70, face="bold", angle=90, hjust=0.5),
          legend.text=element_text(size=60, hjust=0),
          legend.position = "left",
          legend.direction = "vertical",
          legend.justification = "right",
          legend.key.height =unit(0.15, "npc"),
          legend.key.width =unit(0.03, "npc")) +
    scale_y_discrete(position="left") +
    guides(fill=guide_colorbar(title=bquote(paste("Differential log"[.(2)]," Abundance")),
                               title.position="left", ticks=F, label.hjust = 0))
  
  jpeg(paste0("./coreplots/facetupdown/", folder, "/", filename, ".jpeg"), width=w, height=h)
  print(p)
  dev.off()

  write.csv(final, paste0("./coreplots/facetupdown/", folder, "/", filename, ".csv"))
  
  return(final)
}
