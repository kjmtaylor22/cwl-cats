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
      mutate(value=log2(value+1)) %>%
      group_by_(group, subgroup, "variable") %>%
      summarize_at(vars(value), funs(mean, var, length)) %>% as.data.frame()
    names(tran)[1:3] <- c("group", "subgroup", "variable")
    
    grps <- as.character(unique(tran$subgroup))

    mn <- reshape2::dcast(tran, formula= group + variable ~ subgroup, value.var="mean") 
    vr <- reshape2::dcast(tran, formula= group + variable ~ subgroup, value.var="var") 
    df <- reshape2::dcast(tran, formula= group + variable ~ subgroup, value.var="length")
    
    #this data is not log-transformed --> the log transformation will come later
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
      #mutate(value=log2(value+1)) %>%
      group_by_(group, subgroup, "variable") %>%
      summarize_at(vars(value), funs(mean)) %>% as.data.frame()
    names(tran)[1:3] <- c("group", "subgroup", "variable")
        
    tran2 <- reshape2::dcast(tran, formula= group + variable ~ subgroup, value.var="value")
    
    final <- data.frame(group=tran2$group, variable=tran2$variable)
    
    for (i in 1:length(grps)){
      #final[[grps[i]]] <- tran2[,grps[i]]-tran2[,control.grp]
      final[[grps[i]]] <- log2((tran2[,grps[i]]+1)/(tran2[,control.grp]+1))

    }

    com <- combn(1:length(grps), 2)
    out <- data.frame()
    for (i in 1:dim(com)[2]){
      new <- mutate(mn, t=abs(mn[,grps[com[,i][1]]]-mn[,grps[com[,i][2]]])/sqrt((vr[,grps[com[,i][1]]]/(df[,grps[com[,i][1]]]-1))+(vr[,grps[com[,i][2]]]/(df[,grps[com[,i][2]]]-1)))) %>%
        mutate(comp= paste(grps[com[,i][1]],grps[com[,i][2]], sep="-"),
               p=2*pt(t, df=sum(df[,grps[com[,i][1]]],df[,grps[com[,i][2]]],-2), lower=FALSE)) 
      out <- rbind(out, new)
    }
    
    out <- mutate(out, p.adj=p.adjust(p, method = "holm")) %>% subset(subset=p.adj<0.1)
    
    
    twist <- out %>% .[grep(control.grp, .$comp),c("group","variable","comp","p.adj")] 
    twist$comp <- gsub(paste0(control.grp, "-"),"", twist$comp)
    names(twist) <- c(group, level, subgroup, "p.adj")

    final <- reshape2::melt(final) %>%
      `colnames<-` (c(group, level, subgroup, "Log2 Abundance")) %>%
      left_join(twist) %>% mutate(., signif="")
  
    final$signif[final$p.adj<0.05] <- "*"
    
    final[,subgroup] <- as.factor(final[,subgroup])
    final[,subgroup] <- factor(final[,subgroup], levels=levels(meta[,subgroup]))
    #levels(final[,subgroup])[1] <- "Control"
    
    
    
    tax$genus <- gsub("-", "\\.", tax$genus)
    
    final <- data.frame(final, phylum=tax$phylum[match(final[,level], tax[,level])])
    
    final[,level] <- as.character(final[,level])
    
    final <- final[order(final$phylum, final[,level], decreasing=T),]
    
    final[,level] <- factor(final[,level], levels=rev(unique(final[,level])))
  
  
  if (set.lim==F){
    lims <- c(-max(abs(final$Log2.Abundance)),
              max(abs(final$Log2.Abundance)))
  } else {lims <- c(-set.lim, set.lim)}
  
  p <- ggplot(final, aes(eval(parse(text=subgroup)), eval(parse(text=level)))) +
    geom_raster(aes(fill=Log2.Abundance)) +
    geom_text(aes(label=signif), size=30, nudge_y = -0.3) +
    facet_wrap(vars(eval(parse(text=group))), dir="h", nrow=1, strip.position = "top") +
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
  
  sink(paste0("./coreplots/facetupdown/", folder, "/", filename, ".txt"))
  print(out)
  sink()
  
  write.csv(final, paste0("./coreplots/facetupdown/", folder, "/", filename, ".csv"))
  
  return(final)
}