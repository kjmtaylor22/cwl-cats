cyt.regression <- function(comm, tax, meta, select, cyt, group, subgroup, stdv=NULL){
  
  library(dplyr)
  library(ggplot2)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/cytokines/")==FALSE){
    dir.create("./coreplots/cytokines/")
  }
  if(dir.exists(paste0("./coreplots/cytokines/", cyt))==FALSE){
    dir.create(paste0("./coreplots/cytokines/", cyt))
  }
  
  g <- unique(meta[,group])
  sg <- unique(meta[,subgroup])
  
  assoc.list2 <- list()
  
  for (i in g){
    assoc.list2[[i]] <- select
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n+1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  clr <- gg_color_hue(length(select))
  names(clr) <- sort(select)
  names(clr) <- gsub("-", "\\.", names(clr))
  
  shp <- c(7, 18, 8, 6, 0, 17)
  names(shp) <- sg

  assoc.table <- list()
  
  for (i in dim(tax)[2]){
    if(all(select%in%tax[,i])){find <- names(tax)[i];break}else{next}
  }
  
  for (i in g){
    sub <- comm[row.names(comm)%in%meta$SampleID[meta[,group]==i],
                names(comm)%in%tax$tag[tax[,find]%in%assoc.list2[[i]]]] %>%
      
      t() %>%  
      
      data.frame(tag=row.names(.), .) %>%  
      
      left_join(tax[,c("tag", find)]) %>%
      
      .[,-1] %>% 
        
      group_by_(find) %>% 
      
      summarize_all(funs(sum)) %>% 
      
      as.data.frame() %>% 
      
      `rownames<-` (as.character(.[,find])) %>% 
      
      .[,-1] %>% 
      
      t() %>% 
      
      data.frame(SampleID=row.names(.), .) %>% 
      
      reshape2::melt() %>% 
      
      left_join(meta[,c("SampleID", subgroup, cyt)])
    
    if (!is.null(stdv)){sub <- sub[which(sub[,cyt]<(stdv*sd(sub[,cyt], na.rm=T))),]}
    
    if (dim(sub)[1]!=0){
      assoc.table[[i]] <- sub
    }
    
  }
  

  assoc.df <- data.frame()
  for (i in names(assoc.table)){
    p <- ggplot(assoc.table[[i]], aes(eval(parse(text=cyt)), log2(value+1), shape=eval(parse(text=subgroup)), color=variable)) + 
      geom_point(size=2) +
      facet_wrap(vars(variable), ncol=4) +
      scale_color_manual(values=clr) +
      scale_shape_manual(values=shp) +
      theme_bw() + 
      theme(panel.grid.minor.x = element_line(color="white"),
            panel.grid.minor.y = element_line(color="white"),
            panel.grid.major.x = element_line(linetype = "dashed", size=0.5),
            panel.grid.major.y = element_line(linetype = "dashed", size=0.5)) +
      stat_smooth(aes(eval(parse(text=cyt)), log2(value+1)), method="lm", se=T, size=0.5, inherit.aes = F) +
      ggpmisc::stat_poly_eq(inherit.aes = F, formula=y ~ x, parse=T, label.x = "left", label.y="top", size=2,
                            aes(x=eval(parse(text=cyt)), y=log2(value+1), label=paste(stat(eq.label), stat(rr.label), sep="~~~~"))) +
      ggpmisc::stat_fit_glance(inherit.aes=F, method = 'lm', method.args = list(formula = y ~ x), label.x = 'right', label.y = "top", size = 2,
                               aes(x=eval(parse(text=cyt)), y=log2(value+1), label = paste0("p = ", signif(..p.value.., digits = 2)))) +
      guides(color=F, shape=guide_legend(title=subgroup)) +
      xlab(bquote(paste("log"[.(2)]," at 5 dpi"))) + 
      ylab(bquote(paste("log"[.(2)],"(Abundance)")))
    
    tiff(file=paste0("./coreplots/cytokines/", cyt, "/", i, length(select), ".tiff"), width=2500, 
         height=(450*ceiling(length(unlist(unique(assoc.table[[i]]["variable"])))/4)), pointsize=10, res=250)  
    print(p)
    dev.off()
  
    tmp <- data.frame(assoc.table[[i]], cyt=cyt)
    names(tmp) <- c("SampleID","taxon","abund",subgroup,"FC","cyt")
    
    assoc.df <- rbind(assoc.df, tmp)
    
  }
  
  return(assoc.df)
}