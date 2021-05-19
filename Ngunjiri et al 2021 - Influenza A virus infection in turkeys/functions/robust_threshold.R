robust.threshold <- function(comm, meta, group, subgroup, select){

  meta <- mutate(meta, !!"group.subgroup" := paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep="."))
  
  grouped <- data.frame(SampleID=row.names(comm), comm) %>% left_join(meta[,c("SampleID", "group.subgroup")]) %>%
    .[,-1] %>% group_by(group.subgroup) %>% summarize_all(funs(mean)) %>% as.data.frame() %>% `row.names<-` (.$group.subgroup) %>%
    .[,-1] %>% apply(MARGIN=1, FUN=ra) %>% t() %>% as.data.frame()
  
  ls.length <- function(x){return(length(unlist(x,F)))}
  
  new <- data.frame()
  for (k in seq(0.05,1,0.05)){
    
    common <- core.id(comm[,select], meta, group=group, subgroup=subgroup, margin = k)

    y <- c()
    for (j in names(unlist(common, F))){ 
      
     y <- c(y, sum(grouped[which(row.names(grouped)==j), unlist(common, F)[[j]]]))
      
    }
    
    o <- data.frame(group.subgroup=names(unlist(common,F)), margin=k, rel.ab=y,
                    count=unlist(lapply(unlist(common,F), FUN=ls.length)))
    
    new <- rbind(new, o)
    
  }
  
  new <- left_join(new, unique(meta[,c(group, subgroup, "group.subgroup")]))
  
  jpeg("coreplots/predominantasvs.jpeg", width=2200, height=3200, res=300)
  gridExtra::grid.arrange(
    ggplot(new, aes(margin, count, color=Infection, shape=Infection))+
      facet_wrap(vars(BodySite), dir = "v") +  geom_point(size=3) + geom_line()+ theme_bw() + 
      theme(strip.text = element_text(face="bold", size=14), legend.text = element_text(size=12))+
      scale_shape_manual(values=c(7,8,18)) +
      ylab("Number of ASVs Detected") +
      xlab(" Probability of Species Detection \n (Percentage of Host Individuals) ") +
      scale_x_continuous(breaks = seq(0,1, 0.1), labels = paste0(seq(0,100, 10), "%"), limits=c(0,1)),
    ggplot(new, aes(margin, rel.ab*100, color=Infection, shape=Infection))+
      facet_wrap(vars(BodySite), dir = "v") +  geom_point(size=3) + geom_line()+ theme_bw() + 
      theme(strip.text = element_text(face="bold", size=14), legend.text = element_text(size=12))+
      scale_shape_manual(values=c(7,8,18)) +
      ylab("Total Relative Abundance (%)") +
      xlab(" Probability of Species Detection \n (Percentage of Host Individuals) ") +
      scale_x_continuous(breaks = seq(0,1, 0.1), labels = paste0(seq(0,100, 10), "%"), limits=c(0,1))
  )
  dev.off()
  
  return(new)

}