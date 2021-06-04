
library(dplyr) # it's useful to have these two packages loaded from the start
library(ggplot2)


## run initialization (see initialization template script for details)

source("initialize.R")

bibtex::write.bib(c("vegan","stats","pairwiseAdonis","GUniFrac","psych"), file = "rpackages.bib")

## load necessary data files and functions

load("feature_table.RD") # initialized feature table

load("metadata.RD") # formatted metadata

load("taxonomy.RD") # initialized and parsed taxonomy

load("GenUniFrac_unweighted.RD") # unweighted UniFrac distance matrix

load("GenUniFrac_weighted.RD") # weighted UniFrac distance matrix

source("functions.R") # access to all of the functions used in this script

first.str <- function(x){
  y <- strsplit(as.character(x), split=" ")
  z <- vector("character", length = length(y))
  for (i in 1:length(z)){z[i]<-y[[i]][1]}
  return(z)
}

tax <- data.frame(tax, taxonomy2=first.str(tax$taxonomy)) # change lowest taxonomic level to minimum genus


## identify predominant taxa

comm5 <- comm[row.names(comm)%in%meta$SampleID[meta$DPI=="5DPC" & meta$Infection!="NEG"],] # subset community table to 5dpi samples and remove contaminants

common5.m <- core.id(comm5[row.names(comm5)%in%meta$SampleID[meta$Infection=="Mock"],], meta, tax, "BodySite", "Infection", 1) # isolate ASVs in 100% of Mock birds

mb5.m <- core.stack(comm5, common5.m, tax, meta, "BodySite", "Infection", "tag.name",
                  0.01, 1, "5dpiMOCKtoALL", core.only=T, fixed=T, landscape=F, view.all=T, facet=F) # project those ASVs across other treatments and set abundance threshold to 1%

common5.c <- core.id(comm5[row.names(comm5)%in%meta$SampleID[meta$Infection=="CKPA"],], meta, tax, "BodySite", "Infection", 1) # isolate ASVs in 100% of CKPA birds

mb5.c <- core.stack(comm5, common5.c, tax, meta, "BodySite", "Infection", "tag.name",
                  0.01, 1, "5dpiCKPAtoALL", core.only=T, fixed=T, landscape=F, view.all=T, facet=F) # project those ASVs across other treatments and set abundance threshold to 1%

common5.t <- core.id(comm5[row.names(comm5)%in%meta$SampleID[meta$Infection=="TKMN"],], meta, tax, "BodySite", "Infection", 1) # isolate ASVs in 100% of TKMN birds

mb5.t <- core.stack(comm5, common5.t, tax, meta, "BodySite", "Infection", "tag.name",
                  0.01, 1, "5dpiTKMNtoALL", core.only=T, fixed=T, landscape=F, view.all=T, facet=F) # project those ASVs across other treatments and set abundance threshold to 1%


mb.asv <- unique(c(as.character(mb5.m$otu), as.character(mb5.c$otu), as.character(mb5.t$otu)))

taxes <- unique(as.character(tax$taxonomy2[match(mb.asv, tax$tag.name)])) # get genus names for the predominant ASVs from all three subsets
save(taxes, file="predominanttaxa.RD") ## save this file for later

asvs <- tax$tag[tax$taxonomy2%in%taxes]
save(asvs, file="predominantasvs.RD") ## save this file for later


## alpha diversity (observed ASVs) for each body site

comm14 <- comm[row.names(comm)%in%meta$SampleID[meta$DPI=="14DPC" & meta$Infection!="NEG"],] # subset community table to 5dpi samples and remove contaminants

richness <- rbind(data.frame(rbind( # bind together alpha diversity output tables for 5dpi
                                   data.frame(alpha.boxplots(comm5, meta, "SW", "Infection", F), BodySite="NAS"),
                                   data.frame(alpha.boxplots(comm5, meta, "TW", "Infection", F), BodySite="TRA"),
                                   data.frame(alpha.boxplots(comm5, meta, "LL", "Infection", F), BodySite="LRT"),
                                   data.frame(alpha.boxplots(comm5, meta, "CE", "Infection", F), BodySite="CEC"),
                                   data.frame(alpha.boxplots(comm5, meta, "IL", "Infection", F), BodySite="ILE")), dpi="5DPC"),
                  data.frame(rbind(data.frame(alpha.boxplots(comm14, meta, "SW", "Infection", F), BodySite="NAS"),
                                   data.frame(alpha.boxplots(comm14, meta, "TW", "Infection", F), BodySite="TRA"), # bind together alpha diversity output tables for 5dpi
                                   data.frame(alpha.boxplots(comm14, meta, "LL", "Infection", F), BodySite="LRT"),
                                   data.frame(alpha.boxplots(comm14, meta, "CE", "Infection", F), BodySite="CEC"),
                                   data.frame(alpha.boxplots(comm14, meta, "IL", "Infection", F), BodySite="ILE")), dpi="14DPC"))

rich <- ggplot(subset(richness, subset=dpi=="5DPC"), aes(group.id, rich, fill=group.id)) + # print the graph to the device - make sure your variable names are correct!
          facet_grid(dpi ~ BodySite) + geom_boxplot(show.legend=F) + ggthemes::theme_few() + ylab("Observed ASVs") + xlab("") + 
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title=element_text(size=10),
                axis.text.y = element_text(angle=90, size=7, hjust=0.25, vjust=0.5, face="bold"), 
                panel.spacing.x = unit(0.8, "lines"), plot.margin = margin(0,3,0,0), strip.text = element_blank()) +
          scale_y_continuous(limits=c(0,150))


## beta diversity plots for 5dpi samples for each body site

bdist5 <- bdistWd[row.names(bdistWd)%in%meta$SampleID[meta$DPI=="5DPC"], # subset weighted UniFrac distance matrix to only 5dpi samples
                  row.names(bdistWd)%in%meta$SampleID[meta$DPI=="5DPC"]]

bdist5.wd <- list() # initialize list of weighted UniFrac distance matrices subsetted by body site
for (i in c("SW", "TW", "LL", "CE", "IL")){bdist5.wd[[i]] <- bdist5[grep(i, row.names(bdist5)), grep(i, row.names(bdist5))]}

save(bdist5.wd, file="GUniFrac_weighted_5dpi.RD") ## save this file

pcoas <- rbind(data.frame(pcoa.centroids(bdist5.wd$SW, meta, "Infection", "default", F), BodySite="NAS"), # draw PCOA and output statistical analysis for each body site matrix
               data.frame(pcoa.centroids(bdist5.wd$TW, meta, "Infection", "default", F), BodySite="TRA"),
               data.frame(pcoa.centroids(bdist5.wd$LL, meta, "Infection", "default", F), BodySite="LRT"),
               data.frame(pcoa.centroids(bdist5.wd$CE, meta, "Infection", "default", F), BodySite="CEC"),
               data.frame(pcoa.centroids(bdist5.wd$IL, meta, "Infection", "default", F), BodySite="ILE")) %>% data.frame(DPI="5DPC")


## beta diversity plots for 14dpi samples for each body site

bdist14 <- bdistWd[row.names(bdistWd)%in%meta$SampleID[meta$DPI=="14DPC"], # subset weighted UniFrac distance matrix to only 14dpi samples 
                   row.names(bdistWd)%in%meta$SampleID[meta$DPI=="14DPC"]]

bdist14.wd <- list() # initialize list of weighted UniFrac distance matrices subsetted by body site
for (i in c("SW", "TW", "LL", "CE", "IL")){bdist14.wd[[i]] <- bdist14[grep(i, row.names(bdist14)), grep(i, row.names(bdist14))]} # perform subsetting

save(bdist14.wd, file="GUniFrac_weighted_14dpi.RD") ## save this file

pcoas <- rbind(pcoas,
               rbind(data.frame(pcoa.centroids(bdist14.wd$SW, meta, "Infection", "default", F), BodySite="NAS"), # draw PCOA and output statistical analysis for each body site matrix
                     data.frame(pcoa.centroids(bdist14.wd$TW, meta, "Infection", "default", F), BodySite="TRA"),
                     data.frame(pcoa.centroids(bdist14.wd$LL, meta, "Infection", "default", F), BodySite="LRT"),
                     data.frame(pcoa.centroids(bdist14.wd$CE, meta, "Infection", "default", F), BodySite="CEC"),
                     data.frame(pcoa.centroids(bdist14.wd$IL, meta, "Infection", "default", F), BodySite="ILE")) %>% data.frame(DPI="14DPC"))

pc <- ggplot(subset(pcoas, subset=DPI=="5DPC")) + 
  facet_wrap(vars(DPI, BodySite), ncol=5, scales = "free") + # print the graph to the device - make sure your variable names are correct!
  geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
  geom_vline(xintercept=0, color="grey80", linetype="dotdash") +
  geom_segment(aes(x=Axis.1, xend=centroid.1, y=Axis.2, yend=centroid.2, color=factor), alpha=0.5, size=0.5) +
  geom_point(aes(Axis.1, Axis.2, fill=factor), alpha=0.5, shape=21, size=1.5) + 
  geom_point(aes(centroid.1, centroid.2, fill=factor), alpha=0.75, shape=24, size=2) + 
  ggthemes::theme_few() + ylab("PC2") + xlab("PC1") + 
  theme(axis.text.x = element_text(size=7, face="bold"), strip.text = element_blank(), plot.margin = margin(0,3,0,0),
        axis.text.y = element_text(angle=90, size=7, hjust=0.5, vjust=0.5, face="bold"), axis.title=element_text(size=10),
        legend.title = element_blank(), legend.position = "bottom", legend.justification = "left")

## analysis of 16S rRNA copy distribution per ng of sample DNA


copies16 <- droplevels(meta[meta$Infection!="NEG"&meta$DPI!="13DPC", c("SampleID", "ratio_16s", "DPI", "BodySite", "Infection","BodySiteInf")])
copies16$BodySite <- factor(copies16$BodySite, labels=c("Nasal","Trachea","LRT","Cecum","Ileum"))

copies <- ggplot(subset(copies16, subset=DPI=="5DPC"), aes(Infection, ratio_16s, fill=Infection)) + # print the graph to the device - make sure your variable names are correct!
              facet_grid(DPI ~ BodySite) + geom_boxplot(show.legend=F) + ggthemes::theme_few() + ylab("Copies/ng of DNA") + xlab("") + 
              theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x = unit(0.8, "lines"),
                    axis.text.y = element_text(angle=90, size=7, hjust=0.25, vjust=0.5, face="bold"), axis.title=element_text(size=10),
                    strip.text.y = element_blank(), strip.text.x = element_text(face="bold", size=14), plot.margin=margin(0,3,0,0)) +
              scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)),
                            limits = c(10^-2, 10^6))

sink("coreplots/16SCopies.txt") # create and open a .txt file
pairwise.wilcox.test(x=copies16[copies16$DPI=="14DPC", "ratio_16s"], # print this output to the text file
                     g=copies16[copies16$DPI=="14DPC", "BodySiteInf"],
                     p.adjust.method = "holm")
pairwise.wilcox.test(x=copies16[copies16$DPI=="5DPC", "ratio_16s"], # print this output to the text file
                     g=copies16[copies16$DPI=="5DPC", "BodySiteInf"],
                     p.adjust.method = "holm")
sink() # close the file

jpeg("coreplots/Figure3.jpeg", width=4300, height=2900, res=600)
egg::ggarrange(copies, rich, pc, ncol=1, padding=unit(0.1, "line"))
dev.off()

## quantify changes in abundance for predominant genera

dysb.5 <- (apply(comm5, 1, ra)*100) %>% t(.) %>% .[,asvs] %>%
  facet.updown(tax, meta, "BodySite", "Infection", "Mock", # draw plot and output statistics for 5dpi data
               filename="5dpi_from5dpiMOCK_ra", level="taxonomy2", w=2180, h=1500, set.lim = 6) %>%
  group_by(group, Infection) %>% summarize_at(vars(Log2.Abundance), funs(mean)) %>% .[-which(.[,"Infection"]=="Mock"),] # this small object shows the mean increase/decrease for affected genera

dysb.14 <- (apply(comm14, 1, ra)*100) %>% t(.) %>% .[,asvs] %>%
  facet.updown(tax, meta, "BodySite", "Infection", "Mock", # draw plot and output statistics for 14dpi data
               filename="14dpi_from5dpiMOCK_ra", level="taxonomy2", w=2180, h=1500, set.lim = 6) %>%
  group_by(group, Infection) %>% summarize_at(vars(Log2.Abundance), funs(mean)) %>% .[-which(.[,"Infection"]=="Mock"),] # this small object shows the mean increase/decrease for affected genera

inf.5 <- (apply(comm5, 1, ra)*100) %>% t(.) %>% .[row.names(.)%in%meta$SampleID[meta$Infection!="Mock"],asvs] %>%
  facet.updown(tax, meta, "BodySite", "Infection", "TKMN", # draw plot and output statistics for 5dpi data
               filename="5dpi_from5dpiINF_ra", level="taxonomy2", w=2180, h=1500, set.lim = 6) 

inf.14 <- (apply(comm14, 1, ra)*100) %>% t(.) %>% .[row.names(.)%in%meta$SampleID[meta$Infection!="Mock"],asvs] %>%
  facet.updown(tax, meta, "BodySite", "Infection", "Mock", # draw plot and output statistics for 14dpi data
               filename="14dpi_from5dpiINF_ra", level="taxonomy2", w=2180, h=1500, set.lim = 6)

## visualize correlation of predominant genera to gene expression

cyts <- c("UT_IFN.a", "UT_IFN.b","UT_IFN.g","UT_IFN.l","UT_OAS","UT_Mx","UT_IL.6","UT_LITAF") # names of metadata columns containing gene expression data
ct.cyts <- c("CT_IFN.g","CT_IFN.l", "CT_OAS", "CT_Mx","CT_IL4","CT_IL5","CT_IL6","CT_IL10","CT_IL12","CT_LITAF") # names of metadata columns containing gene expression data

for (c in cyts){ # go through cyts iteratively; they will save to separate folders
  cyt.regression(comm5[row.names(comm5) %in% meta$SampleID[meta$BodySite == "TRA"],], # plots show scatterplots and regressions with equations 
                 tax, meta[meta$DPI == "5DPC" & meta$BodySite == "TRA" & meta$Infection != "NEG",], select=taxes, c, "BodySite", "Infection")
}



mar <- robust.threshold(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite!="LRT"],], meta, "BodySite", "Infection", asvs) # check robustness of margin choice

common <- core.id(comm5[,asvs], meta, tax, "BodySite", "Infection", 0.01) # presence-absence distribution of predominant genera

concat <- list("CEC"=list("Mock"=unique(unlist(common$CEC))), "ILE"=list("Mock"=unique(unlist(common$ILE))), # have to trick the function into using the same ASV list for all treatments
               "NAS"=list("Mock"=unique(unlist(common$NAS))), "TRA"=list("Mock"=unique(unlist(common$TRA))))

mb <- core.stack(comm5, concat, tax, meta, "BodySite", "Infection", "taxonomy2", # this data table describes the consensus distribution of predominant genera present in each site
                 0.005, 0.01, "_predominanttaxa", core.only=T, fixed=T, landscape=T, view.all=T, facet=T)

regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite=="TRA" & meta$Infection!="NEG"], ], # this produces the correlation heatmap summarizing the regressions discussed above
                    tax, meta[meta$BodySite=="TRA" & meta$DPI=="5DPC",], as.character(unique(mb$otu[mb$groupvar=="TRA"])), 
                    cyts, "BodySite", 10, "UT_cytokines_TRAonly", "cor", 1200, 1200, F, 1, NULL)
regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite=="TRA" & !meta$Infection%in%c("NEG","TKMN")], ], # this produces the correlation heatmap summarizing the regressions discussed above
             tax, meta[meta$BodySite=="TRA" & meta$DPI=="5DPC",], as.character(unique(mb$otu[mb$groupvar=="TRA"])), 
             cyts, "BodySite", 10, "UT_cytokines_TRAonly_CKPA", "cor", 1200, 1200, F, 1, NULL)
regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite=="TRA" & !meta$Infection%in%c("NEG","CKPA")], ], # this produces the correlation heatmap summarizing the regressions discussed above
             tax, meta[meta$BodySite=="TRA" & meta$DPI=="5DPC",], as.character(unique(mb$otu[mb$groupvar=="TRA"])), 
             cyts, "BodySite", 10, "UT_cytokines_TRAonly_TKMN", "cor", 1200, 1200, F, 1, NULL)

pull.sig <- read.csv("coreplots/cytokines/heatmapdata_UT_cytokines_TRAonly.csv")

write.csv(pull.sig %>% reshape2::dcast(taxonomy2 ~ factor, value.var = "value"), "../TableS1.csv")

sig.only <- data.frame()
for (c in cyts[c(2,3,4,5,7)]){ # go through cyts iteratively; they will save to separate folders
  sig.only <- rbind(sig.only,
    cyt.regression(comm5[row.names(comm5) %in% meta$SampleID[meta$BodySite == "TRA"],], # plots show scatterplots and regressions with equations 
                   tax, meta[meta$DPI == "5DPC" & meta$BodySite == "TRA" & meta$Infection != "NEG",], 
                   select=taxes[c(5,13,18,21:23)], c, "BodySite", "Infection")
  )
}


ggplot(sig.only, aes(abund, FC, color=Infection, shape=Infection)) + # print the graph to the device - make sure your variable names are correct!
  facet_grid(cyt ~ taxon) + geom_point() + ggthemes::theme_few() + ylab("Fold Change") + xlab("log2 abund at 5dpi") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x = unit(0.8, "lines"),
        axis.text.y = element_text(angle=90, size=7, hjust=0.25, vjust=0.5, face="bold"))


## visualize correlation of predominant genera to viral titer

sites <- c("NAS", "TRA", "CEC", "ILE") # names of sites with viral titer data attached; LRT is exempt (see Materials & Methods)

for (s in sites){ # these plots will all go in the same folder
  cyt.regression(comm5[row.names(comm5) %in% meta$SampleID[meta$BodySite == s],], # plots show scatterplots and regressions with equations  
                 tax, meta[meta$DPI == "5DPC" & meta$BodySite == s & meta$Infection != "NEG",], select=taxes, "EID50.mL", "BodySite", "Infection")
  cyt.regression(comm5[row.names(comm5) %in% meta$SampleID[meta$BodySite == s],], # plots show scatterplots and regressions with equations  
                 tax, meta[meta$DPI == "5DPC" & meta$BodySite == s & meta$Infection != "NEG",], select=taxes, "BodyWt_g", "BodySite", "Infection")
}

regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite!="LRT"], asvs], tax[tax$tag%in%asvs,], 
             meta[meta$BodySite!="LRT" & meta$DPI=="5DPC",], 
             taxes, "EID50.mL", "BodySite", 10, "ALL_titers", "cor", 975, 1700, T, 1, NULL) # this produces the correlation heatmap summarizing the regressions discussed above
regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite!="LRT"], asvs], tax[tax$tag%in%asvs,],
             meta[meta$BodySite!="LRT" & meta$DPI=="5DPC" & !meta$Infection%in%c("NEG","TKMN"),], 
             taxes, "EID50.mL", "BodySite", 10, "ALL_titers_CKPA", "cor", 975, 1700, T, 1, NULL)
regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite!="LRT"], asvs], tax[tax$tag%in%asvs,],
             meta[meta$BodySite!="LRT" & meta$DPI=="5DPC" & !meta$Infection%in%c("NEG","CKPA"),], 
             taxes, "EID50.mL", "BodySite", 10, "ALL_titers_TKMN", "cor", 975, 1700, T, 1, NULL)


regr.heatmap(comm5[, asvs], tax[tax$tag%in%asvs,], meta[meta$DPI=="5DPC"&meta$Infection!="NEG",], 
             taxes, "BodyWt_g", "BodySite", NULL, "ALL_weights", "cor", 1025, 1700, T, 1, NULL)  # this produces the correlation heatmap summarizing the regressions discussed above
regr.heatmap(comm5[, asvs], tax[tax$tag%in%asvs,], meta[meta$DPI=="5DPC" & !meta$Infection%in%c("NEG","TKMN"),], 
             taxes, "BodyWt_g", "BodySite", NULL, "ALL_weights_CKPA", "cor", 1025, 1700, T, 1, NULL)
regr.heatmap(comm5[, asvs], tax[tax$tag%in%asvs,], meta[meta$DPI=="5DPC" & !meta$Infection%in%c("NEG","CKPA"),], 
             taxes, "BodyWt_g", "BodySite", NULL, "ALL_weights_TKMN", "cor", 1025, 1700, T, 1, NULL)


meta.cyt <- meta[meta$BodySite=="TRA" & meta$DPI=="5DPC" & meta$Infection!="NEG", c("SampleID", cyts)] %>%
  reshape2::melt() %>% left_join(meta[,c("SampleID","EID50.mL","Infection")]) %>%
  mutate(., CKPA=value, TKMN=value)
meta.cyt$TKMN[meta.cyt$Infection=="CKPA"] <- NA
meta.cyt$CKPA[meta.cyt$Infection=="TKMN"] <- NA


ggplot(meta.cyt, aes(EID50.mL, value, shape=Infection, color=Infection)) + 
  geom_point(size=2) +
  facet_wrap(vars(variable), ncol=4) +
  scale_shape_manual(values=c(7,18,8)) +
  theme_bw() +   ylab("Fold Change") +
  theme(panel.grid.minor.x = element_line(color="white"),
        panel.grid.minor.y = element_line(color="white"),
        panel.grid.major.x = element_line(linetype = "dashed", size=0.5),
        panel.grid.major.y = element_line(linetype = "dashed", size=0.5)) +
  stat_smooth(aes(EID50.mL, CKPA), method="lm", se=T, size=0.5, inherit.aes = F, color="blue") +
  ggpmisc::stat_poly_eq(inherit.aes = F, formula=y ~ x, parse=T, label.x = "left", label.y=1, size=2.5,
                        aes(x=EID50.mL, y=CKPA, label=paste(stat(eq.label), stat(rr.label), sep="~~~~")), color="blue") +
  stat_smooth(aes(EID50.mL, TKMN), method="lm", se=T, size=0.5, inherit.aes = F, color="chartreuse4") +
  ggpmisc::stat_poly_eq(inherit.aes = F, formula=y ~ x, parse=T, label.x = "left", label.y=0.95, size=2.5,
                        aes(x=EID50.mL, y=TKMN, label=paste(stat(eq.label), stat(rr.label), sep="~~~~")), color="chartreuse4") 
ggsave("coreplots/correlation_EID50_cytokines.jpeg", width=9, height=6, dpi=300)


write.csv(read.csv("coreplots/cytokines/heatmapdata_ALL_weights.csv") %>% reshape2::dcast(taxonomy2 ~ group, value.var = "value"),
          "coreplots/correlation_weights.csv")
## produce NMDS of common predominant genera per body site and k-means clusters

ls <- list() # initialize list to store NMDS and k-means output
for (i in sites){ # iterate across body sites
  comm.rmv <- comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite==i], asvs] %>% # community table made up of ASVs ascribed to predominant genera
    t() %>% data.frame(tag=row.names(.)) %>% right_join(tax[,c("tag", "taxonomy2")], .) %>% .[,-1] %>% # match genus names to ASV tags
    group_by(taxonomy2) %>% summarize_all(funs(sum)) %>% as.data.frame() %>% # add read counts of that genus together
    .[which(.[,"taxonomy2"]%in%as.character(unique(mb$otu[mb$groupvar==i]))),] %>% `row.names<-` (.[,"taxonomy2"]) %>% .[,-1] # subset genera down to those predominant in the current body site 
  k <- kmeans(vegan::metaMDS(t(comm.rmv),"bray", trymax = 1200, autotransform = F)[[35]], 3, 100, 20) # perform NMDS and get k-means clusters out of MDS1 and MDS2
  m <- cbind(vegan::metaMDS(t(comm.rmv),"bray", trymax = 1200, autotransform = F)[[35]], k$cluster) %>% data.frame(taxonomy2=row.names(.)) # bind NMDS axis data to k-means clusters
  ls[[i]] <- m # save to list
}

titers <- read.csv("coreplots/cytokines/heatmapdata_ALL_titers.csv") # this data table was output by the heatmap function

titers$taxonomy2 <- gsub("\\.", "-", titers$taxonomy2) # must adjust the characters so genus names match between tables

kmeans <- left_join(titers, rbind(data.frame(ls$NAS, group="NAS"), data.frame(ls$TRA, group="TRA"), # join the viral titer correlation data to k-means cluster 
                                  data.frame(ls$CEC, group="CEC"), data.frame(ls$ILE, group="ILE"))) %>% .[!is.na(.[,"V3"]),] # genera not common in a particular site are removed
  
k1 <- kmeans %>% group_by(group, V3) %>% summarize_at(vars(value), funs(mean, length)) # summarize correlation data for each k-means group

k1 <- k1[order(k1$group, k1$mean),] %>% data.frame(reorder=rep(1:3,4)) # reorder group designations per body site based on mean correlation value for that group

kmeans <- left_join(kmeans, k1[,c(1,2,4,5)]) # join new group designations to the kmeans table

kmeans$reorder <- factor(kmeans$reorder, labels=c("C", "B", "A")) # rename the grouping factor

kmeans <- mutate(kmeans, GroupGen=paste(reorder, taxonomy2, sep="_")) %>% .[order(.$group, .$reorder, .$taxonomy2),] %>%
  data.frame(Y=c(1:dim(kmeans)[1])) %>% group_by(group, reorder) %>% mutate(., mean=mean(Y)) # this is for reshaping the graph to be more aesthetic

titers <- reshape2::dcast(kmeans, taxonomy2 ~ group, value.var = "value") # reshape the table into wide form

write.csv(titers, "coreplots/correlation viraltiters.csv") # save as output table

save(kmeans, file="kmeans_clusters.RD") # save for later

jpeg("coreplots/titer_correlationchart.jpeg", width=2000, height=2000, res=300) # open the graphics device
ggplot(kmeans, aes(Y, value, fill=reorder, group=taxonomy2, label=taxonomy2)) + facet_wrap(vars(group), nrow=2, dir="v", scales="free_y") + # print graph to device
  geom_col(position=position_dodge(width=0.95), color="white", width=0.95) +
  geom_text(aes(label=signif, hjust=scales::rescale(value/abs(value))), size=8, nudge_x = -0.3) +
  geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
  geom_hline(yintercept=-0.5, color="grey80", linetype="dotted") +
  geom_hline(yintercept=0.5, color="grey80", linetype="dotted") +
  scale_x_continuous(breaks=kmeans$mean, labels=kmeans$reorder) +
  geom_text(aes(y=0, hjust=scales::rescale(value/abs(value))), 
            vjust=0.25, position=position_dodge(width=0.95), size=3) +
  ylab("Spearman correlation\n(relative abundance vs. viral titer)") + xlab("k-means group\n(NMDS-derived)") +
  coord_flip() + guides(fill=F) + ggthemes::theme_few() + theme(axis.text.y = element_text(face="bold", size=12))
dev.off() # close the device


## class plots (Bacilli and Gammaproteobacteria)

rel.comm <- apply(comm, MARGIN=1, FUN=ra) %>% t() %>% as.data.frame()
save(rel.comm, file="relativeabundancecomm.RD")

dpi5 <- rel.comm[row.names(rel.comm)%in%meta$SampleID[meta$DPI=="5DPC" & meta$Infection!="NEG"], match(asvs,names(rel.comm))] %>%
  data.frame(SampleID=row.names(.)) %>% reshape2::melt(id.vars="SampleID", variable.name="tag") %>% 
  left_join(meta[,c("SampleID", "BodySite", "Infection")]) %>% left_join(tax[,c("tag", "genus", "class")]) %>%
  group_by(SampleID, BodySite, Infection, class) %>% summarize_at(vars(value), funs(sum)) %>% as.data.frame() %>%
  .[which(.[,"class"]%in%c("Bacilli", "Gammaproteobacteria")),] %>%
  mutate(., groupnames=paste(BodySite, Infection, sep="_"))
jpeg("./coreplots/class_plot_ALL_5dpi.jpeg", width=4500, height=2200, res=600)
ggplot(dpi5) + geom_jitter(aes(Infection, value*100, color=Infection), width=0.2) + guides(color=F) +
  facet_grid(class ~ BodySite) + ggthemes::theme_few() + ylab("Relative Abundance (%)") + xlab("") + ylim(0,100) +
  stat_summary(fun.y=mean, geom="errorbar", aes(Infection, value*100, ymax=..y.., ymin=..y..), width=.75, linetype="solid", size=.8)
dev.off()
sink("./coreplots/class_plot_ALL_5dpi.txt")
pairwise.wilcox.test(dpi5$value[dpi5$class=="Bacilli"], dpi5$groupnames[dpi5$class=="Bacilli"], "holm")
pairwise.wilcox.test(dpi5$value[dpi5$class=="Gammaproteobacteria"], dpi5$groupnames[dpi5$class=="Gammaproteobacteria"], "holm")
sink()

dpi14 <- rel.comm[row.names(rel.comm)%in%meta$SampleID[meta$DPI=="14DPC" & meta$Infection!="NEG"], match(asvs,names(rel.comm))] %>%
  data.frame(SampleID=row.names(.)) %>% reshape2::melt(id.vars="SampleID", variable.name="tag") %>% 
  left_join(meta[,c("SampleID", "BodySite", "Infection")]) %>% left_join(tax[,c("tag", "genus", "class")]) %>%
  group_by(SampleID, BodySite, Infection, class) %>% summarize_at(vars(value), funs(sum)) %>% as.data.frame() %>%
  .[which(.[,"class"]%in%c("Bacilli", "Gammaproteobacteria")),] %>%
  mutate(., groupnames=paste(BodySite, Infection, sep="_"))
jpeg("./coreplots/class_plot_ALL_14dpi.jpeg", width=4500, height=2200, res=600)
ggplot(dpi14) + geom_jitter(aes(Infection, value*100, color=Infection), width=0.2) + guides(color=F) +
  facet_grid(class ~ BodySite) + ggthemes::theme_few() + ylab("Relative Abundance (%)") + xlab("") + ylim(0,100) +
  stat_summary(fun=mean, geom="errorbar", aes(Infection, value*100, ymax=..y.., ymin=..y..), width=.75, linetype="solid", size=.8)
dev.off()
sink("./coreplots/class_plot_ALL_14dpi.txt")
pairwise.wilcox.test(dpi14$value[dpi14$class=="Bacilli"], dpi14$groupnames[dpi14$class=="Bacilli"], "holm")
pairwise.wilcox.test(dpi14$value[dpi14$class=="Gammaproteobacteria"], dpi14$groupnames[dpi14$class=="Gammaproteobacteria"], "holm")
sink()
