
library(dplyr) # it's useful to have these two packages loaded from the start
library(ggplot2)


## run initialization (see initialization template script for details)

setwd(".")

source("initialize.R")

bibtex::write.bib(c("vegan","stats","pairwiseAdonis","GUniFrac","psych"), file = "rpackages.bib")

## load necessary data files and functions

load("feature_table.RD") # initialized feature table

load("metadata.RD") # formatted metadata

load("taxonomy.RD") # initialized and parsed taxonomy

load("GenUniFrac_unweighted.RD") # unweighted UniFrac distance matrix

load("GenUniFrac_weighted.RD") # weighted UniFrac distance matrix

source("functions.R") # access to all of the functions used in this script


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

gens <- unique(as.character(tax$genus[match(unique(c(as.character(mb5.m$otu), as.character(mb5.c$otu), as.character(mb5.t$otu))), tax$tag.name)])) # get genus names for the predominant ASVs from all three subsets

any(gens=="") # if any predominant ASVs have not been assigned a genus, remove them from this list
gens <- gens[-which(gens=="")]

save(gens, file="predominantgenera.RD") ## save this file for later

asvs <- tax$tag[tax$genus%in%gens] # find all ASVs ascribed to the named predominant genera

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

jpeg("coreplots/richness.jpeg", width=5200, height=1800, res=600) # open the .jpeg device
ggplot(richness, aes(group.id, rich, fill=group.id)) + # print the graph to the device - make sure your variable names are correct!
  facet_grid(dpi ~ BodySite) + geom_boxplot() + ggthemes::theme_few() + ylab("Observed ASVs") + xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x = unit(0.8, "lines"),
        axis.text.y = element_text(angle=90, size=7, hjust=0.25, vjust=0.5, face="bold"))
dev.off() # close the device


## beta diversity plots for 5dpi samples for each body site

bdist5 <- bdistWd[row.names(bdistWd)%in%meta$SampleID[meta$DPI=="5DPC"], # subset weighted UniFrac distance matrix to only 5dpi samples
                  row.names(bdistWd)%in%meta$SampleID[meta$DPI=="5DPC"]]

bdist5.wd <- list() # initialize list of weighted UniFrac distance matrices subsetted by body site
for (i in c("SW", "TW", "LL", "CE", "IL")){bdist5.wd[[i]] <- bdist5[grep(i, row.names(bdist5)), grep(i, row.names(bdist5))]}

save(bdist5.wd, file="GUniFrac_weighted_5dpi.RD") ## save this file

pcoa.centroids(bdist5.wd$SW, meta, "Infection", "default", F) # draw PCOA and output statistical analysis for each body site matrix
pcoa.centroids(bdist5.wd$TW, meta, "Infection", "default", F)
pcoa.centroids(bdist5.wd$LL, meta, "Infection", "default", F)
pcoa.centroids(bdist5.wd$CE, meta, "Infection", "default", F)
pcoa.centroids(bdist5.wd$IL, meta, "Infection", "default", F)


## beta diversity plots for 14dpi samples for each body site

bdist14 <- bdistWd[row.names(bdistWd)%in%meta$SampleID[meta$DPI=="14DPC"], # subset weighted UniFrac distance matrix to only 14dpi samples 
                   row.names(bdistWd)%in%meta$SampleID[meta$DPI=="14DPC"]]

bdist14.wd <- list() # initialize list of weighted UniFrac distance matrices subsetted by body site
for (i in c("SW", "TW", "LL", "CE", "IL")){bdist14.wd[[i]] <- bdist14[grep(i, row.names(bdist14)), grep(i, row.names(bdist14))]} # perform subsetting

save(bdist14.wd, file="GUniFrac_weighted_14dpi.RD") ## save this file

pcoa.centroids(bdist14.wd$SW, meta, "Infection", "default", F) # draw PCOA and output statistical analysis for each body site matrix
pcoa.centroids(bdist14.wd$TW, meta, "Infection", "default", F)
pcoa.centroids(bdist14.wd$LL, meta, "Infection", "default", F)
pcoa.centroids(bdist14.wd$CE, meta, "Infection", "default", F)
pcoa.centroids(bdist14.wd$IL, meta, "Infection", "default", F)


## analysis of 16S rRNA copy distribution per ng of sample DNA


jpeg("coreplots/16SCopies.jpeg", width=5200, height=1800, res=600) # open the .jpeg device
ggplot(droplevels(meta[meta$Infection!="NEG"&meta$DPI!="13DPC",]), aes(Infection, ratio_16s, fill=Infection)) + # print the graph to the device - make sure your variable names are correct!
  facet_grid(DPI ~ BodySite) + geom_boxplot() + ggthemes::theme_few() + ylab("Copies/ng of DNA") + xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.spacing.x = unit(0.8, "lines"),
                         axis.text.y = element_text(angle=90, size=7, hjust=0.25, vjust=0.5, face="bold")) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))
dev.off() # close the device

sink("coreplots/16SCopies.txt") # create and open a .txt file
pairwise.wilcox.test(x=meta[meta$DPI=="14DPC" & meta$Infection!="NEG", "ratio"], # print this output to the text file
                     g=meta[meta$DPI=="14DPC" & meta$Infection!="NEG", "BodySiteInf"],
                     p.adjust.method = "holm")
pairwise.wilcox.test(x=meta[meta$DPI=="5DPC" & meta$Infection!="NEG", "ratio"], # print this output to the text file
                     g=meta[meta$DPI=="5DPC" & meta$Infection!="NEG", "BodySiteInf"],
                     p.adjust.method = "holm")
sink() # close the file


## quantify changes in abundance for predominant genera

dysb.5 <- facet.updown(comm5[,asvs], tax, meta, "BodySite", "Infection", "Mock", # draw plot and output statistics for 5dpi data
                     filename="5dpi_from5dpiMOCK", level="genus", w=2180, h=1500, set.lim = 10) %>%
  group_by(group, Infection) %>% summarize_at(vars(Log2.Abundance), funs(mean)) %>% .[-which(.[,"Infection"]=="Mock"),] # this small object shows the mean increase/decrease for affected genera

dysb.14 <- facet.updown(comm14[,asvs], tax, meta, "BodySite", "Infection", "Mock", # draw plot and output statistics for 14dpi data
                     filename="14dpi_from5dpiMOCK", level="genus", w=2180, h=1500, set.lim = 10) %>%
  group_by(group, Infection) %>% summarize_at(vars(Log2.Abundance), funs(mean)) %>% .[-which(.[,"Infection"]=="Mock"),] # this small object shows the mean increase/decrease for affected genera

ra.comm <- apply(comm5, 1, ra)*100
facet.updown(t(ra.comm)[,asvs], tax, meta, "BodySite", "Infection", "Mock", # draw plot and output statistics for 5dpi data
             filename="5dpi_from5dpiMOCK_ra", level="genus", w=2180, h=1500)
ra.comm <- apply(comm14, 1, ra)*100
facet.updown(t(ra.comm)[,asvs], tax, meta, "BodySite", "Infection", "Mock", # draw plot and output statistics for 5dpi data
             filename="14dpi_from5dpiMOCK_ra", level="genus", w=2180, h=1500)

## visualize correlation of predominant genera to gene expression

cyts <- c("UT_IFN.a", "UT_IFN.b","UT_IFN.g","UT_IFN.l","UT_OAS","UT_Mx","UT_IL.6","UT_LITAF") # names of metadata columns containing gene expression data
ct.cyts <- c("CT_IFN.g","CT_IFN.l", "CT_OAS", "CT_Mx","CT_IL4","CT_IL5","CT_IL6","CT_IL10","CT_IL12","CT_LITAF") # names of metadata columns containing gene expression data

for (c in cyts){ # go through cyts iteratively; they will save to separate folders
  cyt.regression(comm5[row.names(comm5) %in% meta$SampleID[meta$BodySite == "TRA"],], # plots show scatterplots and regressions with equations 
                 tax, meta[meta$DPI == "5DPC" & meta$BodySite == "TRA" & meta$Infection != "NEG",], select=gens, c, "BodySite", "Infection")
}

mar <- robust.threshold(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite!="LRT"],], meta, "BodySite", "Infection", asvs) # check robustness of margin choice

common <- core.id(comm5[,asvs], meta, tax, "BodySite", "Infection", 0.01) # presence-absence distribution of predominant genera

concat <- list("CEC"=list("Mock"=unique(unlist(common$CEC))), "ILE"=list("Mock"=unique(unlist(common$ILE))), # have to trick the function into using the same ASV list for all treatments
               "NAS"=list("Mock"=unique(unlist(common$NAS))), "TRA"=list("Mock"=unique(unlist(common$TRA))))

mb <- core.stack(comm5, concat, tax, meta, "BodySite", "Infection", "genus", # this data table describes the consensus distribution of predominant genera present in each site
                 0.005, 0.01, "_predominantgenera", core.only=T, fixed=T, landscape=T, view.all=T, facet=T)

regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite=="TRA" & meta$Infection!="NEG"], ], # this produces the correlation heatmap summarizing the regressions discussed above
                    tax, meta[meta$BodySite=="TRA" & meta$DPI=="5DPC",], as.character(unique(mb$otu[mb$groupvar=="TRA"])), 
                    cyts, "BodySite", 10, "UT_cytokines_TRAonly", "cor", 1200, 1050, F, 1, NULL)
regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite=="CEC" & meta$Infection!="NEG"], ], # this produces the correlation heatmap summarizing the regressions discussed above
             tax, meta[meta$BodySite=="CEC" & meta$DPI=="5DPC",], as.character(unique(mb$otu[mb$groupvar=="CEC"])), 
             ct.cyts, "BodySite", 10, "CT_cytokines_CEConly", "cor", 1200, 1050, F, 1, NULL)

write.csv(read.csv("coreplots/cytokines/heatmapdata_UT_cytokines_TRAonly.csv") %>% 
            reshape2::dcast(genus ~ factor, value.var = "value"), "../TableS1.csv")

## visualize correlation of predominant genera to viral titer

sites <- c("NAS", "TRA", "CEC", "ILE") # names of sites with viral titer data attached; LRT is exempt (see Materials & Methods)

for (s in sites){ # these plots will all go in the same folder
  cyt.regression(comm5[row.names(comm5) %in% meta$SampleID[meta$BodySite == s],], # plots show scatterplots and regressions with equations  
                 tax, meta[meta$DPI == "5DPC" & meta$BodySite == s & meta$Infection != "NEG",], select=gens, "EID50.mL", "BodySite", "Infection")
}

regr.heatmap(comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite!="LRT"], asvs], tax[tax$tag%in%asvs,], meta[meta$BodySite!="LRT" & meta$DPI=="5DPC",], 
             gens, "EID50.mL", "BodySite", 10, "ALL_titers", "cor", 975, 1450, T, 1, NULL)  # this produces the correlation heatmap summarizing the regressions discussed above
regr.heatmap(comm5[, asvs], tax[tax$tag%in%asvs,], meta[meta$DPI=="5DPC"&meta$Infection!="NEG",], 
             gens, "BodyWt_g", "BodySite", NULL, "ALL_weights", "cor", 975, 1450, T, 1, NULL)  # this produces the correlation heatmap summarizing the regressions discussed above

write.csv(read.csv("coreplots/cytokines/heatmapdata_ALL_weights.csv") %>% reshape2::dcast(genus ~ group, value.var = "value"),
          "coreplots/correlation_weights.csv")
## produce NMDS of common predominant genera per body site and k-means clusters

ls <- list() # initialize list to store NMDS and k-means output
for (i in sites){ # iterate across body sites
  comm.rmv <- comm5[row.names(comm5)%in%meta$SampleID[meta$BodySite==i], asvs] %>% # community table made up of ASVs ascribed to predominant genera
    t() %>% data.frame(tag=row.names(.)) %>% right_join(tax[,c("tag", "genus")], .) %>% .[,-1] %>% # match genus names to ASV tags
    group_by(genus) %>% summarize_all(funs(sum)) %>% as.data.frame() %>% # add read counts of that genus together
    .[which(.[,"genus"]%in%as.character(unique(mb$otu[mb$groupvar==i]))),] %>% `row.names<-` (.[,"genus"]) %>% .[,-1] # subset genera down to those predominant in the current body site 
  k <- kmeans(vegan::metaMDS(t(comm.rmv),"bray", trymax = 1200, autotransform = F)[[35]], 3, 100, 20) # perform NMDS and get k-means clusters out of MDS1 and MDS2
  m <- cbind(vegan::metaMDS(t(comm.rmv),"bray", trymax = 1200, autotransform = F)[[35]], k$cluster) %>% data.frame(genus=row.names(.)) # bind NMDS axis data to k-means clusters
  ls[[i]] <- m # save to list
}

titers <- read.csv("coreplots/cytokines/heatmapdata_ALL_titers.csv") # this data table was output by the heatmap function

titers$genus <- gsub("\\.", "-", titers$genus) # must adjust the characters so genus names match between tables

kmeans <- left_join(titers, rbind(data.frame(ls$NAS, group="NAS"), data.frame(ls$TRA, group="TRA"), # join the viral titer correlation data to k-means cluster 
                                  data.frame(ls$CEC, group="CEC"), data.frame(ls$ILE, group="ILE"))) %>% .[!is.na(.[,"V3"]),] # genera not common in a particular site are removed
  
k1 <- kmeans %>% group_by(group, V3) %>% summarize_at(vars(value), funs(mean, length)) # summarize correlation data for each k-means group

k1 <- k1[order(k1$group, k1$mean),] %>% data.frame(reorder=rep(1:3,4)) # reorder group designations per body site based on mean correlation value for that group

kmeans <- left_join(kmeans, k1[,c(1,2,4,5)]) # join new group designations to the kmeans table

kmeans$reorder <- factor(kmeans$reorder, labels=c("C", "B", "A")) # rename the grouping factor

kmeans <- mutate(kmeans, GroupGen=paste(reorder, genus, sep="_")) %>% .[order(.$group, .$reorder, .$genus),] %>%
  data.frame(Y=c(1:dim(kmeans)[1])) %>% group_by(group, reorder) %>% mutate(., mean=mean(Y)) # this is for reshaping the graph to be more aesthetic

titers <- reshape2::dcast(kmeans, genus ~ group, value.var = "value") # reshape the table into wide form

write.csv(titers, "coreplots/correlation viraltiters.csv") # save as output table

save(kmeans, file="kmeans_clusters.RD") # save for later

jpeg("coreplots/genus-titer_correlationchart.jpeg", width=2000, height=2000, res=300) # open the graphics device
ggplot(kmeans, aes(Y, value, fill=reorder, group=genus, label=genus)) + facet_wrap(vars(group), nrow=2, dir="v", scales="free_y") + # print graph to device
  geom_col(position=position_dodge(width=0.95), color="white", width=0.95) +
  geom_text(aes(label=signif, hjust=scales::rescale(value/abs(value))), size=8, nudge_x = -0.3) +
  geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
  geom_hline(yintercept=-0.5, color="grey80", linetype="dotted") +
  geom_hline(yintercept=0.5, color="grey80", linetype="dotted") +
  scale_x_continuous(breaks=kmeans$mean, labels=kmeans$reorder) +
  geom_text(aes(y=0, hjust=scales::rescale(value/abs(value))), 
            vjust=0.25, position=position_dodge(width=0.95), size=3) +
  ylab("Pearson's correlation coefficient\n(relative abundance vs. viral titer)") + xlab("k-means group\n(NMDS-derived)") +
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
