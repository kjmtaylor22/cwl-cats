setwd("GitHub/Taylor et al 2020 - Respiratory and gut microbiota in commercial turkey flocks/")

load("./feature_table.RD")
load("./metadata.RD")
load("./taxonomy.RD")
load("./GenUniFrac_unweighted.RD")

source("functions.R")

library(dplyr)
library(ggplot2)


## beta diversity plots & statistics (Figures 2, 3, S4 & S7, Table 2)

bdist <- list() # initialize list of unweighted UniFrac distance matrices subsetted by body site
for (i in c("S", "T", "C", "I")){bdist[[i]] <- bdistUW[grep(i, row.names(bdistUW)), grep(i, row.names(bdistUW))]}

save(bdist, file="GUniFrac_byBodySite.RD") ## save this file

pcoa.centroids(bdist$S, meta, "Flock", "paired", T, "Barn") # draw PCOA and output statistical analysis for each body site matrix
pcoa.centroids(bdist$`T`, meta, "Flock", "paired", T, "Barn") ## looking at how the flocks group 
pcoa.centroids(bdist$C, meta, "Flock", "paired", T, "Barn")
pcoa.centroids(bdist$I, meta, "Flock", "paired", T, "Barn")

pcoa.centroids(bdist$S, meta, "Age", "rainbow", T, NULL)  # draw PCOA and output statistical analysis for each body site matrix
pcoa.centroids(bdist$`T`, meta, "Age", "rainbow", T, NULL) ## looking at how composition shifts with age
pcoa.centroids(bdist$C, meta, "Age", "rainbow", T, NULL)
pcoa.centroids(bdist$I, meta, "Age", "rainbow", T, NULL)

pcoa.centroids(bdistUW, meta, "BodySite", "Set1", F, NULL) ## providing centroids for body sites
pcoa.centroids(bdistUW, meta, "BodySiteFlock2", "paired", T, "BodySite") ## looking at proximity between flocks for each body site
pcoa.centroids(bdistUW, meta, "Age", "rainbow", T, NULL) ## looking at how composition shifts with age, regardless of body site


## shared taxa plot (Figures 4, S2 & S5)

venn <- mbiom.venn(comm, meta, "BodySite", tax, file="byBodySite", tax.grp="tag")
shared <- shared.taxa(venn, "byBodySite", T)
mbiom.bar(shared, comm, NULL, tax, "order", meta, "BodySite", "Age", "byBodySite/AgeFixed", T, F, set.y=T)
mbiom.bar(shared, comm, NULL, tax, "genus", meta, "BodySite", "Age", "byBodySite/AgeGenus", T, T)


## alpha diversity (observed ASVs) for each body site (Table 1, Figure S6)

f1 <- comm[row.names(comm)%in%meta$SampleID[meta$Flock2=="F1"],]
f2 <- comm[row.names(comm)%in%meta$SampleID[meta$Flock2=="F2"],]

alpha <- rbind(alpha.boxplots(f1, meta, "T", "Age", F, type="point"), # bind together alpha diversity output tables for 5dpi
               alpha.boxplots(f1, meta, "S", "Age", F, type="point"),
               alpha.boxplots(f1, meta, "C", "Age", F, type="point"),
               alpha.boxplots(f1, meta, "I", "Age", F, type="point"),
               alpha.boxplots(f2, meta, "T", "Age", F, type="point"), # bind together alpha diversity output tables for 5dpi
               alpha.boxplots(f2, meta, "S", "Age", F, type="point"),
               alpha.boxplots(f2, meta, "C", "Age", F, type="point"),
               alpha.boxplots(f2, meta, "I", "Age", F, type="point")) %>%
  left_join(meta[,c("SampleID", "Flock2", "BodySite")], by=c("sample"="SampleID"))

jpeg("coreplots/alphadiversity.jpeg", width=4200, height=2100, res=350) # open the .jpeg device
gridExtra::grid.arrange(nrow=1,
                        ggplot(alpha, aes(Flock2, rich, fill=group.id)) + # print the graph to the device - make sure your variable names are correct!
                          facet_wrap(vars(BodySite), scales="free_y") + geom_boxplot(outlier.shape = 1, outlier.size = 1) + 
                          geom_vline(xintercept=1.5) + ylab("Number of ASVs per Sample") + xlab("") + ggthemes::theme_few(base_size = 14) + 
                          guides(fill=guide_legend(nrow=1)) + coord_cartesian(xlim=c(1.2,1.8)) +
                          theme(axis.ticks.x = element_blank(), 
                                panel.spacing.x = unit(0.8, "lines"), legend.position="bottom",
                                axis.text.y = element_text(angle=90, hjust=0.5, vjust=0.5, face="bold")),
                        ggplot(alpha, aes(Flock2, even, fill=group.id)) + # print the graph to the device - make sure your variable names are correct!
                          facet_wrap(vars(BodySite), scales="free_y") + geom_boxplot(outlier.shape = 1, outlier.size = 1) + 
                          geom_vline(xintercept=1.5) + ggthemes::theme_few(base_size = 14) + 
                          ylab("Uniformity of ASV Abundance (within a sample)") + xlab("") +
                          guides(fill=guide_legend(nrow=1)) + coord_cartesian(xlim=c(1.2,1.8)) +
                          theme(axis.ticks.x = element_blank(),  
                                panel.spacing.x = unit(0.8, "lines"), legend.position="bottom",
                                axis.text.y = element_text(angle=90, hjust=0.5, vjust=0.5, face="bold"))
)
dev.off() # close the device


## antibodies and weight gain trajectory plot (Figures 1 & S3D)

wt <- read.csv("Weight Trajectories.csv")

antibodies <- read.csv("Turkey Antibody Analysis.csv")[,c(2,3,5:8)] %>%
  group_by(Age, Flock) %>% summarize_all(funs(mean)) %>% reshape2::melt(id.vars=c("Age", "Flock"))

tiff(file="coreplots/flockhealthstatus.tiff", width=3500, height=2250, res=300)
egg::ggarrange(heights=c(3,1,2),
               ggplot(wt, aes(Age, Wt, color=Flock, linetype=Flock), size=3) +
                 geom_point(color="white") + ggalt::geom_xspline(size=1.5) + theme_minimal() + 
                 scale_color_manual(values=c("#0072B2","#640364", "red")) + scale_linetype_manual(values=c("solid","solid","31")) +
                 scale_x_continuous(breaks=c(1,3,5,8,12,16), labels=c(1,3,5,8,12,16), minor_breaks = seq(0,16,1)) +
                 theme(panel.border = element_rect(fill=NA, color="black", size=1),
                       panel.grid.minor.x=element_line(linetype="dashed"),
                       axis.text.y=element_text(color="black", size=17, face="bold"), axis.text.x = element_blank(),
                       axis.title.y = element_text(face="bold", size=20),
                       axis.ticks.x = element_blank(), axis.title.x = element_blank(),
                       legend.text = element_text(size=14), legend.position = c(0.1, 0.6),
                       legend.title=element_blank(), legend.background = element_rect(fill="white", color="black")) +
                 xlab("Age (weeks)") + ylab("Weight (g)") + coord_cartesian(xlim=c(1,16), ylim=c(0,17500)),
               ggplot(wt[wt$Vac_grp!="",], size=3) +
                 geom_point(aes(x=Vac_age, y=y, fill=Vac_grp), shape=23, size=5, color="black", show.legend = c(shape=F, fill=T, color=F)) +
                 scale_x_continuous(breaks=c(1,3,5,8,12,16), labels=c(1,3,5,8,12,16), minor_breaks = seq(0,16,1)) +
                 scale_fill_manual(values=c("#0072B2","#640364"), na.translate=F) + theme_minimal() + 
                 theme(panel.border = element_rect(fill=NA, color="black", size=1), plot.margin = unit(c(2,1,1,0), "lines"),
                       axis.text.x=element_text(color="black", size=17, face="bold"), axis.title.x = element_text(face="bold", size=20),
                       axis.text.y=element_blank(), axis.title.y=element_text(face="bold", size=15),
                       axis.ticks.x = element_line(color="black"), axis.ticks.length.x.top = unit(1,"lines"),
                       panel.grid.minor.x=element_line(linetype="dashed"), panel.grid.minor.y=element_blank(),
                       panel.grid.major.y=element_blank(), legend.text = element_text(size=14), legend.position = c(0.08,0.5),
                       legend.title=element_blank(), legend.background = element_rect(fill="white", color="black")) +
                 xlab("Age (weeks)") + ylab("Vaccination \n schedule") + coord_cartesian(xlim=c(1,16), ylim=c(-500, 1750)),
               ggplot(antibodies, aes(as.factor(Age), value, fill=Flock)) +
                 facet_wrap(vars(variable), ncol=4, dir="h", scales="free_y")+
                 geom_bar(stat="identity", width=.7, position="dodge") + theme_minimal() +
                 scale_x_discrete(labels=c("1","3","5","8","12","16")) +
                 scale_fill_manual(values=c("#0072B2","#640364", "red")) + 
                 theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
                       legend.position="right", legend.justification="right",
                       panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
                       panel.grid.major.y = element_line(colour="grey90", linetype="36"),
                       panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
                       strip.text = element_text(size = 25, face="bold"),
                       plot.title = element_text(size = 25, face="bold"),
                       axis.title = element_text(size = 20, face="bold"),
                       axis.text.x = element_text(color = "black", size = 20),
                       axis.text.y = element_text(color = "black", size = 15, face="bold"),
                       legend.title = element_blank(), legend.text = element_text(size = 25)) +
                 guides(fill=guide_legend(keywidth=1, keyheight=1, ncol=1)) +
                 ylab("ELISA Titer") + xlab("Age (in weeks)")
)
dev.off()


wt.fc <- meta[meta$BodySite=="CECUM",] %>%
  left_join(wt[wt$Flock=="Hybrid", 1:4], by=c("Age_weeks"="Age")) %>%
  mutate(d_wt=log2(Bird_weight.g./Wt))

jpeg("coreplots/weight_foldchangeoverHybrid_Flock.jpeg", width=1600, height=1550, res=300)
ggplot(wt.fc) + geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_boxplot(aes(Age, d_wt, fill=Flock2), show.legend = F, size=0.5) + xlab("") + 
  scale_fill_manual(values=c("#0072B2","#640364")) + ylab("Fold change over expected Hybrid weight") + 
  ggthemes::theme_few(base_size = 15) + theme(panel.background = element_rect(size=1.2, color="black")) + 
  geom_vline(xintercept=seq(1.5,5.5,1))
dev.off()

c(t.test(d_wt ~ Flock2, data=wt.fc[wt.fc$Age=="01W",]) %>% .[["p.value"]],
  t.test(d_wt ~ Flock2, data=wt.fc[wt.fc$Age=="03W",]) %>% .[["p.value"]],
  t.test(d_wt ~ Flock2, data=wt.fc[wt.fc$Age=="05W",]) %>% .[["p.value"]],
  t.test(d_wt ~ Flock2, data=wt.fc[wt.fc$Age=="08W",]) %>% .[["p.value"]],
  t.test(d_wt ~ Flock2, data=wt.fc[wt.fc$Age=="12W",]) %>% .[["p.value"]],
  t.test(d_wt ~ Flock2, data=wt.fc[wt.fc$Age=="16W",]) %>% .[["p.value"]]) %>%
  p.adjust(., "holm")



## emergent taxa plot (Figure 5)

venn <- mbiom.venn(comm, meta, "Flock", tax, file="byFlocks", tax.grp="tag.name")
shared <- shared.taxa(venn, "byFlocks", T)
mbiom.bar(shared, comm, NULL, tax, "order", meta, "Flock", "Age", "byFlocks/5A-5B", T, F)

f1 <- meta$SampleID[meta$Flock=="F1G"]
halp1 <- xlsx::read.xlsx("coreplots/venn/byFlocks/byFlocks.xlsx", 2) %>% 
  left_join(tax[,c("tag", "tag.name")], by=c("F1G"="tag.name"))
step1 <- core.id(comm[row.names(comm)%in%f1, halp1$tag], meta, tax, "BodySite", "Age", 0.01)
step2 <- mbiom.venn(comm[row.names(comm)%in%f1, halp1$tag], meta,
                    "BodySite", tax, "byFlocks/5C-F1NovelInGrowOut", "tag") %>%
  mbiom.bar(., comm[row.names(comm)%in%f1,], NULL, tax, 'taxonomy', meta, 
                   "BodySite", "Age", "byFlocks/F1NovelInGrowOut", legend=T, compare.otus = T) %>% 
  group_by(Group, tag, Taxonomy) %>% summarize_at(vars(value), funs(sum)) %>%
  write.csv(., "coreplots/venn/byFlocks/F1NovelInGrowOut/byBodySite.csv")

f2 <- meta$SampleID[meta$Flock=="F2G"]
halp2 <- xlsx::read.xlsx("coreplots/venn/byFlocks/byFlocks.xlsx", 4) %>% 
  left_join(tax[,c("tag", "tag.name")], by=c("F2G"="tag.name"))
step1 <- core.id(comm[row.names(comm)%in%f2, halp2$tag], meta, tax, "BodySite", "Age", 0.01)
step2 <- mbiom.venn(comm[row.names(comm)%in%f2, halp2$tag], meta,
                    "BodySite", tax, "byFlocks/5C-F2NovelInGrowOut", "tag") %>%
  mbiom.bar(., comm[row.names(comm)%in%f2,], NULL, tax, 'taxonomy', meta, 
                   "BodySite", "Age", "byFlocks/F2NovelInGrowOut", legend=T, compare.otus = T) %>%
  group_by(Group, tag, Taxonomy) %>% summarize_at(vars(value), funs(sum)) %>%
  write.csv(., "coreplots/venn/byFlocks/F2NovelInGrowOut/byBodySite.csv")


## phylum-level profiles (Figure S1)

phyla <- data.frame(sample=row.names(comm), 
                    firmicutes=rowSums(comm[,which(tax$phylum=="Firmicutes")]), 
                    bacteroidetes=rowSums(comm[,which(tax$phylum=="Bacteroidetes")]), 
                    actinobacteria=rowSums(comm[,which(tax$phylum=="Actinobacteria")]), 
                    proteobacteria=rowSums(comm[,which(tax$phylum=="Proteobacteria")]), 
                    other=rowSums(comm[,which(!tax$phylum%in%c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria"))])) %>%
  left_join(meta[,c("SampleID", "BodySite", "AgeBodySite",  "Age", "Barn")], by=c("sample"="SampleID")) %>% 
  reshape2::melt() %>% group_by(BodySite, AgeBodySite,  Age, Barn, variable) %>% 
  summarize(mean=mean(value)) %>% mutate(BodySiteFlockAge=paste(BodySite,  sep="_"))

phyla <- plyr::ddply(phyla, "AgeBodySite", transform, percent_abund = mean / sum(mean) * 100) %>%
  .[order(.$mean, decreasing=T),]

phyla$variable <- factor(phyla$variable, levels=rev(c("firmicutes", "bacteroidetes", "actinobacteria", "proteobacteria", "other")))

tiff(file="coreplots/phylumprofiles.tiff", width=700, height=700, pointsize=12)
ggplot(phyla, aes(Age, percent_abund, fill=variable)) + theme_minimal() +
  facet_grid(BodySite~Barn, scales="free_x") + geom_bar(stat="identity", width=.96) +
  scale_fill_manual(values=c("#666666", "#00FF00", "#FF0033", "#FF9933", "#0072B2"))+
  theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
        panel.grid.major.x = element_line(color="grey50", linetype="dotted"), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey90", linetype="36"), panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 25, face="bold"), plot.title = element_text(size = 25, face="bold"),
        legend.position="right", legend.justification="right", axis.title = element_text(size = 20, face="bold"),
        axis.text.x = element_text(color = "black", size = 20, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(color = "black", size = 15, face="bold"),
        legend.title = element_blank(), legend.text = element_text(size = 25)) +
  guides(fill=guide_legend(keywidth=1, keyheight=1, ncol=1)) +
  coord_cartesian(ylim=c(5,96)) + ylab("Relative Abundance (%)")
dev.off()


## select taxa effects (Figure S3A & S3B)

pull <- tax$tag[tax$taxonomy%in%c("Candidatus Arthromitus", "Mycoplasma", "Ornithobacterium rhinotracheale")]

select <- apply(comm, MARGIN=1, FUN=ra) %>% t() %>% .[,pull] %>%
  data.frame(SampleID=row.names(.)) %>% reshape2::melt() %>%
  left_join(meta[,c("SampleID", "BodySite", "Flock2", "Age")]) %>% 
  left_join(tax[,c("tag","taxonomy")], by=c("variable"="tag")) %>%
  group_by(taxonomy, Age, Flock2, BodySite, SampleID) %>% summarize_at(vars(value), funs(sum)) %>%
  group_by(taxonomy, Age, Flock2, BodySite) %>% summarize_at(vars(value), funs(mean))
  
jpeg(file="coreplots/selecttaxa.tiff", width=6000, height=4000, res=300)
ggplot(select, aes(Age, value*100, fill=Flock2)) +
  facet_wrap(vars(taxonomy, BodySite), ncol=4, dir="h", scales="free_y")+
  geom_bar(stat="identity", width=.7, position="dodge") + theme_minimal() +
  scale_fill_manual(values=c("#0072B2","#640364")) + 
  theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
        legend.position="right", legend.justification="right",
        panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
        panel.grid.major.y = element_line(colour="grey90", linetype="36"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 25, face="bold"),
        plot.title = element_text(size = 25, face="bold"),
        axis.title = element_text(size = 20, face="bold"),
        axis.text.x = element_text(color = "black", size = 20, angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(color = "black", size = 20),
        legend.title = element_blank(), legend.text = element_text(size = 25)) +
  guides(fill=guide_legend(keywidth=1, keyheight=1, ncol=1)) +
  ylab("Relative Abundance (%)") + xlab("Age (in weeks)")
dev.off()


## select taxa effects (Figure S3C)

Myco.ORT <- rbind(
  data.frame(SampleID=row.names(comm), M=rowSums(comm[,tax$tag[tax$genus=="Mycoplasma"]]),
             ORT=comm[,tax$tag[tax$genus=="Ornithobacterium"]]) %>%
    right_join(meta[meta$BodySite=="NASAL",]) %>% 
    left_join(wt[wt$Flock=="Hybrid", 1:4], by=c("Age_weeks"="Age")) %>%
    mutate(d_wt=Bird_weight.g./Wt, M.presence=M>0, ORT.presence=ORT>0),
  data.frame(SampleID=row.names(comm), M=rowSums(comm[,tax$tag[tax$genus=="Mycoplasma"]]),
             ORT=comm[,tax$tag[tax$genus=="Ornithobacterium"]]) %>%
    right_join(meta[meta$BodySite=="TRACHEA",]) %>% 
    left_join(wt[wt$Flock=="Hybrid", 1:4], by=c("Age_weeks"="Age")) %>%
    mutate(d_wt=Bird_weight.g./Wt, M.presence=M>0, ORT.presence=ORT>0)
)

jpeg("coreplots/FoldChange_M-ORT.jpeg", width=1200, height=3000, res=300)
egg::ggarrange(nrow=2,
ggplot(Myco.ORT) + geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_boxplot(aes(M.presence, log2(d_wt), fill=M.presence), size=0.5) + xlab("") + 
  facet_wrap(vars(BodySite), nrow=1, scales="free_y") +
  ylab("Fold change over expected Hybrid weight") + 
  ggthemes::theme_few(base_size = 15) + 
  theme(panel.background = element_rect(size=1.2, color="black"), legend.position = "bottom",
        legend.direction = "horizontal", axis.ticks.x = element_blank(), axis.text.x = element_blank()),
ggplot(Myco.ORT) + geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_boxplot(aes(ORT.presence, log2(d_wt), fill=ORT.presence), size=0.5) + xlab("") + 
  facet_wrap(vars(BodySite), nrow=1, scales="free_y") +
  ylab("Fold change over expected Hybrid weight") + 
  ggthemes::theme_few(base_size = 15) + 
  theme(panel.background = element_rect(size=1.2, color="black"), legend.position = "bottom",
        legend.direction = "horizontal", axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
)
dev.off()

sink("coreplots/FoldChange_M-ORT.txt")
  t.test(log2(d_wt) ~ M.presence, data=Myco.ORT[Myco.ORT$BodySite=="NASAL",])
  t.test(log2(d_wt) ~ ORT.presence, data=Myco.ORT[Myco.ORT$BodySite=="NASAL",])
  t.test(log2(d_wt) ~ M.presence, data=Myco.ORT[Myco.ORT$BodySite=="TRACHEA",])
  t.test(log2(d_wt) ~ ORT.presence, data=Myco.ORT[Myco.ORT$BodySite=="TRACHEA",])
sink()

# progression of predominant taxa (Figures S8-S11)

common <- core.id(comm, meta, tax, "BodySiteFlock2", "Age", 0.5)
mb1 <- core.stack(comm, common, tax, meta, "BodySiteFlock2", "Age", "tag.name", 0.035, 0.5, "CoreProgression", core.only=F, fixed=F, landscape=T, view.all=F)
mb2 <- core.stack(comm, common, tax, meta, "BodySiteFlock2", "Age", "tag.name", 0.035, 0.5, "StackedCores", core.only=F, fixed=F, landscape=F, view.all=T)

