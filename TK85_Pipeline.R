biom.as.csv("./feature-table.biom")

comm <- as.data.frame(read.csv("feature_table.csv", row.names = 1))
taxonomy1 <- read.csv("taxonomy.csv")
taxonomy2 <- read.csv("taxonomy.tsv", sep="\t")
taxonomy <- left_join(taxonomy1, taxonomy2, by=c("den.otu"="Feature.ID"))
taxonomy <- taxonomy[,c(1,3)]
names(taxonomy)[2] <- "taxonomy"
rm(taxonomy1, taxonomy2)
taxonomy$taxonomy <- as.character(taxonomy$taxonomy)

pull <- c(grep("Chloroplast", taxonomy$taxonomy, ignore.case = T), grep("Mitochondria", taxonomy$taxonomy, ignore.case = T))
taxonomy <- taxonomy[-pull,]
tax <- bact.tax(taxonomy, "silva")
save(tax, file="taxonomy.RD")

meta <- read.csv("./metadata.txt", sep="\t") # replace this path with your own
meta <- meta[which(meta$SampleID%in%names(comm)==T),]
meta <- meta[order(meta$Age, decreasing=F),]
save(meta, file="./metadata.RD")

tree <- ape::read.tree("./tree.nwk") # replace this path with your own
#if tree is unrooted, use phytools::midpoint.root(tree) to root: tree <- phytools::midpoint.root(tree)

prune <- match(tax$den.otu, tree$tip.label)
tips <- tree$tip.label[-prune]

drop.tree <- ape::drop.tip(tree, tip=tips)
comm <- comm[-pull,]
bdistALL <- GUniFrac::GUniFrac(t(comm), drop.tree, 0)

row.names(comm) <- tax$tag
comm <- as.data.frame(t(comm))
#load("./metadata.RD")
comm <- comm[match(meta$SampleID,row.names(comm)),]
save(comm, file="feature_table.RD")
#comm <- GUniFrac::Rarefy(comm, 5000)
#comm <- as.data.frame(comm$otu.tab.rff)
#save(comm, file="feature_table_rare5000.RD")

bdistUW <- bdistALL$unifracs[,,2]
save(bdistUW, file="./GenUniFrac_unweighted.RD")

ape::write.tree(drop.tree, "rep_set_reduced.tre")

## NOTE: You may avoid doing these steps again by using 
load("./feature_table.RD")
load("./metadata.RD")
load("./taxonomy.RD")
load("./GenUniFrac_unweighted.RD")

## divide distance matrix by body site
TW <- bdistUW[grep("T", row.names(bdistUW)), grep("T", row.names(bdistUW))]
CE <- bdistUW[grep("C", row.names(bdistUW)), grep("C", row.names(bdistUW))]
SW <- bdistUW[grep("S", row.names(bdistUW)), grep("S", row.names(bdistUW))]
IL <- bdistUW[grep("I", row.names(bdistUW)), grep("I", row.names(bdistUW))]

source("../functions/pcoa_centroids.R")

## pcoa plots with centroids and connections
pcoa.centroids(TW, meta, "Flock", "paired", T, "Barn") ## looking at how the flocks group 
pcoa.centroids(CE, meta, "Flock", "paired", T, "Barn")
pcoa.centroids(SW, meta, "Flock", "paired", T, "Barn")
pcoa.centroids(IL, meta, "Flock", "paired", T, "Barn")

## PERMANOVA for differences in flock group with body site
pairwiseAdonis::pairwise.adonis(bdistUW, meta$BodySiteFlock) 

pcoa.centroids(TW, meta, "Age", "rainbow", T, NULL) ## looking at how composition shifts with age
pcoa.centroids(CE, meta, "Age", "rainbow", T, NULL)
pcoa.centroids(SW, meta, "Age", "rainbow", T, NULL)
pcoa.centroids(IL, meta, "Age", "rainbow", T, NULL)

## PERMANOVA for differences in flock group with body site
vegan::adonis(TW ~ meta$AgeBodySite[meta$BodySite=="TRACHEA"]) 
vegan::adonis(CE ~ meta$AgeBodySite[meta$BodySite=="CECUM"]) 
vegan::adonis(SW ~ meta$AgeBodySite[meta$BodySite=="NASAL"]) 
vegan::adonis(IL ~ meta$AgeBodySite[meta$BodySite=="ILEUM"]) 


pcoa.plotting(bdistUW, meta, "BodySite", "Set1", "UW Unifrac") ## base pcoa with PERMANOVA output comparing body sites
pcoa.plotting(bdistUW, meta, "Age", "rainbow", "UW Unifrac") ## base pcoa with PERMANOVA output comparing ages
pcoa.centroids(bdistUW, meta, "BodySite", "Set1", F, NULL) ## providing centroids for body sites
pcoa.centroids(bdistUW, meta, "BodySiteFlock2", "paired", T, "BodySite") ## looking at proximity between flocks for each body site
pcoa.centroids(bdistUW, meta, "Age", "rainbow", T, NULL) ## looking at how composition shifts with age, regardless of body site


source("TK85_Functions.R")

## alpha diversity metrics for each body site
alpha.boxplots(comm, meta, "T", "FlockAge", FALSE)
alpha.boxplots(comm, meta, "S", "FlockAge", FALSE)
alpha.boxplots(comm, meta, "C", "FlockAge", FALSE)
alpha.boxplots(comm, meta, "I", "FlockAge", FALSE)

rd2 <- comm[row.names(comm)%in%meta$SampleID[meta$Flock2=="RD2"],]
rd3 <- comm[row.names(comm)%in%meta$SampleID[meta$Flock2=="RD3"],]

## linear trends between alpha diversity and age for each flock
alpha.boxplots(rd2, meta, "T", "Age_weeks", legend = T, type = "point", fixed=list(c(0,4),c(0,200),c(0,.81)))
alpha.boxplots(rd2, meta, "S", "Age_weeks", legend = T, type = "point", fixed=list(c(1,4.5),c(0,225),c(0.25,0.9)))
alpha.boxplots(rd2, meta, "C", "Age_weeks", legend = T, type = "point", fixed=list(c(1,5),c(0,450),c(0.3,0.9)))
alpha.boxplots(rd2, meta, "I", "Age_weeks", legend = T, type = "point", fixed=list(c(0,5.25),c(0,375),c(0.1,.91)))

## alpha diversity boxplots for each flock by age with statistcal comparisons between age groups
alpha.boxplots(rd3, meta, "T", "Age_weeks", legend = T, type = "point", fixed=list(c(0,4),c(0,200),c(0,.81)))
alpha.boxplots(rd3, meta, "S", "Age_weeks", legend = T, type = "point", fixed=list(c(1,4.5),c(0,225),c(0.25,0.9)))
alpha.boxplots(rd3, meta, "C", "Age_weeks", legend = T, type = "point", fixed=list(c(1,5),c(0,450),c(0.3,0.9)))
alpha.boxplots(rd3, meta, "I", "Age_weeks", legend = T, type = "point", fixed=list(c(0,5.25),c(0,375),c(0.1,.91)))

## predominant taxa for each body site and flock
common <- core.id(comm, meta, tax, "BodySiteFlock2", "Age", 0.5)
mb <- core.stack(comm, common, tax, meta, "BodySiteFlock2", "Age", "tag.name", 0.035, 0.5, "totalCoreProgression", core.only=F, fixed=F, landscape=F, view.all=F)
mb <- core.stack(comm, common, tax, meta, "BodySiteFlock2", "Age", "tag.name", 0.035, 0.5, "totalStackedCores", core.only=F, fixed=F, landscape=F, view.all=T)

## 
venn <- mbiom.venn(comm, meta, "Flock", tax, file="byFlocks", tax.grp="tag")
shared <- shared.taxa(venn, "byFlocks", T)
mbiom.bar(shared, comm, NULL, tax, "order", meta, "Flock", "Age", "byFlocks", T, F)

venn <- mbiom.venn(comm, meta, "BodySite", tax, file="byBodySite", tax.grp="tag")
shared <- shared.taxa(venn, "byBodySite", T)
mbiom.bar(shared, comm, NULL, tax, "order", meta, "BodySite", "Age", "byBodySiteAgeFixed", T, F, set.y=T)
mbiom.bar(shared, comm, NULL, tax, "genus", meta, "BodySite", "Age", "byBodySiteAgegenus", T, T)

lacto <- tax$tag[tax$genus%in%c("Candidatus")]
venn3 <- mbiom.venn(comm[grep("I", row.names(comm)),lacto], meta, "Flock2", tax, file="CandidatusByFlock2", tax.grp="tag")
shared3 <- shared.taxa(venn3, "CandidatusByFlock2")
mbiom.bar(shared3, comm[grep("I", row.names(comm)),], lacto, tax, "taxonomy", meta, "Flock2", "Age", "CandidatusByFlock2", T, T)

halp <- xlsx::read.xlsx("coreplots/venn/byFlocks/byFlocks.xlsx", 4) %>% 
  left_join(tax[,c(3,5)], by=c("RD3.DR"="tag.name"))
step1 <- core.id(comm[,halp$tag], meta, tax, "BodySite", "Age", 0.01)
step2 <- mbiom.venn(comm[,halp$tag], meta,"BodySite", tax, "RD3-DR", "tag.name")
step3 <- shared.taxa(step2, "RD3-DR")
step4 <- mbiom.bar(step2, comm[row.names(comm)%in%meta$SampleID[meta$Flock=="RD3-DR"],], NULL, tax, 
                   'taxonomy', meta, "BodySite", "Age", "RD3-DR", legend=T, compare.otus = T) %>%
  group_by(Group, Subgroup, Taxonomy) %>% summarize_at(vars(value), funs(sum)) %>%
  group_by(Group, Taxonomy) %>% summarize_at(vars(value), funs(mean)) %>% write.csv("coreplots/venn/RD3-DR/meanTaxonomy.csv")
for (i in 1:length(step2)){
  write.xlsx(step2[[i]], file="coreplots/venn/RD3-DR/bodysites.xlsx", sheetName = names(step2)[i], append=T)
}

halp <- xlsx::read.xlsx("coreplots/venn/byFlocks/byFlocks.xlsx", 2) %>% 
  left_join(tax[,c(3,5)], by=c("RD2.JF"="tag.name"))
step1 <- core.id(comm[,halp$tag], meta, tax, "BodySite", "Age", 0.01)
step2 <- mbiom.venn(comm[,halp$tag], meta,"BodySite", tax, "RD2.JF", "tag")
step3 <- shared.taxa(step2, "RD2.JF")
step4 <- mbiom.bar(step2, comm[row.names(comm)%in%meta$SampleID[meta$Flock=="RD2-JF"],], NULL, tax, 
                   'taxonomy', meta, "BodySite", "Age", "RD2.JF", legend=T, compare.otus = T) %>%
  group_by(Group, Subgroup, Taxonomy) %>% summarize_at(vars(value), funs(sum)) %>%
  group_by(Group, Taxonomy) %>% summarize_at(vars(value), funs(mean)) %>% write.csv("coreplots/venn/RD2.JF/meanTaxonomy.csv")
for (i in 1:length(step2)){
  write.xlsx(step2[[i]], "coreplots/venn/RD2.JF/bodysites.xlsx", sheetName = names(step2)[i], append=T)
}

genus <- read.csv("coreplots/venn/byBodySiteAgegenus/dataframe.csv") %>% 
  group_by(Group, Subgroup, Taxonomy) %>% summarize_at(vars(value), funs(sum)) %>% 
  .[.["Taxonomy"]!="",] %>% group_by(Group, Subgroup) %>% 
  summarize(max=max(value), tax=Taxonomy[which(value==max(value))]) %>% 
  write.csv("coreplots/venn/byBodySiteAgegenus/maxgenusperbodysite.csv")

firm <- which(tax$phylum=="Firmicutes")
bact <- which(tax$phylum=="Bacteroidetes")
act <- which(tax$phylum=="Actinobacteria")
prot <- which(tax$phylum=="Proteobacteria")
other <- which(!tax$phylum%in%c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria"))

comm_f <- rowSums(comm[,firm])
comm_b <- rowSums(comm[,bact])
comm_a <- rowSums(comm[,act])
comm_p <- rowSums(comm[,prot])
comm_o <- rowSums(comm[,other])


f.b <- data.frame(sample=row.names(comm), 
                  firmicutes=rowSums(comm[,firm]), 
                  bacteroidetes=rowSums(comm[,bact]), 
                  actinobacteria=rowSums(comm[,act]), 
                  proteobacteria=rowSums(comm[,prot]), 
                  other=rowSums(comm[,other])) %>%
  dplyr::left_join(meta[,c(1,6,8,10,12)], by=c("sample"="SampleID")) %>%
  melt()

f.b2 <- f.b[-grep("^B", f.b$sample),-c(1)] %>%
  group_by(BodySite, AgeBodySite,  Age, Barn, variable) %>%
  summarize(mean=mean(value)) %>%
  mutate(BodySiteFlockAge=paste(BodySite,  sep="_"))

f.b3 <- plyr::ddply(f.b2, "AgeBodySite", transform, percent_abund = mean / sum(mean) * 100)
f.b3 <- f.b3[order(f.b3$mean, decreasing=T),]
f.b3$variable <- factor(f.b3$variable, levels=rev(c("firmicutes", "bacteroidetes", "actinobacteria", "proteobacteria", "other")))

d <- ggplot(f.b3[f.b3$BodySite!="TT",], aes(Age, percent_abund, fill=variable)) +
  facet_grid(BodySite~Barn, scales="free_x")+
  geom_bar(stat="identity", width=.96) +
  theme_minimal() +
  scale_fill_manual(values=rev(c("#330033", "#0072B2", "#FFCC00", "#56B4E9","#00FF00")))+
  theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
        legend.position="right",
        legend.justification="right",
        panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey90", linetype="36"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 25, face="bold"),
        plot.title = element_text(size = 25, face="bold"),
        axis.title = element_text(size = 20, face="bold"),
        axis.text.x = element_text(color = "black", size = 20, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(color = "black", size = 15, face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 25)) +
  guides(fill=guide_legend(keywidth=1, keyheight=1, ncol=1)) +
  coord_cartesian(ylim=c(5,96)) +
  ylab("Relative Abundance (%)")
tiff(file="coreplots/phylumprofiles.tiff", width=700, height=700, pointsize=12)
d
dev.off()



wt <- read.csv("weights2.csv")
wt$Flock <- factor(wt$Flock, levels=c("Hybrid", "RD2", "RD3"))
wt$Vac_grp <- factor(wt$Vac_grp, levels=c("RD2", "RD3"))
wt$Vac <- factor(wt$Vac, levels=c("Albendazole","CTC","Di-Methox","HE Vac","M-Ninevax","Megan-Egg" ))
g <- ggplot(wt, aes(Age, Wt, color=Flock, linetype=Flock), size=3) +
  geom_point(color="white")+
  ggalt::geom_xspline(size=1.5)+
  scale_color_manual(values=colors[c(8,9,11)])+
  scale_linetype_manual(values=c("31","solid","solid")) +
  #geom_point(aes(x=Vac_age, y=y, fill=Vac_grp), shape=23, size=5, color="black", show.legend = c(shape=F, fill=T, color=F)) +
  scale_x_continuous(breaks=c(1,3,5,8,12,16), labels=c(1,3,5,8,12,16), minor_breaks = seq(0,16,1)) +
  #scale_fill_manual(values=c("tomato", "limegreen", "deepskyblue"), na.translate=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA, color="black", size=1),
        axis.text.y=element_text(size=17, face="bold"),
        axis.title.y = element_text(face="bold", size=20),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.x=element_line(linetype="dashed"),
        legend.text = element_text(size=14),
        legend.position = c(0.1, 0.6),
        legend.title=element_blank(),
        legend.background = element_rect(fill="white", color="black")) +
  xlab("Age (weeks)")+
  ylab("Weight (g)") +
  coord_cartesian(xlim=c(1,16), ylim=c(0,17500))

h <- ggplot(wt, aes(Age, Wt, color=Flock, linetype=Flock), size=3) +
  geom_point(aes(x=Vac_age, y=y, fill=Vac_grp), shape=23, size=5, color="black", show.legend = c(shape=F, fill=T, color=F)) +
  scale_x_continuous(breaks=c(1,3,5,8,12,16), labels=c(1,3,5,8,12,16), minor_breaks = seq(0,16,1)) +
  scale_fill_manual(values=c("limegreen", "deepskyblue"), na.translate=F) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA, color="black", size=1),
        plot.margin = unit(c(2,1,1,0), "lines"),
        axis.text.x=element_text(size=17, face="bold"),
        axis.title.x = element_text(face="bold", size=20),
        axis.text.y=element_blank(),
        axis.title.y=element_text(face="bold", size=15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.length.x.top = unit(1,"lines"),
        panel.grid.minor.x=element_line(linetype="dashed"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.text = element_text(size=14),
        legend.position = c(0.08,0.5),
        legend.title=element_blank(),
        legend.background = element_rect(fill="white", color="black")) +
  xlab("Age (weeks)")+
  ylab("Vaccination \n schedule") +
  coord_cartesian(xlim=c(1,16), ylim=c(-500, 1750))


tiff(file="weights2.tiff", width=3500, height=1350, res=300)
egg::ggarrange(g,h,heights=c(3,1))
dev.off()


vegan::adonis2(bdistUW ~ BodySite + Age + Flock, data=meta)

sink("weight comparison.txt")
pairwise.t.test(antibody$Weight_g[antibody$Tissue=="CE "], antibody$AgeTissueFlock[antibody$Tissue=="CE "], p.adjust.method = "bonferroni")
sink()



asvs <- read.csv("NicheSpace.csv")
d <- ggplot(asvs, aes(Distribution,NicheSpace*100, group=as.factor(Age),fill=Distribution)) +
  geom_bar(stat="identity", width=0.9, position=position_dodge2(padding=0.2))+
  geom_text(aes(x=Distribution, y=-5, group=as.factor(Age), label=Age), position=position_dodge2(width=0.9, padding=0.3), size=4)+
  geom_vline(aes(xintercept=1.5)) +
  facet_wrap(vars(BodySite), ncol=4, dir="h", scales="free_x") +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA, color="black", size=1),
        plot.margin = unit(c(2,1,1,0), "lines"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y=element_text(face="bold", size=15),
        axis.ticks.x = element_blank(),#element_line(color="black"),
        axis.ticks.length.x.top = unit(1,"lines"),
        strip.text = element_text(size = 25, face="bold"),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_line(linetype="dashed"),
        panel.grid.major.x=element_blank(),
        legend.text = element_text(size=14, margin=margin(0,40,0,0)),
        legend.position = "top",
        legend.justification = "center",
        legend.title=element_blank(),
        legend.background = element_blank())+
  coord_cartesian(ylim=c(-5,105), xlim=c(1.1,1.9)) +
  xlab("Age (in weeks)") +
  ylab("Total Niche Space (%)")


asvs2 <- xlsx::read.xlsx("coreplots/venn/byBodySiteAge/byBodySiteAge.xlsx", 15) 
asvs2 <- tax$tag[match(asvs2$all, tax$tag.name)]
asvs2 <- comm[,asvs2] %>% rowSums() %>%
  data.frame(SampleID=names(.), sumASV=.) %>%
  left_join(meta[,c(1,3,6)]) 
asvs2$sumASV <- asvs2$sumASV/50
asvs2$BodySite <- factor(asvs2$BodySite, levels=c("NASAL", "TRACHEA", "CECUM", "ILEUM"), labels=c("URT"))
asvs3 <- asvs2 %>% group_by(BodySite) %>% summarize_at(vars(sumASV), funs(sd, mean)) 


asvs4 <- asvs2
asvs4$BodySite <- factor(asvs4$BodySite, levels=c("NASAL", "TRACHEA", "CECUM", "ILEUM"), labels=c("URT", "URT", "LIT", "LIT"))
asvs6 <- asvs4 %>% group_by(BodySite, Age_weeks) %>% summarize_at(vars(sumASV), funs(sd, mean)) 
asvs5 <- asvs4 %>% group_by(BodySite) %>% summarize_at(vars(sumASV), funs(sd, mean)) 


f <- ggplot(asvs6, aes(BodySite,mean, group=as.factor(Age_weeks),fill=BodySite)) +
  geom_bar(stat="identity", width=0.9, position=position_dodge2(padding=0.2))+
  geom_text(aes(x=BodySite, y=-5, group=as.factor(Age_weeks), label=Age_weeks), position=position_dodge2(width=0.9, padding=0.3), size=4)+
  geom_vline(aes(xintercept=1.5)) +
  scale_fill_manual(values=c("turquoise", "goldenrod"))+
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA, color="black", size=1),
        plot.margin = unit(c(2,1,1,0), "lines"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y=element_text(face="bold", size=15),
        axis.ticks.x = element_blank(),#element_line(color="black"),
        axis.ticks.length.x.top = unit(1,"lines"),
        strip.text = element_text(size = 25, face="bold"),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_line(linetype="dashed"),
        panel.grid.major.x=element_blank(),
        plot.title = element_text(face="bold", size=20),
        legend.text = element_text(size=14, margin=margin(0,20,0,0)),
        legend.direction="horizontal",
        legend.position = c(0.25, 0.75),
        legend.justification = "center",
        legend.title=element_blank(),
        legend.background = element_rect(fill = "#FFFFFFE6", color="black"))+
  coord_cartesian(ylim=c(-5,105), xlim=c(1.1,1.9)) +
  xlab("Age (in weeks)") +
  ylab("Total Niche Space (%)") +
  labs(title="Total niche space (%) of ASVs observed in all body sites")

g <- ggplot() +
  geom_text(aes(x=0, y=0.023, label="SD:", hjust=0), nudge_x=1, size=5, color="black")+
  geom_text(aes(x=0, y=0.027, label="Mean:", hjust=0), nudge_x=1, size=5, color="black")+
  geom_text(data=asvs5, aes(x=mean, y=0.027, label=c("31.4%", "79.6%"), color=BodySite, hjust=1), nudge_x=-1, size=5, show.legend=F) +
  geom_text(data=asvs5, aes(x=mean, y=0.023, label=signif(sd,3), color=BodySite), nudge_x=-1, hjust=1, size=5, show.legend = F) +
  geom_vline(data=asvs5, aes(xintercept = mean, color=BodySite), linetype="dashed", size=1, show.legend = F) +
  geom_density(data=asvs4, aes(sumASV, fill=BodySite, color=BodySite), alpha=0.3, size=1.5, show.legend = c(fill=T, color=F)) +
  scale_fill_manual(values=c("turquoise", "goldenrod")) +
  scale_color_manual(values=c("turquoise", "goldenrod")) +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA, color="black", size=1),
        plot.margin = unit(c(2.5,1,1,0), "lines"),
        axis.text.x=element_text(face="bold", size=12),
        axis.text.y=element_blank(),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y=element_text(face="bold", size=15),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.length.x.top = unit(1,"lines"),
        strip.text = element_text(size = 25, face="bold"),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_line(linetype="dashed"),
        panel.grid.major.x=element_blank(),
        legend.text = element_text(size=14, margin=margin(0,20,0,0)),
        legend.direction = "horizontal",
        legend.position = c(0.55, 0.7),
        legend.justification = "center",
        legend.title=element_blank(),
        legend.background = element_rect(fill = "#FFFFFFE6", color="black"))+
  #coord_cartesian(xlim=c(4,96), ylim=c(0.001,0.045))+
  xlab("Total niche space (%)") +
  ylab("Frequency of samples") +
  labs(title="")


tiff(file="NicheSpace.tiff", width=3500, height=2700, res=300)
gridExtra::grid.arrange(d,f,g,layout_matrix=rbind(c(1,1,1,1,1,1), c(1,1,1,1,1,1), c(2,2,2,3,3,3), c(2,2,2,3,3,3)))
dev.off()

sink("NicheSpace.txt")
pairwise.t.test(asvs$NicheSpace, interaction(asvs$BodySite, asvs$Distribution), "bonferroni", paired=T)
sink()