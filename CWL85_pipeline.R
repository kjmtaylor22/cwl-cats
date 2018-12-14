## The following example pipeline was written by Kara JM Taylor for use with 
## 16S rRNA metabarcoding data, as output in table format by Qiime 1.8-1.9.
## It was written to produce the data presented in Ngunjiri et al. (2019). 
## "Longitudinal dynamics and interaction of respiratory and gut microbiota in 
## commercial chicken layers across all stages of the farm sequence."

## To perform properly, the pipeline must be used with Qiime formatted data:
## OTUs must have an associated taxonomy, either from the 
## GreenGenes database, or from the Silva release 132 classifier for Qiime. 
## Other table formats and taxonomic classifications from other versions of 
## these databases may not function as expected. 

## Please read the README.txt file prior to use. Descriptions of the functions
## and their options are given in README.txt. The functions themselves are given in
## CWL85_functions.R. 

## These functions may need to be modified to perform correctly with other data.
## Any questions, feel free to contact Kara at taylor.1895@osu.edu

setwd("D:/Manuscripts/GitHub/TestRun/")

#Initial setup: getting all files ready for use
source("CWL85_functions.R")

biom.as.csv("ck_all_rare5000.biom") # replace this path with your own

comm <- as.data.frame(t(read.csv("feature_table.csv", row.names = 1)))
taxonomy <- read.csv("taxonomy.csv", row.names = 1)

tax <- bact.tax(taxonomy, "silva")
save(tax, file="taxonomy.RD")

meta <- read.csv("./chicken-mapping-calypso-V6_allRare.csv") # replace this path with your own
save(meta, file="./metadata.RD")

tree <- ape::read.tree("./rep_set.tre") # replace this path with your own
#if tree is unrooted, use phytools::midpoint.root(tree) to root
prune <- match(tax$den.otu, tree$tip.label)
tips <- tree$tip.label[-prune]

drop.tree <- ape::drop.tip(tree, tip=tips)

bdistALL <- GUniFrac::GUniFrac(comm, drop.tree, 0)

colnames(comm) <- tax$tag
comm <- comm[-c(1:2, 393:409),] #remove unused samples
save(comm, file="feature_table.RD")

save(bdistALL, file="./GenUniFrac_ALL.RD")

bdistUW <- bdistALL$unifracs[,,2]
save(bdistUW, file="./GenUniFrac_unweighted.RD")

ape::write.tree(drop.tree, "rep_set_reduced.tre")

## NOTE: You may avoid doing these steps again by using 
load("./feature_table.RD")
load("./metadata.RD")
load("./taxonomy.RD")
load("./GenUniFrac_ALL.RD")


#Use of files and functions for basic exploratory analysis

all <- bdistUW[-c(1:2),-c(1:2)]
pcoa.plotting(all, meta, "BodySite", "Set1", "GenUniFracUW")
pcoa.plotting(all, meta, "Age", "rainbow", "GenUniFracUW")

common <- core.id(comm, meta, tax, "BodySite", "Age", 0.75)

mb <- core.stack(comm, common, tax, meta, "BodySite", "Age", 
                 "tag.name", 0.035, 0.75, "JustCores", T, fixed=F, landscape=F)
mb <- core.stack(comm, common, tax, meta, "BodySite", "Age", 
                 "tag.name", 0.035, 0.75, "StackedCores", core.only=F, fixed=F, landscape=F, view.all=T)
mb <- core.stack(comm, common, tax, meta, "BodySite", "Age", 
                 "tag.name", 0.035, 0.75, "CoreProgression", core.only=F, fixed=F, landscape=T)

dendro.heatmap(comm, tax, meta, path="./rep_set_reduced.tre", 
               "BodySite", "Age", mb, "tag.name", "coreMAJtag")

venn <- mbiom.venn(comm, meta, "BodySite", tax, file="bodysiteOTU_ALL", tax.grp="tag.name")

shared <- shared.taxa(venn, "bodysiteOTU_ALL")

mbiom.bar(shared, comm, NULL, tax, "taxonomy", meta, "BodySite", "Age", "bodysiteOTU_ALL", F)


#To examine OTU distribution as ascribed to potentially pathogenic genera

paths <- read.csv("./PotentialPathogens_short.csv")
matched <- c()
for (i in 1:length(unique(paths$Genus))){
  pull <- tax$tag[grep(as.character(unique(paths$Genus)[i]), tax$genus)]
  matched <- c(matched, pull[which(!is.na(pull))])
}

matched <- matched[which(!duplicated(matched))]

venn2 <- mbiom.venn(comm[,matched], meta, "BodySite", tax, file="pathogens", tax.grp="tag.name")
shared2 <- shared.taxa(venn2, "pathogens")
mbiom.bar(shared2, comm, matched, tax, "genus", meta, "BodySite", "Age", "pathogens", T)


#To examine OTU distribution as ascribed to the genus Lactobacillus

lacto <- tax$tag[grep("Lactobacillus", tax$taxonomy)]

venn3 <- mbiom.venn(comm[,lacto], meta, "BodySite", tax, file="lactobacillus", tax.grp="tag.name")
shared3 <- shared.taxa(venn3, "lactobacillus")
mbiom.bar(shared3, comm, lacto, tax, "taxonomy", meta, "BodySite", "Age", "lactobacillus", T)
