## need three files output from QIIME2: taxonomy.tsv, feature-table.biom, tree.nwk
## also need your metadata in *.csv format with your sample column titled "SampleID"

source("functions.R")

meta <- read.csv("metadata.csv", stringsAsFactors = T)

biom.as.csv("./feature-table.biom")

comm <- as.data.frame(read.csv("feature_table.csv", row.names = 1))
comm <- comm[,colnames(comm)%in%meta$SampleID]


taxonomy <- left_join(read.csv("taxonomy.csv"), read.csv("taxonomy.tsv", sep="\t"), by=c("den.otu"="Feature.ID"))
taxonomy <- taxonomy[taxonomy$den.otu%in%row.names(comm),]

taxonomy <- taxonomy[,c(1,3)]
names(taxonomy)[2] <- "taxonomy"
taxonomy$taxonomy <- as.character(taxonomy$taxonomy)

tax <- bact.tax(taxonomy, T, "silva")


tree <- ape::read.tree("./tree.nwk") 
#if tree is unrooted, use phytools::midpoint.root(tree) to root: tree <- phytools::midpoint.root(tree)

prune <- match(tax$den.otu, tree$tip.label)
tips <- tree$tip.label[-prune]
drop.tree <- ape::drop.tip(tree, tip=tips)
ape::write.tree(drop.tree, "rep_set_reduced.tre")


comm <- comm[match(tax$den.otu, row.names(comm)),]

bdistALL <- GUniFrac::GUniFrac(t(comm), drop.tree, 1)
save(bdistALL, file="./GenUniFrac.RD")

bdistUW <- bdistALL$unifracs[,,2]
save(bdistUW, file="./GenUniFrac_unweighted.RD")

bdistWd <- bdistALL$unifracs[,,1]
save(bdistWd, file="./GenUniFrac_weighted.RD")


row.names(comm) <- tax$tag
comm <- as.data.frame(t(comm))

## identify and isolate likely respiratory sample contaminants

Negs <- core.id(comm[row.names(comm)%in%meta$SampleID[meta$Infection=="NEG"],], meta, tax, "BodySite", "Infection", 1) ## ASVs found in 100% of negative controls extracted with each body site
mb.n <- core.stack(comm, Negs, tax, meta, "BodySite", "Infection", "tag.name", ## see how the abundance of taxa in negative controls affects compositon of other samples
                   0.005, 1, "NegativeControls", core.only=T, fixed=T, landscape=F, view.all=F, facet=T)

rmv <- Negs$NAS$NEG[Negs$NAS$NEG%in%Negs$LRT$NEG] ## the ones that are shared are more likely to be true reagent contaminants, rather than cross-contaminants
save(rmv, file="contaminants.RD") 

comm <- comm[,which(!names(comm)%in%rmv)] ## remove contaminant ASVs from the feature table
save(comm, file="feature_table.RD")

tax <- tax[which(!tax$tag%in%rmv),] ## remove contaminant ASVs from the taxonomy table
save(tax, file="taxonomy.RD")


meta$Infection <- factor(meta$Infection, levels=c("Mock","TKMN","CKPA","NEG"))
meta$BodySite <- factor(meta$BodySite, levels=c("NAS","TRA","LRT","CEC","ILE"))
meta$DPI <- factor(meta$DPI, levels=c("5DPC","14DPC","13DPC"))
save(meta, file="./metadata.RD")
