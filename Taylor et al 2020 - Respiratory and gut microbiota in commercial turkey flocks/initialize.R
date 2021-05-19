source("functions.R")

("./feature-table.biom")

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

meta <- read.csv("./metadata.csv") # replace this path with your own
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
comm <- comm[match(meta$SampleID,row.names(comm)),]
save(comm, file="feature_table.RD")

bdistUW <- bdistALL$unifracs[,,2]
save(bdistUW, file="./GenUniFrac_unweighted.RD")

ape::write.tree(drop.tree, "rep_set_reduced.tre")