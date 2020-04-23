# cwl-cats
Work for CWL lab

Each branch from 'master' is an individual project, associated with a particular publication. 
The raw sequence data for these projects is save to the Sequence Read Archive (https://www.ncbi.nlm.nih.gov/sra).

The functions used for each project are saved with the project files (so you don't need to find and download them individually). 


The functions named in the script are stored in separate files and must be loaded from wherever they are saved prior to use.
Your working directory need only contain:
	1.	a *.biom feature table (as output by QIIME2)
	2. 	taxonomy.tsv (as output by QIIME2)
	2.	a *.csv containing metadata
	3.	a phylogenetic tree (*.tre) of the representative sequences (as output by Qiime)
	4. 	the *.R file "CWL85_functions.R" in this directory, or each of the following functions individually

The following are instructions for using these functions.

1.	biom.as.csv(path)
	•	“path” is a path to a biom table output by Qiime v1.8/v1.9
	•	Required packages (installed by the function): dplyr, biomformat
	•	Output from this function:
		i.	“feature_table.csv” : community matrix where rows are sample identifiers and columns are OTU identifiers
		ii.	“taxonomy.csv” : two column matrix 
			1.	column 1 (“den.otu”) is the OTU identifier given by Qiime
			2.	column 2 (“taxonomy”) is the string defining the classification in Qiime format

2.	bact.tax(taxonomy, database=NULL)
	•	“taxonomy” is a two-column matrix output by biom.as.csv
	•	“database” is a string, either “green” or “silva”, describing the database used to classify OTUs in Qiime
	•	Required packages (installed by the function): dplyr
	•	Output from this function:
		i.	Table is output both as an R object and an *.RD file saved in working directory
		ii.	The data frame has 12 columns: 
			1.	den.otu: Qiime OTU identifier
			2.	taxonomy: Lowest taxonomy assigned by Qiime
			3.	tag: ex. "otu1”
			4.	otu.name: “Genus species” if classified to species, otherwise “taxonomy (tag)”
			5.	tag.name: “taxonomy (tag)” for all OTUs
			6.	kingdom
			7.	phylum
			8.	class
			9.	order
			10.	family
			11.	genus
			12.	species

3.	pcoa.plotting(dist, meta, group, colors="rainbow", method="", axes=c(1,2,3), fixed=NULL)
	•	“dist” is a distance matrix 
	•	“meta” is the metadata file; make sure there is a SampleID column
	•	“group” is the grouping variable you intend to use to color the points
	•	“colors” is one of “rainbow”, “blind”, “Dark2”, “Set1”, or “reds” and is used for coloring points
	•	“method” says the distance method being used, and is also a tag for the file
	•	“axes” selects the principal coordinates you want to display (select 3)
	•	“fixed” can be used if you want the x and y axes of each plot to match between plots
	•	Required packages (installed by function): dplyr, ape, ggplot2, gridExtra, RColorBrewer, pairwiseAdonis
	•	Output from this function:
		i.	A *.tiff image containing each 2-dimentional plot side by side
		ii.	A text file containing PERMANOVA results for the grouping variable used

4.	core.id(comm, meta, tax, group, subgroup, margin)
	•	“comm” is the community table, samples-by-OTUs
	•	“tax” is the 12-column taxonomy file
	•	“meta” is the metadata file; make sure there is a sample ID column
	•	“group” is the largest grouping variable you intend to used
	•	“subgroup” is a nested grouping variable
	•	“margin” is the proportion of individuals within a subgroup that should have an OTU for it to be considered “core” for that subgroup
	•	Required packages (installed by function): dplyr
	•	Output is a list of groups and subgroups, with core OTUs listed within the subgroups

5.	core.stack(data, list, tax=NULL, factors, group, subgroup, hi.tax="tag", threshold=NULL, margin, date, core.only=FALSE, fixed=FALSE, landscape=TRUE, view.all=FALSE)
	•	“data” = “comm”
	•	“list” is the output list from the core.id function
	•	“tax” is the 12-column taxonomy file
	•	“factors” = “meta”
	•	“group” is the largest grouping variable you intend to used
	•	“subgroup” is a nested grouping variable
	•	“hi.tax” can be used if you want the heatmap to be labeled by something other than “tag”
	•	“threshold” is the minimum relative abundance an OTU must be present at in at least one subgroup for it to be considered “core”
	•	“margin” is the proportion of individuals within a subgroup that should have an OTU for it to be considered “core” for that subgroup; this number is just a file label
	•	“date” is just a file label, and can actually say anything
	•	“core.only” is useful in time-dependent data by specifying whether to show only the core present from birth (TRUE), or those that emerge with time (FALSE)
	•	“fixed” can be used if you want the y-axis to extend to 100
	•	“landscape” can be used if you want the facet grid to be portrait (FALSE) or landscape (TRUE)
	•	“view.all” can be used if you want all core (initial and emergent) to be viewed on the same graph
	•	Required packages (installed by function): dplyr, ggplot2, RColorBrewer
	•	Output from this function:
		i.	A set of *.tiff images for each group 
		ii.	A facet *.tiff image containing all groups
		iii.	*.csv files containing the data used in to generate all images

6.	dendro.heatmap (comm, tax, meta, path="", group, subgroup, core.list=NULL, hi.tax="tag", file="", classification=NULL, trans="log", stop=NULL)
	•	“comm” is the community table, samples-by-OTUs
	•	“tax” is the 12-column taxonomy file
	•	“meta” is the metadata file; make sure there is a SampleID column
	•	“path” is the path to the phylogenetic tree output by Qiime
	•	“group” is the largest grouping variable you intend to use
	•	“subgroup” is a nested grouping variable
	•	“core.list” may be either a numeric or the output data frame from the core.stack function
	•	“hi.tax” can be used if you want the heatmap to be labeled by something other than “tag”
	•	“file” is a tag to add to the output filename to differentiate it from other heatmaps
	•	“classification” can be used if you want to group OTUs by a particular taxonomic level prior to doing the other steps in heatmap assembly; make sure to use one of the column names in “tax”
	•	“trans” is a transformation to use on the data, if desired; suggest using “log”, but may also use “rank”
	•	“stop” may be used if you want to output one of the component figures by itself, prior to assembly; may use “top”, “den”, or “phy” for the subsetted table, cluster dendrogram, and subsetted phylogenetic tree, respectively
	•	Required packages (installed by the function): dplyr, egg, ape, vegan, ggplot2, ggtree, ggdendro, dendextend
	•	Output from this function:
		i.	*.tiff of the assembled heatmap
		ii.	*.csv of the data used

7.	mbiom.venn(mat, meta, group, tax=NULL, file="./", tax.grp=”tag”)
	•	“mat” = “comm”
	•	“meta” is the metadata file; make sure there is a SampleID column
	•	“group” is the grouping variable you intend to use
	•	“tax” is the 12-column taxonomy file
	•	“file” is a tag to add to the output filename to differentiate it from other digrams
	•	“tax.grp” can be used to group OTUs by taxonomic level prior to segregation into Venn sections; make sure the name is a column name in your “tax” file
	•	Required packages (installed by function): dplyr, VennDiagram, RColorBrewer
	•	Output from this function:
		i.	A *.tiff file of the Venn diagram
		ii.	A text file containing the parameters of the diagram
		iii.	Another text file containing a list of the OTUs assigned to each group
		iv.	An R object (list) also containing a list of the OTUs assigned to each group
		
8.	shared.taxa(mbiom, file)
	•	“mbiom” is the output R object (list) from mbiom.venn
	•	“file” is a tag to add to the output filename to match it to its associated diagram
	•	Required packages (installed by function): xlsx
	•	Output from this function:
		i.	And R object (list) with the OTUs assigned to each section of the Venn diagram
		ii.	An Excel workbook containing the same information

9.	mbiom.bar(shared, comm, select, tax, hi.tax, meta, group, subgroup, file=””, legend=F)
	•	“shared” is the output list from shared.taxa
	•	“comm” is the community table, samples-by-OTUs
	•	“select” is an option to subset “comm” by OTU
	•	“tax” is the 12-column taxonomy file
	•	“hi.tax” can be used if you want the chart legend to be labeled by something other than “tag”
	•	“meta” is the metadata file; make sure there is a Sample ID column
	•	“group” is the largest grouping variable you intend to use
	•	“subgroup” is a nested grouping variable
	•	“file” is a tag to add to the output filename to match it to its associated set of charts
	•	“legend” can be used if you want to print the legend or not; would suggest doing it once with and one without, if you want the graphs to be approximately the same size (and you can copy and paste the legends over later)
	•	Required packages (installed by function): dplyr, ggplot2
	•	Output is a folder containing a series of charts representing the contents of each section of the Venn diagram

