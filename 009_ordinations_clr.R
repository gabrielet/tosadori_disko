# source("009_ordinations_clr.R")

# load libraries
library("tidyverse")
library("reshape2")
library("scales")
library("ellipse")
library("ggpubr")
library("ggsignif")
library("ggrepel")
library("ggtext")
library("gridExtra")
library("phyloseq")
library("vegan")
library("RColorBrewer")

# import functions
source("/microbiology/disko2013/code/000_micro_functions_disko2013.R")

######### SELECT EXPERIMENT AND DEFINE PATHS #########
# select exp name
exp_name <- "submission_final_norejects_dada_new_bootstrap"

# select the proper organism
orgn <- "bacteria"

# select the clustering method
if (orgn == "bacteria") {
	clust_method <- "ASV"
} else if (orgn == "fungi") {
	clust_method <- "OTU"
}

# select the taxa assigment method
taxa_algo <- "NBC"

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, "  EXPERIMENT NAME ", exp_name))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/microbiology/disko2013/"
# this is the general path to the experiment
path_to_exp <- paste0(root_path, orgn_dir, "experiments/", exp_name, "/")
# this is the path where counts results will be stored
path_to_taxa <- paste0(path_to_exp, clust_method, "/")
save_img <- paste0(path_to_taxa, "figures/", taxa_algo, "/")
save_phylo <- paste0(path_to_taxa, "phyloseq/", taxa_algo, "/")
save_orderings <- paste0(save_img, "orderings/")
save_orderings_simple <- paste0(save_orderings, "simple/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_orderings_simple), "dir exists!", dir.create(save_orderings_simple, recursive=TRUE))

# print info regarding which directory is analysed, i.e. dada or qiime
print(paste0("the analysis is performed on ", toupper(path_to_exp)))

# load centered-transformed object
centered_phylo_data <- readRDS(paste0(save_phylo, "centered_phylo_data_", exp_name, ".Rds"))

# get metadata
centered_phylo_metadata <- sample_data(centered_phylo_data)

# set levels for Treatment
centered_phylo_metadata$Treatment <- factor(centered_phylo_metadata$Treatment, levels=unique(centered_phylo_metadata$Treatment))
# and for TimePoint
centered_phylo_metadata$TimePoint <- factor(centered_phylo_metadata$TimePoint, levels=unique(centered_phylo_metadata$TimePoint))

############################################################################## PERFORMING ORDINATIONS USING PCoA ON THE FULL DATASET

print("ORDINATIONS STEP ON DNA AND RNA TOGETHER")

# compute PCoA using phyloseq function ordinate()
# using the full dataset, i.e. DNA and RNA, together
invisible(capture.output(computed_ordination_all <- ordinate(centered_phylo_data, "PCoA", "euclidean")))

# get first two vectors from the ordination
data_to_plot <- cbind.data.frame(pcoa_one = computed_ordination_all$vectors[, 1], pcoa_two = computed_ordination_all$vectors[, 2])

# compute variance pcoa
# get eigenvecs
eigenvecs <- computed_ordination_all$values$Relative_eig
# compute fraction variability
frac_var <- eigenvecs[c(1, 2)] / sum(eigenvecs)
# Percent variability, percvar
perc_var <- round(100*frac_var, 2)

# get further data
treatment <- vector()
timepoint <- vector()
dna_rna <- vector()
collectionsite <- vector()

# for each sample in data
for (ind in c(1:nrow(data_to_plot))) {
	treatment <- append(treatment, as.character(centered_phylo_metadata$Treatment[which(centered_phylo_metadata$sampleID==rownames(data_to_plot)[ind])]))
	timepoint <- append(timepoint, as.character(centered_phylo_metadata$TimePoint[which(centered_phylo_metadata$sampleID==rownames(data_to_plot)[ind])]))
	dna_rna <- append(dna_rna, as.character(centered_phylo_metadata$DNA_RNA[which(centered_phylo_metadata$sampleID==rownames(data_to_plot)[ind])]))
	collectionsite <- append(collectionsite, as.character(centered_phylo_metadata$CollectionSite[which(centered_phylo_metadata$sampleID==rownames(data_to_plot)[ind])]))
}
# put all together
data_to_plot <- cbind.data.frame(data_to_plot, treatment, timepoint, dna_rna, collectionsite)
#set levels
data_to_plot$treatment <- ifelse(data_to_plot$treatment=="C", "Control", "Warming")
data_to_plot$treatment <- factor(data_to_plot$treatment, levels=unique(data_to_plot$treatment))
data_to_plot$timepoint <- factor(data_to_plot$timepoint, levels=unique(data_to_plot$timepoint))
data_to_plot$dna_rna <- factor(data_to_plot$dna_rna, levels=unique(data_to_plot$dna_rna))

# define shapes for timepoint. both shapes must be in the same order
timeshape <- c(21, 23, 24)
names(timeshape) <- unique(data_to_plot$timepoint)

# define shapes for xna. both shapes must be in the same order
treat_shape <- c(21, 24)
names(treat_shape) <- unique(data_to_plot$treatment)

# get palette
xna_palette <- brewer.pal(n=4, name = "Set2")[c(1,4)]
names(xna_palette) <- c("DNA", "RNA")

# plot without facets
plT_xna_site <- ggplot(data_to_plot, aes(x=pcoa_one, y=pcoa_two, fill=dna_rna, shape=treatment)) +
	geom_point(size=20, color="black", stroke=3) +
#	geom_text(aes(label=collectionsite), col="black", size=10) +
	scale_fill_manual(values=xna_palette) +
	scale_shape_manual(values=treat_shape) +
#	ggtitle("PCoA analysis") +
	ggtitle("") +
	xlab(paste0("PCoA1 (", perc_var[1], "%)")) +
	ylab(paste0("PCoA2 (", perc_var[2], "%)")) +
	guides(fill=guide_legend(override.aes=list(shape=22, size=20)), shape=guide_legend(override.aes=list(size=20))) +
	ggplot_theme(leg_pos="bottom", ang_le=0)

# and export
export_svg(paste0(save_orderings_simple, "collection_site_GGPLOT_ordination_all_samples_together_treatment_DNA_RNA_", exp_name, "_centered_PCoA_CORRECTED"), plT_xna_site, as_rds=T)

# plot without facets
plT_xna_time <- ggplot(data_to_plot, aes(x=pcoa_one, y=pcoa_two, fill=dna_rna, shape=treatment)) +
	geom_point(size=20, color="black", stroke=3) +
	geom_text(aes(label=timepoint), col="black", size=10) +
	scale_fill_manual(values=xna_palette) +
	scale_shape_manual(values=treat_shape) +
#	ggtitle("PCoA analysis") +
	ggtitle("") +
	xlab(paste0("PCoA1 (", perc_var[1], "%)")) +
	ylab(paste0("PCoA2 (", perc_var[2], "%)")) +
	guides(fill=guide_legend(override.aes=list(shape=22, size=20)), shape=guide_legend(override.aes=list(size=20))) +
	ggplot_theme(leg_pos="bottom", ang_le=0)

# and export
export_svg(paste0(save_orderings_simple, "timepoint_GGPLOT_ordination_all_samples_together_treatment_DNA_RNA_", exp_name, "_centered_PCoA_CORRECTED"), plT_xna_time, as_rds=T)

# plot without facets
plT_xna <- ggplot(data_to_plot, aes(x=pcoa_one, y=pcoa_two, fill=dna_rna, shape=treatment)) +
	geom_point(size=20, color="black", stroke=3) +
#	geom_text(aes(label=collectionsite), col="black", size=10) +
	scale_fill_manual(values=xna_palette) +
	scale_shape_manual(values=treat_shape) +
#	ggtitle("PCoA analysis") +
	ggtitle("") +
	xlab(paste0("PCoA1 (", perc_var[1], "%)")) +
	ylab(paste0("PCoA2 (", perc_var[2], "%)")) +
	guides(fill=guide_legend(override.aes=list(shape=22, size=20)), shape=guide_legend(override.aes=list(size=20))) +
	ggplot_theme(leg_pos="bottom", ang_le=0)

# and export
export_svg(paste0(save_orderings_simple, "GGPLOT_ordination_all_samples_together_treatment_DNA_RNA_", exp_name, "_centered_PCoA_CORRECTED"), plT_xna, as_rds=T)

############################################################################## PERFORMING ORDINATIONS USING PCoA ON DNA OR RNA, SEPARATELY

print("ORDINATIONS STEP ON DNA AND RNA, SEPARATELY")

# define an empty variable to store the ordinations
# this might be used afterwards to evaluate the parameters
# that were computed
computed_ordination_xna <- list()

# define variable to store plots
plT_by_treat_time_site_xna <- list()

# and phylo objs
phylo_xna <- list()

# loop through DNA and RNA
for (xna in c("DNA", "RNA")) {

	# get DNA or RNA samples only
	phylo_xna[[xna]] <- prune_samples(as.character(centered_phylo_metadata$sampleID[which(centered_phylo_metadata$DNA_RNA == xna)]), centered_phylo_data)
	
	# remove ASVs which are always zero in the subset under consideration
	phylo_xna[[xna]] <- prune_taxa(apply(tax_table(phylo_xna[[xna]]), 1, function(x) length(which(x==0))!=length(x)), phylo_xna[[xna]])
	
	# compute PCoA using phyloseq function ordinate()
	# and store the results for both DNA and RNA
	invisible(capture.output(computed_ordination_xna[[xna]] <- ordinate(phylo_xna[[xna]], "PCoA", "euclidean")))

	# compute variance pcoa
	# get eigenvecs
	eigenvecs_xna <- computed_ordination_xna[[xna]]$values$Relative_eig
	# compute fraction variability
	frac_var_xna <- eigenvecs_xna[c(1, 2)] / sum(eigenvecs_xna)
	# Percent variability, percvar
	perc_var_xna <- round(100*frac_var_xna, 2)

	# get first two vectors from the ordination
	data_to_plot <- cbind.data.frame(pcoa_one = computed_ordination_xna[[xna]]$vectors[, 1], pcoa_two = computed_ordination_xna[[xna]]$vectors[, 2])

	# get further data
	treatment <- vector()
	timepoint <- vector()
	dna_rna <- vector()
	collectionsite <- vector()
	
	# get metadata
	meta_data <- sample_data(phylo_xna[[xna]])

	# for each sample in data
	for (ind in c(1:nrow(data_to_plot))) {
		treatment <- append(treatment, as.character(meta_data$Treatment[which(meta_data$sampleID==rownames(data_to_plot)[ind])]))
		timepoint <- append(timepoint, as.character(meta_data$TimePoint[which(meta_data$sampleID==rownames(data_to_plot)[ind])]))
		dna_rna <- append(dna_rna, as.character(meta_data$DNA_RNA[which(meta_data$sampleID==rownames(data_to_plot)[ind])]))
		collectionsite <- append(collectionsite, as.character(meta_data$CollectionSite[which(meta_data$sampleID==rownames(data_to_plot)[ind])]))
	}
	# put all together
	data_to_plot <- cbind.data.frame(data_to_plot, treatment, timepoint, dna_rna, collectionsite)
	#set levels
	data_to_plot$treatment <- ifelse(data_to_plot$treatment=="C", "Control", "Warming")
	data_to_plot$treatment <- factor(data_to_plot$treatment, levels=c("Control", "Warming"))
	data_to_plot$timepoint <- factor(data_to_plot$timepoint, levels=c("June", "July", "August"))

	# plot without facets and make sure to select the proper ordinations
	# and the proper nucleic acid: computed_ordination_xna[[xna]]
	plT_ttcxna <- ggplot(data_to_plot, aes(x=pcoa_one, y=pcoa_two, fill=treatment, shape=timepoint)) +
				geom_point(size=20, color="black", stroke=3) +
				scale_shape_manual(values=timeshape) +
				ggtitle("") +
				xlab(paste0("PCoA1 (", perc_var_xna[1], "%)")) +
				ylab(paste0("PCoA2 (", perc_var_xna[2], "%)")) +
				guides(fill=guide_legend(override.aes=list(shape=22, size=20)), shape=guide_legend(override.aes=list(size=20))) +
				ggplot_theme(leg_pos="bottom", ang_le=0)

	plT_by_treat_time_site_xna[[xna]] <- plT_ttcxna

	# and export
	export_svg(paste0(save_orderings_simple, "GGPLOT_", xna, "_ordination_all_samples_together_treatment_by_timepoint_", exp_name, "_centered_PCoA_CORRECTED"), plT_ttcxna, as_rds=T)
}
