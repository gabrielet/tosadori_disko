# source("013_adonis_clr.R")

# load libraries
library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("ellipse")
library("ggpubr")
library("ggsignif")
library("ggrepel")
library("ggtext")
library("gridExtra")
library("phyloseq")
library("vegan")
library("viridis")

# import functions
source("/mnt/cinqueg/gabriele/work/microbiology/disko2013/code/000_micro_functions_disko2013.R")

######### SELECT EXPERIMENT AND DEFINE PATHS #########
# select exp name
exp_name <- "submission_final_norejects_dada_new_bootstrap" ; orgn <- "bacteria"

# select the clustering method. for bacteria is ASV
if (orgn == "bacteria") {
	clust_method <- "ASV"
# for fungi is OTU
} else if (orgn == "fungi") {
	clust_method <- "OTU"
}

# select the taxa assigment method
taxa_algo <- "NBC" ; boot <- 80

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, "  EXPERIMENT NAME ", exp_name))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/mnt/cinqueg/gabriele/work/microbiology/disko2013/"
# this is the general path to the experiment
path_to_exp <- paste0(root_path, orgn_dir, "experiments/", exp_name, "/")
# this is the path where counts results will be stored
path_to_taxa <- paste0(path_to_exp, clust_method, "/")
save_img <- paste0(path_to_taxa, "figures/", taxa_algo, "/")
save_phylo <- paste0(path_to_taxa, "phyloseq/", taxa_algo, "/")
save_adonis <- paste0(save_img, "adonis/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_adonis), "dir exists!", dir.create(save_adonis, recursive=TRUE))

# print info regarding which directory is analysed, i.e. dada or qiime
print(paste0("the analysis is performed on ", toupper(path_to_exp)))

# load clr-transformed object
centered_phylo_data <- readRDS(paste0(save_phylo, "centered_phylo_data_", exp_name, ".Rds"))

# load metadata table
centered_metadata <- sample_data(centered_phylo_data)

# list environmental pararameters that might require imputation
cols_of_interest <- c("SM_gwater_g_soil", "SM_perc", "pH", "NH4_ug_g_soil", "NO3_ug_g_soil", "PO4_ug_g_soil", "perc_of_N", "perc_of_C", "mcrA", "mxaF", "nifH", "nirS", "nosZ", "TOC_mg_g_soil", "microC_mg_g_soil", "DON_mg_g_soil", "TN_mg_g_soil", "microN_mg_g_soil", "microCN", "soil_d15N", "soil_d13C", "SOM_perc", "soil_water_g", "CNratio")

# make sure the variables that will be used in the formula are
# factors since statistical models require it, to work properly
sample_data(centered_phylo_data)$Treatment <- ifelse(sample_data(centered_phylo_data)$Treatment=="C", "Control", "Warming")
sample_data(centered_phylo_data)$Treatment <- factor(sample_data(centered_phylo_data)$Treatment, levels=unique(sample_data(centered_phylo_data)$Treatment))
sample_data(centered_phylo_data)$TimePoint <- factor(sample_data(centered_phylo_data)$TimePoint, levels=unique(sample_data(centered_phylo_data)$TimePoint))
sample_data(centered_phylo_data)$CollectionSite <- factor(sample_data(centered_phylo_data)$CollectionSite, levels=unique(sample_data(centered_phylo_data)$CollectionSite))

######### ADONIS ANALYSIS TO TEST environmental parameters

# create two phylo objects with only DNA or RNA samples
phylo_xna <- list()
phylo_xna[["DNA"]] <- prune_samples(as.character(sample_data(centered_phylo_data)$sampleID[which(sample_data(centered_phylo_data)$DNA_RNA == "DNA")]), centered_phylo_data)
phylo_xna[["RNA"]] <- prune_samples(as.character(sample_data(centered_phylo_data)$sampleID[which(sample_data(centered_phylo_data)$DNA_RNA == "RNA")]), centered_phylo_data)

# define variables
adonis_test <- list()

####################################################################################### perform adonis using TimePoint as a variable

print("adonis/permanova on DNA and RNA, separately, by treatment on all timepoint")

# loop throught DNA/RNA
for (xna in names(phylo_xna)) {

#	xna <- "DNA"

	# remove ASVs which are always zero for the subset under consideration
	phylo_xna[[xna]] <- prune_taxa(apply(tax_table(phylo_xna[[xna]]), 1, function(x) length(which(x==0))!=length(x)), phylo_xna[[xna]])
#	phylo_xna[[xna]] <- prune_taxa(taxa_sums(phylo_xna[[xna]]) > 0, phylo_xna[[xna]])

	# get taxa table and transpose it in to obtain a table which is
	# sample X species as described in the Adonis documentation
	transposed_taxa <- t(otu_table(phylo_xna[[xna]]))

	# get metadata
	adonis_meta <- as.data.frame(as.matrix(sample_data(phylo_xna[[xna]])))

	# make sure the variables in the formula are factors since
	# statistical models require it, to work properly
	adonis_meta$Treatment <- factor(adonis_meta$Treatment, levels=unique(adonis_meta$Treatment))
	adonis_meta$Treatment <- relevel(adonis_meta$Treatment, ref="Control")
	adonis_meta$TimePoint <- factor(adonis_meta$TimePoint, levels=unique(adonis_meta$TimePoint))
	adonis_meta$TimePoint <- relevel(adonis_meta$TimePoint, ref="July")
	adonis_meta$CollectionSite <- factor(adonis_meta$CollectionSite, levels=unique(adonis_meta$CollectionSite))

	# set seed
	set.seed(131)

	# using adonis2 to perform the analysis

	# set number of permutations
	permuts <- how(nperm = 999)
	
	# set blocking variable, i.e. the sampling sites
	adonis_meta$CollectionSite <- as.factor(adonis_meta$CollectionSite)

	# run adonis2 on euclidean-curtis distance using CollectionSite
	# as blocking factor
	adonis_test[[xna]] <- adonis2(transposed_taxa~CollectionSite+Treatment*TimePoint, data=adonis_meta, permutations=permuts, method="euclidean")

	# run betadisper for treatment
	dispersion_treat <- betadisper(vegdist(transposed_taxa, method="euclidean"), group=adonis_meta$Treatment)
	# plot function
	box_treat <- function() { boxplot(dispersion_treat) }
	# export plot
	export_figs_tabs(paste0(save_adonis, "box_treat_", xna, "_", exp_name), box_treat(), base=TRUE, width=168*2, height=168*3)
	# plot function
	disper_treat <- function() { plot(dispersion_treat, hull=TRUE, ellipse=FALSE) }
	# export plot
	export_figs_tabs(paste0(save_adonis, "disper_treat_", xna, "_", exp_name), disper_treat(), base=TRUE, width=168*2, height=168*3)

	# run betadisper for timepoint
	dispersion_time <- betadisper(vegdist(transposed_taxa, method="euclidean"), group=adonis_meta$TimePoint)
	# plot function
	box_time <- function() { boxplot(dispersion_time) }
	# export plot
	export_figs_tabs(paste0(save_adonis, "box_time_", xna, "_", exp_name ), box_time(), base=TRUE, width=168*2, height=168*3)
	# plot function
	disper_time <- function() { plot(dispersion_time, hull=TRUE, ellipse=FALSE) }
	# export plot
	export_figs_tabs(paste0(save_adonis, "disper_time_", xna, "_", exp_name ), disper_time(), base=TRUE, width=168*2, height=168*3)

	sink(paste0(save_adonis, "adonis_with_interaction_", xna, ".txt"))

	cat("\n")
	cat(paste0("Adonis test for ", xna, "~CollectionSite+Treatment*TimePoint"))
	cat("\n")
	print(adonis_test[[xna]])

	cat("\n")
	cat(paste0("Permutest test for betadisper on ", xna, " and Treatment"))
	cat("\n")

	print(permutest(dispersion_treat))

	cat("\n")
	cat(paste0("Anova test for betadisper on ", xna, " and Treatment"))
	cat("\n")

	print(anova(dispersion_treat))

	cat("\n")
	cat(paste0("Permutest test for betadisper on ", xna, " and TimePoint"))
	cat("\n")

	print(permutest(dispersion_time))

	cat("\n")
	cat(paste0("Anova test for betadisper on ", xna, " and TimePoint"))
	cat("\n")

	print(anova(dispersion_time))

	sink()
}

####################################################################################### perform adonis DNA versus RNA

print("adonis/permanova DNA versus RNA")

# perform adonis on DNA versus RNA
transposed_taxa <- t(otu_table(centered_phylo_data))

# get metadata
adonis_meta <- as(sample_data(centered_phylo_data), "data.frame")

# make sure the variables in the formula are factors since
# statistical models require it, to work properly
adonis_meta$Treatment <- factor(adonis_meta$Treatment, levels=unique(adonis_meta$Treatment))
adonis_meta$TimePoint <- factor(adonis_meta$TimePoint, levels=unique(adonis_meta$TimePoint))
adonis_meta$CollectionSite <- factor(adonis_meta$CollectionSite, levels=unique(adonis_meta$CollectionSite))

# set reference level
adonis_meta$DNA_RNA <- relevel(as.factor(adonis_meta$DNA_RNA), ref="DNA")
adonis_meta$Treatment <- relevel(adonis_meta$Treatment, ref="Control")
adonis_meta$TimePoint <- relevel(adonis_meta$TimePoint, ref="June")

# set seed
set.seed(131)

# using adonis2 to perform the analysis

# set number of permutations
permuts <- how(nperm = 999)

# run adonis2
adonis_test[["DNA_RNA"]] <- adonis2(transposed_taxa~CollectionSite+DNA_RNA, data=adonis_meta, permutations=permuts, method="euclidean", by="terms")

# run betadisper for treatment
dispersion_dna_rna <- betadisper(vegdist(transposed_taxa, method="euclidean"), group=adonis_meta$Treatment)
# plot function
box_dna_rna <- function() { boxplot(dispersion_dna_rna) }
# export plot
export_figs_tabs(paste0(save_adonis, "box_dna_rna_", exp_name), box_dna_rna(), base=TRUE, width=168*2, height=168*3)
# plot function
disper_dna_rna <- function() { plot(dispersion_dna_rna, hull=TRUE, ellipse=FALSE) }
# export plot
export_figs_tabs(paste0(save_adonis, "disper_dna_rna_", exp_name), disper_dna_rna(), base=TRUE, width=168*2, height=168*3)

sink(paste0(save_adonis, "adonis_DNA_RNA.txt"))

# print results
cat("\n")
cat(paste0("Adonis test for DNA versus RNA"))
cat("\n")
print(adonis_test[["DNA_RNA"]])

cat("\n")
cat(paste0("Permutest test for betadisper on dna and rna"))
cat("\n")

print(permutest(dispersion_dna_rna))

cat("\n")
cat(paste0("Anova test for betadisper on dna and rna"))
cat("\n")

print(anova(dispersion_dna_rna))

sink()

# run adonis2
adonis_test[["t_by_t_DNA_RNA"]] <- adonis2(transposed_taxa~CollectionSite+DNA_RNA+Treatment*TimePoint, data=adonis_meta, permutations=permuts, method="euclidean", by="terms")

# run betadisper for treatment
dispersion_full_model <- betadisper(vegdist(transposed_taxa, method="euclidean"), group=adonis_meta$Treatment)
# plot function
box_full_model <- function() { boxplot(dispersion_full_model) }
# export plot
export_figs_tabs(paste0(save_adonis, "box_full_model_", exp_name), box_full_model(), base=TRUE, width=168*2, height=168*3)
# plot function
disper_full_model <- function() { plot(dispersion_full_model, hull=TRUE, ellipse=FALSE) }
# export plot
export_figs_tabs(paste0(save_adonis, "disper_full_model_", exp_name), disper_full_model(), base=TRUE, width=168*2, height=168*3)

sink(paste0(save_adonis, "adonis_treat_by_time_DNA_RNA.txt"))

# print results
cat("\n")
cat(paste0("Adonis test for DNA versus RNA, comprising treatment and time"))
cat("\n")
print(adonis_test[["t_by_t_DNA_RNA"]])

cat("\n")
cat(paste0("Permutest test for betadisper on full model"))
cat("\n")

print(permutest(dispersion_full_model))

cat("\n")
cat(paste0("Anova test for betadisper on full model"))
cat("\n")

print(anova(dispersion_full_model))

sink()

# export results
saveRDS(adonis_test, paste0(save_adonis, "adonis_test_clr_", exp_name, "_FIXED_SITE_CORRECTED.Rds"))

