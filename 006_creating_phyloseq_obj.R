# source("006_creating_phyloseq_obj_bacteria_dada.R")

# load libraries
library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("ellipse")
library("ggpubr")
library("ggsignif")
library("ggrepel")
library("ggtext")
library("gridExtra")
library("phyloseq")
library("vegan")

# import functions
source("/microbiology/disko2013/code/000_micro_functions_disko2013.R")

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
taxa_algo <- "NBC"

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, "  EXPERIMENT NAME ", exp_name, " CLUSTERED WITH ", clust_method, " TAXA ASSIGNED WITH ", taxa_algo))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/microbiology/disko2013/"
# this is the general path to the experiment
path_to_exp <- paste0(root_path, orgn_dir, "experiments/", exp_name, "/")
# this is the path where counts results will be stored
path_to_taxa <- paste0(path_to_exp, clust_method, "/")
save_img <- paste0(path_to_taxa, "figures/", taxa_algo, "/")
save_phylo <- paste0(path_to_taxa, "phyloseq_renamed_test/", taxa_algo, "/")

# this is the path to the outputs of dada
path_to_counts <- paste0(path_to_exp, "dada_files/tables/")

# if directories do not exist, create it
ifelse(dir.exists(save_img), TRUE, dir.create(save_img, recursive=T))
ifelse(dir.exists(save_phylo), TRUE, dir.create(save_phylo, recursive=T))

# and to taxonomy
path_to_taxonomy <- paste0(path_to_exp, clust_method, "/taxonomy/")

########################### select threshold for filtering bad assignments ###########################

print("IMPORTING TAXONOMY AND COUNTS TABLES")

# NBC bootstrap
boot_filter <- 80

# import counts
dada_counts <- readRDS(paste0(path_to_counts, "counts_table_ASV.Rds"))

# rename rownames to match metadata sample names
rownames(dada_counts) <- sapply(strsplit(rownames(dada_counts), "_"), "[", 1)

# taxonomy created through naive bayes classifier
dada_taxonomy <- readRDS(paste0(path_to_taxonomy, "taxonomy_table_boot_", boot_filter, "_", taxa_algo, "_", clust_method, ".Rds"))

######### CLEAN TAXA TABLE AND UNIDENTIFIED ORGANISMS #########

print("CLEANING ASVs/OTUs TABLE")

# subset taxonomy table
taxa_tab_clean <- dada_taxonomy[which(rownames(dada_taxonomy) %in% colnames(dada_counts)), ]

if (orgn == "bacteria") {
	# get bacteria and archaea only from taxa table
	bacteria_only <- grep("Bacteria|Archaea", taxa_tab_clean[, "Kingdom"])
	# if we have at least one bacteria
	if (length(bacteria_only) > 1) {
		# get it
		taxa_tab_clean <- taxa_tab_clean[bacteria_only, ]
	}
} else if (orgn == "fungi") {
	# get fungi only from taxa table
	fungi_only <- grep("Fungi", taxa_tab_clean[, "Kingdom"])
	# if we have at least one fungus
	if (length(fungi_only) > 1) {
		# get it
		taxa_tab_clean <- taxa_tab_clean[fungi_only, ]
	}
}

print(paste0("ASVs/OTUs we are actually using after selecting fungi only ", nrow(taxa_tab_clean), " which correspond to ", length(unique(taxa_tab_clean[, "Genus"])), " Genus"))

# transpose dada_counts for good
dada_counts <- t(dada_counts)

# remove undefined ASVs using rownames from counts
dada_counts <- dada_counts[which(rownames(dada_counts) %in% rownames(taxa_tab_clean)), ]

# print some info for the paper
print("SOME INFO FOR THE PAPER AFTER THE CLEANING")

print(paste0("DNA AVERAGE READS ", ceiling(mean(colSums(dada_counts)[grep("AD006", colnames(dada_counts))]))))
print(paste0("RNA AVERAGE READS ", ceiling(mean(colSums(dada_counts)[grep("AD012", colnames(dada_counts))]))))
print(paste0("DNA TOTAL READS ", sum(colSums(dada_counts)[grep("AD006", colnames(dada_counts))])))
print(paste0("RNA TOTAL READS ", sum(colSums(dada_counts)[grep("AD012", colnames(dada_counts))])))

# search for mitochondria and chloroplasts, in the bacterial taxa
is_mito <- vector()
is_chloro <- vector()

# loop through rows in taxa_tab_clean. one row is one taxa
# we have a total of i rows
for (i in c(1:nrow(taxa_tab_clean))) {
	# loop through cols. on col is on taxonomic level
	# we have j taxonomic levels. i.e. j will range in [1-7]
	for (j in c(1:ncol(taxa_tab_clean))) {
		# create label specifying the taxonomic level we are currently analysing
		tlabel <- tolower(substr(colnames(taxa_tab_clean)[j], 1, 1))
		# if NA, then unassigned, uncultured or unknown is found
		if (length(grep("NA|nassigned|[u|U]ncultured|nknown|nidentified", taxa_tab_clean[i, j])) > 0 | is.na(taxa_tab_clean[i, j])) {
			# assign the Unidentified label at the specific, undidentified
			# taxonomic level
			taxa_tab_clean[i, j] <- paste0(tlabel, "__Unidentified")
		# if mitochondria is found
		} else if (length(grep("itochond", taxa_tab_clean[i, j])) > 0) {
			# store its position, i.e. the row number
			is_mito <- append(is_mito, i)
		# if chloroplast is found
		} else if (length(grep("hloropla", taxa_tab_clean[i, j])) > 0) {
			# store its position, i.e. the row number
			is_chloro <- append(is_chloro, i)
		# if a suspicious case is found
		} else {
			# add tlabel if the taxa_algo was NBC
			if (taxa_algo=="NBC" | taxa_algo=="BLAST") {
				taxa_tab_clean[i, j] <- paste0(tlabel, "__", taxa_tab_clean[i, j])
			}
		}
	}
}

# remove mito and chloro from taxa_tab_clean, all at once
to_be_removed <- c(unique(is_mito), unique(is_chloro))
if (length(to_be_removed) > 1) {
	print(paste0(length(unique(is_mito)), " mitochondria found. removing them!"))
	print(paste0(length(unique(is_chloro)), " chloroplasts found. removing them!"))
	# remove rows that should be removed
	taxa_tab_clean <- taxa_tab_clean[-to_be_removed, ]
}

# remove non-bacterial ASVs using rownames from counts
# counts are now clean and ready to use
dada_counts <- dada_counts[which(rownames(dada_counts) %in% rownames(taxa_tab_clean)), ]

# remove weird taxonomy levels, i.e. g__g__, or similar
if (orgn == "fungi") {
	taxa_tab_clean <- gsub("k__k__", "k__", taxa_tab_clean)
	taxa_tab_clean <- gsub("p__p__", "p__", taxa_tab_clean)
	taxa_tab_clean <- gsub("c__c__", "c__", taxa_tab_clean)
	taxa_tab_clean <- gsub("o__o__", "o__", taxa_tab_clean)
	taxa_tab_clean <- gsub("f__f__", "f__", taxa_tab_clean)
	taxa_tab_clean <- gsub("g__g__", "g__", taxa_tab_clean)
	taxa_tab_clean <- gsub("s__s__", "s__", taxa_tab_clean)
}

# export the full taxa table, now cleaned
write.table(taxa_tab_clean, paste0(save_phylo, "taxa_tab_clean_DNA_and_RNA_", exp_name, ".csv"), row.names=F, col.names=TRUE, sep="\t")

# print out percentage of unidentified
for (i in 1:6) {
	perc <- round((length(grep("Unidenti", taxa_tab_clean[, i]))/nrow(taxa_tab_clean))*100, digits=1)
	print(paste0("the percentage of unidentified at ", colnames(taxa_tab_clean)[i], " level is equal to ", perc))
}

######### GET COUNTS AND CLEAN IT BEFORE CREATING PHYLOSEQ OBJS THE CREATE PHYLO OBJS #########

print("CREATING PHYLOSEQ OBJECTS")

# load DNA metadata, the new version from Jana
metadata_dna <- as.data.frame(read.csv("/microbiology/disko2013/metadata/csv_original/Disko2013_warming_control_metadata_updated_removed_few_vars.csv", sep="\t", header=T, stringsAsFactors=TRUE))

# if orgn is fungi, modify the labels of the sample names
if (orgn == "fungi") {
	metadata_dna$sampleID <- gsub("B0", "F0", metadata_dna$sampleID)
	metadata_dna$sampleID <- gsub("B1", "F1", metadata_dna$sampleID)
}

# metadata for RNA samples are also needed. to create them
# start by creating a copy of the data from DNA samples, since
# the parameters are the same for both DNA and RNA
metadata_rna <- metadata_dna

# create new names for RNA samples, what was called "AD006"
# now become "AD012"
metadata_rna$sampleID <- as.factor(gsub("AD006", "AD012", metadata_dna$sampleID))
metadata_rna$DNA_RNA <- as.factor(rep("RNA", length(metadata_dna$DNA_RNA)))

# combine the DNA and RNA tables in a single table with shared cols
combined_exps_metadata <- rbind.data.frame(metadata_dna, metadata_rna)

# drop useless columns
combined_exps_metadata$mergingIDs <- NULL
combined_exps_metadata$sample <- NULL
combined_exps_metadata$Bacteria_Fungi <- NULL
combined_exps_metadata$ligation_tag <- NULL
combined_exps_metadata$dual_labeled_primer_tag <- NULL
combined_exps_metadata$TreatBYTime <- NULL

# add a column with the ratio between %C and %N
combined_exps_metadata$CNratio <- combined_exps_metadata$perc_of_C / combined_exps_metadata$perc_of_N

# recompute the unit
# before: copy numbers/mg of fresh soil
# after: copy numbers/g of dry soil

# compute the new unit: copy/mg_fresh x 1000 = copy/g_fresh
# then: 		copy/g_fresh x g_fresh / g_dry = copy
combined_exps_metadata$mcrA <- (combined_exps_metadata$mcrA * 1000) * (combined_exps_metadata$soil_fresh_weight_g/combined_exps_metadata$soil_dry_weight_g)
combined_exps_metadata$mxaF <- (combined_exps_metadata$mxaF * 1000) * (combined_exps_metadata$soil_fresh_weight_g/combined_exps_metadata$soil_dry_weight_g)
combined_exps_metadata$nifH <- (combined_exps_metadata$nifH * 1000) * (combined_exps_metadata$soil_fresh_weight_g/combined_exps_metadata$soil_dry_weight_g)
combined_exps_metadata$nirS <- (combined_exps_metadata$nirS * 1000) * (combined_exps_metadata$soil_fresh_weight_g/combined_exps_metadata$soil_dry_weight_g)
combined_exps_metadata$nosZ <- (combined_exps_metadata$nosZ * 1000) * (combined_exps_metadata$soil_fresh_weight_g/combined_exps_metadata$soil_dry_weight_g)

# before constructing a pyhyloseq object, are the following requirements met?

# 1. taxaTable is a matrix. Must specify if the species are rows or columns,
# i.e. taxa_are_rows=TRUE/FALSE
# 2. metaData can be data.frame or whatever. Its rownames must match sample_names(taxaTable)
# 3. Taxa is a character matrix. Its rownames must match taxa_names(taxaTable)

# get DNA and RNA metadata
meta_phylo <- sample_data(combined_exps_metadata)
rownames(meta_phylo) <- meta_phylo$sampleID

# get taxa table, i.e. taxa_tab_clean
taxa_phylo <- tax_table(taxa_tab_clean)

# remove any ASV with fewer 0.1% of the mean number of reads per sample per library size
dada_counts_clean <- otu_table(dada_counts[(rowSums(dada_counts) > round(mean(colSums(dada_counts))*0.001)), ], taxa_are_rows=TRUE)

# apply transformation to raw counts
# 1. NO TRANSFORMATION
# 2. HELLINGER TRANSFORMATION ON ABUNDANCES
# 3. CENTERED LOG RATIO

# 1. build phyloseq object with raw data
phylo_data <- phyloseq(otu_table(dada_counts_clean, taxa_are_rows=TRUE), meta_phylo, taxa_phylo)

# remove taxa that are lost because of the subsampling:
# all the treatments were used to build the ASV table but
# we only kept Warming and Controls, hence some taxa that
# were only found in, for instance, Shrub-Removal samples
# are now all zero and we need to remove them.
# as mentioned here https://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html
# indeed `any(taxa_sums(phylo_data) == 0)` being TRUE is
# something to investigate
phylo_data <- prune_taxa(taxa_sums(phylo_data) > 0, phylo_data)

# 2. transform counts using hellinger transformation
# https://github.com/joey711/phyloseq/issues/585
hellinger_phylo_data <- transform_sample_counts(phylo_data, function(x) sqrt(x / sum(x)))

# 3. apply robust centered log ratio transformation using microbiome package function
centered_counts_phylo <- microbiome::transform(otu_table(phylo_data), 'clr')

# 4. build phyloseq object with CLRed counts
# NOTE THAT FILTERING WITH CENTERED DATA IS TRICKY BE CAREFUL HERE.
# also the n. of taxa here is not trimmed so it is expected this obj
# is different from the hellinger and untransformed counts
centered_phylo_data <- phyloseq(centered_counts_phylo, meta_phylo, taxa_phylo)

############ EXPORTING PHYLOSEQ OBJECTS ############

# export non-transformed object
saveRDS(phylo_data, paste0(save_phylo, "phylo_data_", exp_name, ".Rds"))

# export hellinger transformed object
saveRDS(hellinger_phylo_data, paste0(save_phylo, "hellinger_phylo_data_", exp_name, ".Rds"))

# export CLR transformed object
saveRDS(centered_phylo_data, paste0(save_phylo, "centered_phylo_data_", exp_name, ".Rds"))

# plot distribution of total reads per sample using the non-transformed data
# see https://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html

# create data frame
read_sums <- data.frame(nreads=sort(taxa_sums(phylo_data), TRUE), sorted=1:ntaxa(phylo_data), type="OTUs")
read_sums <- rbind(read_sums, data.frame(nreads=sort(sample_sums(phylo_data), TRUE), sorted=1:nsamples(phylo_data), type="Samples"))

# get the plot ready
plot_read_num <- ggplot(read_sums, aes(x=sorted, y=nreads)) +
	geom_bar(stat="identity") +
	ggtitle("Total number of reads") +
	scale_y_log10() +
	facet_wrap(~type, 1, scales="free") +
	theme(text=element_text(size=15, face="bold"), legend.position="right") +
	theme_bw()

# plot
export_svg(paste0(save_img, "read_distribution_", exp_name, ".svg"), plot_read_num)

######### CREATING QPCR DATA TABLE AND PHYLOSEQ OBJECT #########

# NOT THAT THE qPCR phyloseq object will contain less taxa
# just because it was created using those taxa found for DNA
# only. So, if phylo_data is filtered on DNA, then the dimensions
# of phylo_data and phylo_data_qpcr must match

print("CREATING QPCR PHYLOSEQ OBJECTS")

# extract info from the original metadata, regarding genes expression
qpcr_full_original <- metadata_dna[, c("sampleID", "DNA_RNA", "Treatment", "TimePoint", "CollectionSite", "mcrA", "mxaF", "nifH", "nirS", "nosZ", "soil_fresh_weight_g", "soil_dry_weight_g")]

# get samples for Warming and Control only
qpcr_metadata <- qpcr_full_original[which(qpcr_full_original$sampleID %in% meta_phylo$sampleID), ]

# recompute the unit
# before: copy numbers/mg of fresh soil
# after: copy numbers/g of dry soil

# compute the new unit: copy/mg_fresh x 1000 = copy/g_fresh
# then: 		copy/g_fresh x g_fresh / g_dry = copy
qpcr_metadata$mcrA <- (qpcr_metadata$mcrA * 1000) * (qpcr_metadata$soil_fresh_weight_g/qpcr_metadata$soil_dry_weight_g)
qpcr_metadata$mxaF <- (qpcr_metadata$mxaF * 1000) * (qpcr_metadata$soil_fresh_weight_g/qpcr_metadata$soil_dry_weight_g)
qpcr_metadata$nifH <- (qpcr_metadata$nifH * 1000) * (qpcr_metadata$soil_fresh_weight_g/qpcr_metadata$soil_dry_weight_g)
qpcr_metadata$nirS <- (qpcr_metadata$nirS * 1000) * (qpcr_metadata$soil_fresh_weight_g/qpcr_metadata$soil_dry_weight_g)
qpcr_metadata$nosZ <- (qpcr_metadata$nosZ * 1000) * (qpcr_metadata$soil_fresh_weight_g/qpcr_metadata$soil_dry_weight_g)

# get combined metadata
meta_phylo_qpcr <- sample_data(qpcr_metadata)
rownames(meta_phylo_qpcr) <- meta_phylo_qpcr$sampleID

# apply transformation to raw counts
# 1. NO TRANSFORMATION
# 2. HELLINGER TRANSFORMATION ON ABUNDANCES
# 3. CENTERED LOG RATIO

# 1 build phyloseq object, using dada_counts_clean
phylo_data_qpcr <- phyloseq(otu_table(dada_counts_clean, taxa_are_rows=TRUE), meta_phylo_qpcr, taxa_phylo)

# remove taxa that are lost because of the subsampling:
# all the samples were used to build the ASV table but here
# we only kept Warming and Controls, hence some taxa that
# were only found in, for instance, Shrub-Removal samples
# are now all zero and we need to remove them.
# as mentioned here https://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html
# indeed `any(taxa_sums(phylo_data) == 0)` being TRUE is
# something to investigate
phylo_data_qpcr <- prune_taxa(taxa_sums(phylo_data_qpcr) > 0, phylo_data_qpcr)

# IMPORTANT: the number of taxa for qPCR is expected to be lower than
# the n. of taxa in the phylo_data objects. This is because for qPCR
# we are using DNA-taxa only while phylo_data also comprise RNA-taxa
# so this object should be same as phylo_data_qpcr:
# try_me <-prune_samples(as.character(sample_data(phylo_data)$sampleID[which(sample_data(phylo_data)$DNA_RNA == "DNA")]), phylo_data); prune_taxa(taxa_sums(try_me) > 0, try_me)

# 1.1 transform counts using hellinger transformation
# https://github.com/joey711/phyloseq/issues/585
hellinger_phylo_qpcr <- transform_sample_counts(phylo_data_qpcr, function(x) sqrt(x / sum(x)))

# 2 transform counts using centered log ratio transformation
# use counts already centered in the phylo_data creation
centered_counts_qpcr <- microbiome::transform(otu_table(phylo_data_qpcr), 'clr')

# 2.1 build phyloseq object with CLRed counts
# NOTE THAT FILTERING WITH CENTERED DATA IS TRICKY. BE CAREFUL HERE.
# also the n. of taxa here is not trimmed so it is expected this obj
# is different from the hellinger and untransformed counts
centered_phylo_qpcr <- phyloseq(centered_counts_qpcr, meta_phylo_qpcr, taxa_phylo)

############ EXPORTING ############

# export non-transformed object
saveRDS(phylo_data_qpcr, paste0(save_phylo, "phylo_qpcr_", exp_name, ".Rds"))

# export hellinger transformed object
saveRDS(hellinger_phylo_qpcr, paste0(save_phylo, "hellinger_phylo_qpcr_", exp_name, ".Rds"))

# export CLR transformed object
saveRDS(centered_phylo_qpcr, paste0(save_phylo, "centered_phylo_qpcr_", exp_name, ".Rds"))

# PRINT OUT SOME IMPORTANT INFO
cat("\n")
cat("IMPORTANT NOTE: data are exported as raw counts, hellinger and centered-log-ratio transformed!")
cat("\n")

