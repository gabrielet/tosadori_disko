# source("011_differential_abundance_linda.R")

# load libraries
library("tidyverse")
library("reshape2")
library("ggrepel")
library("ggtext")
library("MicrobiomeStat")
library("phyloseq")
library("viridis")
library("scales")
library("RColorBrewer")

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
#taxa_algo <- "idtaxa" ; boot <- 80

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, " EXPERIMENT NAME ", exp_name))

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
save_diff_abundance <- paste0(save_img, "differential_abundance_new_tables_jan_2025/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_diff_abundance), "dir exists!", dir.create(save_diff_abundance, recursive=TRUE))

# print info regarding which directory is analysed, i.e. dada or qiime
print(paste0("the analysis is performed on ", toupper(path_to_exp)))

# load non-transformed object
phylo_data <- readRDS(paste0(save_phylo, "phylo_data_", exp_name, ".Rds"))

# get metadata
phylo_metadata <- sample_data(phylo_data)

######### DIFFERENTIAL ABUNDANCE

print("DIFFERENTIAL ABUNDANCE TREATMENT BY TIME BY DNA AND RNA")

# set interesting levels
interesting_levels <- c("Phylum", "Order", "Genus")

#################################################################### PERFORM ANALYSIS ON phylo_xna TO TEST DNA AND RNA SEPARATELY

# create a list of phylo objects containing DNA, RNA, or ALL samples
phylo_xna <- list()
phylo_xna[["DNA"]] <- prune_samples(as.character(phylo_metadata$sampleID[which(phylo_metadata$DNA_RNA == "DNA")]), phylo_data)
phylo_xna[["RNA"]] <- prune_samples(as.character(phylo_metadata$sampleID[which(phylo_metadata$DNA_RNA == "RNA")]), phylo_data)

# set threshold for significance
signif_t <- 0.05
# set log-fold-change threshold
lfc_thresh <- 1
# set p-value threshold
p_thresh <- abs(log10(signif_t))

# set abundance threshold
abundance_threshold <- 0

# plot padjusted or pvalue
pvaladj <- "padjusted"

# get palette
xna_palette <- brewer.pal(n=4, name = "Set2")[c(1,4)]
names(xna_palette) <- c("DNA", "RNA")

# storage for all results
all_res <- list()

# to throw an error instead of a warning set warn=2. normal behaviour is warn=0
# options(warn=2)

# loop through DNA and RNA
for (xna in names(phylo_xna)) {

	print(paste0("Analysing ", xna))

#	xna <- "DNA" ; taxa_level <- "Order"

	# remove ASVs which are always zero for the subset under consideration
	phylo_xna[[xna]] <- prune_taxa(taxa_sums(phylo_xna[[xna]]) > 0, phylo_xna[[xna]])

	# set variables of interests as factor
	# first level will be the one used by ANCOMBC as a reference level for the analysis
	sample_data(phylo_xna[[xna]])$Treatment <- factor(sample_data(phylo_xna[[xna]])$Treatment, levels=c("C", "W"))
	# set levels to compare to
	ref_lev_treat <- "C"
	sample_data(phylo_xna[[xna]])$Treatment <- relevel(sample_data(phylo_xna[[xna]])$Treatment, ref_lev_treat)
	sample_data(phylo_xna[[xna]])$TimePoint <- factor(sample_data(phylo_xna[[xna]])$TimePoint, levels=c("June", "July", "August"))
	# set levels to compare to
	ref_lev_time <- "July"
	sample_data(phylo_xna[[xna]])$TimePoint <- relevel(sample_data(phylo_xna[[xna]])$TimePoint, ref_lev_time)
	sample_data(phylo_xna[[xna]])$CollectionSite <- factor(sample_data(phylo_xna[[xna]])$CollectionSite, levels=unique(sample_data(phylo_xna[[xna]])$CollectionSite))

	# combine treatment and timepoint into a single variable
	treat_by_time <- with(sample_data(phylo_xna[[xna]]), interaction(Treatment, TimePoint))
	sample_data(phylo_xna[[xna]])$treat_by_time <- treat_by_time

	# allocate vector for future table
	all_significant_only <- vector()

	# loop through different taxa levels
	for (taxa_level in interesting_levels) {

		print(paste("treat by time at ", taxa_level))

		# agglomerate phylo
		if (taxa_level=="ASV") {
			# leave the table as is
			phylo_agg <- phylo_xna[[xna]]
		} else {
			# else aggregate
			phylo_agg <- tax_glom(phylo_xna[[xna]], taxrank=taxa_level)
		}

		# compute abundances
		phylo_abundances <- otu_table(transform_sample_counts(phylo_agg, function(x) x / sum(x)))

		# threshold for DA analysis
		above_thres_taxa <- rownames(phylo_abundances[apply(phylo_abundances, 1, mean) >= abundance_threshold, ]) ; updown <- "above"

		# get info using the thresholds on the abundances
		otu_data <- as(otu_table(phylo_agg), "matrix")[above_thres_taxa, ]
		# get metadata
		metadata <- as(sample_data(phylo_agg), "data.frame")

		# set seed
		set.seed(131)

		# run linda
		linda_outs <- linda(otu_data, metadata, formula = '~Treatment+TimePoint+(1|CollectionSite)', alpha=signif_t, prev.filter=0.06, n.cores=4) ; withint <- ""

		# show the available tables generated by linda
		# print(names(linda_outs$output))

		# get results from all the tables in linda_out
		for (which_tab in names(linda_outs$output)) {

#			which_tab <- names(linda_outs$output)[3]

			# get table
			volcano_data <- linda_outs$output[[which_tab]]#[linda_outs$output[[which_tab]]$reject==T, ]

			# if the test is against timepoint
			if(length(grep("Time", which_tab))==1) {
				# create label describing timepoint, and taxa level
				lbl <- paste0(updown, "_", abundance_threshold, "_", xna, "_", taxa_level, "_ref_", ref_lev_time, "_test_", which_tab)
			} else {
				# otherwise it is between control and warming
				lbl <- paste0(updown, "_", abundance_threshold, "_", xna, "_", taxa_level, "_ref_", ref_lev_treat, "_test_", which_tab)
			}

			# store results
			all_res[[lbl]] <- linda_outs$output[[which_tab]]

			# get taxa names for plotting purposes
			if (taxa_level=="ASV") {
				# if taxa_level is set at ASV level, get Genus level names
				volcano_data$taxa <- as(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), "Genus"], "character")
				volcano_data$table_taxa <- apply(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), c("Phylum", "Class", "Order", "Family", "Genus")], 1, paste, collapse="_")
			} else {
				# get taxa names at correct taxa level
				volcano_data$taxa <- as(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), taxa_level], "character")
				volcano_data$table_taxa <- apply(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), c("Phylum", "Class", "Order", "Family", "Genus")], 1, paste, collapse="_")
			}

			# make unique
			volcano_data$taxa <- make.unique(volcano_data$taxa)

			# if any result is available, plot
			if (any(volcano_data$log2FoldChange >= lfc_thresh & volcano_data$padj <= signif_t) | any(volcano_data$log2FoldChange <= -lfc_thresh & volcano_data$padj <= signif_t)) {

				# add to table with all significant
				significant_ones <- volcano_data[which(volcano_data$log2FoldChange >= lfc_thresh & volcano_data$padj <= signif_t | volcano_data$log2FoldChange <= -lfc_thresh & volcano_data$padj <= signif_t), ]
				# add significant ones
				all_significant_only <- rbind.data.frame(all_significant_only, cbind.data.frame(taxa_level, which_tab, significant_ones))

				# export table with volcano info
				write.table(cbind.data.frame(ASV=rownames(volcano_data), volcano_data[, c(9, 2:6)]), paste0(save_diff_abundance, "volcano_data_", withint, exp_name, "_", lbl, "_", pvaladj, "_WITH_TAXA_LEVS.csv"), row.names=F, col.names=T, sep="\t", quote=F)

				# if padj is plotted
				if(pvaladj=="padjusted") {

					# apply threshold
					volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= lfc_thresh & padj <= signif_t,"up", ifelse(log2FoldChange <= -lfc_thresh & padj <= signif_t, "down", "other")))
					volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= lfc_thresh & padj <= signif_t, 15, ifelse(log2FoldChange <= -lfc_thresh & padj <= signif_t, 15, 2.5)))
					# compute log of padj
					volcano_data$log_p <- abs(log10(volcano_data$padj))
				} else {
					# use pvalues

					# apply threshold
					volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= lfc_thresh & pvalue <= signif_t,"up", ifelse(log2FoldChange <= -lfc_thresh & pvalue <= signif_t, "down", "other")))
					volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= lfc_thresh & pvalue <= signif_t, 15, ifelse(log2FoldChange <= -lfc_thresh & pvalue <= signif_t, 15, 2.5)))
					# compute log of pvalue
					volcano_data$log_p <- abs(log10(volcano_data$pvalue))
				}

				# unlabel non significant dots, if any exist
				if (any(volcano_data$threshold=="other")) {
					volcano_data$taxa[which(volcano_data$threshold=="other")] <- ""
				}

				# remove taxa label
				volcano_data$taxa <- gsub(paste0(tolower(substr(taxa_level, 1, 1)), "__"), "", volcano_data$taxa)

				# reorganise for barplots
				significant_ones$threshold <- ifelse(significant_ones$log2FoldChange >= 0, "up", "down")
				significant_ones$taxa <- gsub(paste0(tolower(substr(taxa_level, 1, 1)), "__"), "", significant_ones$taxa)
				#significant_ones$taxa <- make.unique(significant_ones$taxa)
				significant_ones$taxa <- factor(significant_ones$taxa, levels=significant_ones$taxa[order(significant_ones$log2FoldChange)])

			} else {
				print("nothing significant!")
			}
		}
	}
	# modify label
	lbl <- paste0(updown, "_", abundance_threshold, "_", xna, "_ref_", ref_lev_time, "_test_", which_tab)

	if (length(all_significant_only) > 0) {
		print("EXPORTING TABLE ONE")

		# get table ready for export
#		all_significant_only$taxa_level
		all_significant_only$table_taxa <- gsub("_", ", ", gsub("[a-z]__", "", all_significant_only$table_taxa))

		all_significant_only$taxa <- NULL
		all_significant_only$baseMean <- NULL
		all_significant_only$lfcSE <- NULL
		all_significant_only$reject <- NULL
		all_significant_only$df<- NULL
		all_significant_only$log2FoldChange <- round(all_significant_only$log2FoldChange, digits = 4)
		all_significant_only$stat <- NULL
		all_significant_only$pvalue <- NULL
		all_significant_only$padj <- round(all_significant_only$padj, digits = 3)

		# move cols for publication-ready version
		all_significant_only_final <- all_significant_only[, c("which_tab", "taxa_level", "table_taxa", "log2FoldChange", "padj")]

		# export table with volcano info
		write.table(all_significant_only_final, paste0(save_diff_abundance, "WHAT_all_significant_only_", withint, exp_name, "_", lbl, "_", pvaladj, "_WITH_TAXA_LEVS.csv"), row.names=F, col.names=T, sep="\t", quote=F)
	}
}

################################################### PERFORM ANALYSIS ON phylo_data TO TEST DNA VS RNA

print("DNA vs RNA")

# allocate vector for future table
all_significant_only_drna <- vector()

# loop through different taxa levels
for (taxa_level in interesting_levels) {

#	taxa_level <- "Phylum"

	# set variables of interests as factor
	# first level will be the one used by ANCOMBC as a reference level for the analysis
	sample_data(phylo_data)$Treatment <- factor(sample_data(phylo_data)$Treatment, levels=c("C", "W"))
	sample_data(phylo_data)$TimePoint <- factor(sample_data(phylo_data)$TimePoint, levels=c("June", "July", "August"))
	sample_data(phylo_data)$CollectionSite <- factor(sample_data(phylo_data)$CollectionSite, levels=unique(sample_data(phylo_data)$CollectionSite))
	# set levels to compare to
	ref_lev_dnarna <- "DNA"
	sample_data(phylo_data)$DNA_RNA <- relevel(sample_data(phylo_data)$DNA_RNA, ref_lev_dnarna)

	print(paste0("treat by time at ", taxa_level))

	# agglomerate phylo
	if (taxa_level=="ASV") {
		# leave the table as is
		phylo_agg <- phylo_data
	} else {
		# else aggregate
		phylo_agg <- tax_glom(phylo_data, taxrank=taxa_level)
	}

	# get info
	otu_data <- as(otu_table(phylo_agg), "matrix")
	metadata <- as(sample_data(phylo_agg), "data.frame")

	# set seed
	set.seed(131)

	# this is to evaluate which value to assign to prev.filter
	# it basically counts the percentage of zero by row, i.e. by taxa
	# rowSums(otu_data==0)/ncol(otu_data)*100

	# run linda
	linda_outs <- linda(otu_data, metadata, formula = '~DNA_RNA+(1|CollectionSite)', alpha = signif_t, prev.filter=0.03, feature.dat.type="count", zero.handling="pseudo-count", n.cores=4)

	# get results
	volcano_data <- linda_outs$output[[1]]#[linda_outs$output$reject==T, ]

	# create label describing timepoint, and taxa level
	lbl <- paste0(taxa_level, "_ref_DNA_test_RNA_", abundance_threshold)

	# store results
	all_res[[lbl]] <- linda_outs$output[[which_tab]]

	# get taxa names for plotting purposes
	if (taxa_level=="ASV") {
		# if taxa_level is set at ASV level, get Genus level names
		volcano_data$taxa <- as(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), "Genus"], "character")
		volcano_data$table_taxa <- apply(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), c("Phylum", "Class", "Order", "Family", "Genus")], 1, paste, collapse="_")
	} else {
		# get taxa names at correct taxa level
		volcano_data$taxa <- as(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), taxa_level], "character")
		volcano_data$table_taxa <- apply(tax_table(phylo_agg)[which(rownames(tax_table(phylo_agg)) %in% rownames(volcano_data)), c("Phylum", "Class", "Order", "Family", "Genus")], 1, paste, collapse="_")
	}
	# make unique
	volcano_data$taxa <- make.unique(volcano_data$taxa)

	# if any result is available, plot
	if (any(volcano_data$log2FoldChange >= lfc_thresh & volcano_data$padj <= signif_t) | any(volcano_data$log2FoldChange <= -lfc_thresh & volcano_data$padj <= signif_t)) {

		# add to table with all significant
		significant_ones_xna <- volcano_data[which(volcano_data$log2FoldChange >= lfc_thresh & volcano_data$padj <= signif_t | volcano_data$log2FoldChange <= -lfc_thresh & volcano_data$padj <= signif_t), ]
		# add significant ones
		all_significant_only_drna <- rbind.data.frame(all_significant_only_drna, cbind.data.frame(taxa_level, significant_ones_xna))

		# export table with volcano info
		write.table(cbind.data.frame(ASV=rownames(volcano_data), volcano_data[, c(9, 2:6)]), paste0(save_diff_abundance, "volcano_data_", exp_name, "_", lbl, "_", pvaladj, ".csv"), row.names=F, col.names=T, sep="\t", quote=F)

		# if padj is plotted
		if(pvaladj=="padjusted") {

			# apply threshold
			volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= lfc_thresh & padj <= signif_t,"up", ifelse(log2FoldChange <= -lfc_thresh & padj <= signif_t, "down", "other")))
			volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= lfc_thresh & padj <= signif_t, 15, ifelse(log2FoldChange <= -lfc_thresh & padj <= signif_t, 15, 5)))

			# compute log of padj
			volcano_data$log_p <- abs(log10(volcano_data$padj))
		} else {
			# plot pvalues
			# apply threshold
			volcano_data <- volcano_data %>% mutate(threshold = ifelse(log2FoldChange >= lfc_thresh & pvalue <= signif_t,"up", ifelse(log2FoldChange <= -lfc_thresh & pvalue <= signif_t, "down", "other")))
			volcano_data <- volcano_data %>% mutate(dsize=ifelse(log2FoldChange >= lfc_thresh & pvalue <= signif_t, 15, ifelse(log2FoldChange <= -lfc_thresh & pvalue <= signif_t, 15, 5)))

			# compute log of pvalue
			volcano_data$log_p <- abs(log10(volcano_data$pvalue))
		}

		# unlabel non significant dots, if any exist
		if (any(volcano_data$threshold=="other")) {
			volcano_data$taxa[which(volcano_data$threshold=="other")] <- ""
		}

		# remove taxa label
		volcano_data$taxa <- gsub(paste0(tolower(substr(taxa_level, 1, 1)), "__"), "", volcano_data$taxa)

		# reorganise for barplots
		significant_ones_xna$threshold <- ifelse(significant_ones_xna$log2FoldChange >= 0, "up", "down")
		significant_ones_xna$taxa <- gsub(paste0(tolower(substr(taxa_level, 1, 1)), "__"), "", significant_ones_xna$taxa)
		#significant_ones_xna$taxa <- make.unique(significant_ones_xna$taxa)
		significant_ones_xna$taxa <- factor(significant_ones_xna$taxa, levels=significant_ones_xna$taxa[order(significant_ones_xna$log2FoldChange)])

	} else {
		print("nothing significant!")
	}
}

# modify label for all_significant_only_drna
lbl <- paste0("ref_DNA_test_RNA_", abundance_threshold)

print("EXPORTING TABLE TWO")

# get table ready for export
all_significant_only_drna$table_taxa <- gsub("_", ", ", gsub("[a-z]__", "", all_significant_only_drna$table_taxa))

# adding column for taxa names
all_significant_only_drna$t_names <- ""
all_significant_only_drna$t_names[grep("Unidentified", all_significant_only_drna$taxa)] <- paste0("Unidentified ", all_significant_only_drna$taxa_level[grep("Unidentified", all_significant_only_drna$taxa)])
all_significant_only_drna$t_names[grep("Unidentified", all_significant_only_drna$taxa, invert=T)] <- gsub("[p|o|g]__", "", all_significant_only_drna$taxa[grep("Unidentified", all_significant_only_drna$taxa, invert=T)])

all_significant_only_drna$taxa <- NULL
all_significant_only_drna$baseMean <- NULL
all_significant_only_drna$lfcSE <- NULL
all_significant_only_drna$reject <- NULL
all_significant_only_drna$df<- NULL
all_significant_only_drna$log2FoldChange <- round(all_significant_only_drna$log2FoldChange, digits = 4)
all_significant_only_drna$stat <- NULL
all_significant_only_drna$pvalue <- NULL
all_significant_only_drna$padj <- round(all_significant_only_drna$padj, digits = 3)

# move cols for publication-ready version
all_significant_only_drna_final <- all_significant_only_drna[, c("t_names", "table_taxa", "log2FoldChange", "padj")]

# export table with volcano info
write.table(all_significant_only_drna_final, paste0(save_diff_abundance, "WHAT_all_significant_only_drna_", exp_name, "_", lbl, "_", pvaladj, ".csv"), row.names=F, col.names=T, sep="\t", quote=F)

print("LinDA analysis is done, export rds")

# export all the data thata were computed
saveRDS(all_res, paste0(save_diff_abundance, "WHAT_linda_out_", ref_lev_time, "_", ref_lev_treat, "_treat_by_time_and_site_threshold_", lfc_thresh, "_", pvaladj, ".Rds"))

