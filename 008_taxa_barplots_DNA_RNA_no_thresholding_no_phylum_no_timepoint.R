# source("008_taxa_barplots_DNA_RNA_no_thresholding_no_phylum_no_timepoint.R")

# load libraries
library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("phyloseq")
library("viridis")
library("scales")
library("plotrix")

# import functions
source("/mnt/cinqueg/gabriele/work/microbiology/disko2013/code/000_micro_functions_disko2013.R")

######### SELECT EXPERIMENT AND DEFINE PATHS #########
# select exp name
exp_name <- "submission_final_norejects_dada_new_bootstrap" ; orgn <- "bacteria"

# select the clustering method
if (orgn == "bacteria") {
	clust_method <- "ASV"
} else if (orgn == "fungi") {
	clust_method <- "OTU"
}

# select the taxa assigment method
taxa_algo <- "NBC"

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, "  EXPERIMENT NAME ", exp_name, " CLUSTERED WITH ", clust_method, " TAXA ASSIGNED WITH ", taxa_algo))

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
save_taxa_figs <- paste0(save_img, "taxaFigs/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_taxa_figs), "dir exists!", dir.create(save_taxa_figs, recursive=TRUE))

# load non-transformed object
phylo_data <- readRDS(paste0(save_phylo, "phylo_data_", exp_name, ".Rds"))

# get metadata table
phylo_metadata <- sample_data(phylo_data)

# loop through taxa levels, to aggregate appropriately
# get best 10 taxa
n_best <- c(6, 15, 15)
names(n_best) <- c("Phylum", "Order", "Genus")

######### PLOTTING BEST SCORING TAXA

# loop through taxa_levels
for (taxa_level in c("Phylum", "Order")) {

	# aggregate taxa at a specific taxonomy level
	phylo_aggregated <- tax_glom(phylo_data, taxrank=taxa_level)

	# get metadata for phylo_aggregated
	phylo_aggregated_meta <- sample_data(phylo_aggregated)

	# set levels to have controls first
	sample_data(phylo_aggregated)$Treatment <- factor(sample_data(phylo_aggregated)$Treatment, levels=c("C", "W"))

	# compute sums to obtain the contribution of all sampling sites
	# to a specific taxa, then order from most abundant to less abundant
	xna_sums <- sort(rowSums(otu_table(phylo_aggregated)), decreasing=T)

	# get ordered taxa names, using factors
	taxa_names <- tax_table(phylo_aggregated)[names(xna_sums), taxa_level]

	# find unidenti, best scoring taxa, and other taxa
	unidenti <- grep("Unidentified", taxa_names)

	# define unidenti_taxa as NULL to test whether some taxa were
	# actually unidentified
	unidenti_taxa <- NULL

	# same for other and best
	other_taxa <- NULL
	best_taxa <- NULL

	# test unidentified
	if (length(unidenti)>0){
		# get taxa_names
		unidenti_taxa <- taxa_names[unidenti]
		# remove unidentified taxa
		taxa_names <- taxa_names[-unidenti]
	}
	# if n_best is more than available taxa
	if (n_best[taxa_level] > nrow(taxa_names)) {
		# throw an error
		stop(paste0("SELECTED N_BEST ARE MORE THAN THE AVAILABLE TAXA, WHICH ARE ", nrow(taxa_names)))
	} else {
		# otherwise get best taxa
		best_taxa <- taxa_names[1:n_best[taxa_level]]
	}
	# if n_best is less than available taxa
	if (n_best[taxa_level] < nrow(taxa_names)) {
		# get other taxa
		other_taxa <- taxa_names[(n_best[taxa_level]+1):nrow(taxa_names)]
	} else {
		other_taxa <- ""
	}

	################################################################################################

	# define some empty lists to store all results
	storing_all_xnas <- vector()

	# loop through treatment
	for (xna in unique(sample_data(phylo_aggregated)$DNA_RNA)) {

#		xna <- "DNA"

		print(paste0("Analysing ", xna))

		# select samples belonging to a treatment
		phylo_xna <- prune_samples(as.character(sample_data(phylo_aggregated)$sampleID[which(sample_data(phylo_aggregated)$DNA_RNA == xna)]), phylo_aggregated)

		# remove ASVs which are always zero for the subset under consideration
		phylo_xna <- prune_taxa(taxa_sums(phylo_xna) > 0, phylo_xna)

		# compute relative abundance
		phylo_xna <- transform_sample_counts(phylo_xna, function(x) x / sum(x))

		# define empty to store single xna
		storing_single_xna <- as.data.frame(matrix(ncol=8))
		colnames(storing_single_xna) <- c("T_level", "DNA_RNA", "Relative", "SE", "SE_up", "SE_down")

		# loop through taxa types
		for (taxa_type in c("best_taxa", "other_taxa", "unidenti_taxa")) {

#			taxa_type <- "other_taxa"

			print(paste0("Analysing ", taxa_type))

			taxa_names_sub <- ""

			if (taxa_type=="best_taxa") {
				# get best taxa
				treat_sums <- rowMeans(otu_table(phylo_xna)[which(tax_table(phylo_xna)[, taxa_level] %in% get(taxa_type)), ])
				# compute standard errors
				all_se <- apply(otu_table(phylo_xna)[which(tax_table(phylo_xna)[, taxa_level] %in% get(taxa_type)), ], 1, std.error)
				standerr_up <- treat_sums + all_se
				standerr_down <- treat_sums - all_se
				# get taxa names
				taxa_names_sub <- tax_table(phylo_xna)[names(treat_sums), taxa_level]

				# put columns together for phylo_xna samples and compute
				# the percentage of contribution each taxa is providing to
				# the total, to finally compare the taxa. by doing so, we
				# make sure the total will always sum up to 1, i.e. 100%
				storing_single_xna <- cbind.data.frame(T_level=taxa_names_sub, DNA_RNA=xna, Relative=treat_sums, SE=all_se, SE_up=standerr_up, SE_down=standerr_down)
				rownames(storing_single_xna) <- NULL
				colnames(storing_single_xna) <- c("T_level", "DNA_RNA", "Relative", "SE", "SE_up", "SE_down")

				# sort table by abundance
				storing_single_xna <- storing_single_xna[order(storing_single_xna$Relative, decreasing=T), ]
			} else {
				# get taxa names
				if (taxa_type=="other_taxa") {
					taxa_names_sub <- "Other"
				} else {
					taxa_names_sub <- "Unidentified"
				}
				# compute sums to obtain the contribution of all sampling sites
				# to a specific taxa, then order from most abundant to less abundant
				# here, Camelia would compute the average from the sites (see the
				# file with the example computations)
				# test if the taxa are found in the subset
				if (any(tax_table(phylo_xna)[, taxa_level] %in% get(taxa_type))) {
					# get taxa sums
					treat_sums <- sum(rowMeans(otu_table(phylo_xna)[which(tax_table(phylo_xna)[, taxa_level] %in% get(taxa_type)), ]))
					# compute standard errors
					all_se <- std.error(c(otu_table(phylo_xna)[which(tax_table(phylo_xna)[, taxa_level] %in% get(taxa_type)), ]))
					standerr_up <- treat_sums + all_se
					standerr_down <- treat_sums - all_se

					# put columns together for phylo_xna samples and compute
					# the percentage of contribution each taxa is providing to
					# the total, to finally compare the taxa. by doing so, we
					# make sure the total will always sum up to 1, i.e. 100%
					storing_single_xna <- cbind.data.frame(T_level=taxa_names_sub, DNA_RNA=xna, Relative=treat_sums, SE=all_se, SE_up=standerr_up, SE_down=standerr_down)
					rownames(storing_single_xna) <- NULL
					colnames(storing_single_xna) <- c("T_level", "DNA_RNA", "Relative", "SE", "SE_up", "SE_down")
				} else {
					storing_single_xna <- NULL
				}
			}
			# add to storing_all_xnas
			storing_all_xnas <- rbind.data.frame(storing_all_xnas, storing_single_xna)
		}
	}

	# remove taxa_level label
	which_level <- paste0(tolower(substr(taxa_level, 1, 1)), "__")
	storing_all_xnas$T_level <- gsub(which_level, "", storing_all_xnas$T_level)

	# order levels
	last_two <- as.character(c(storing_all_xnas$T_level[!grepl("Other|Unidenti", storing_all_xnas$T_level)], storing_all_xnas$T_level[grepl("Other|Unidenti", storing_all_xnas$T_level)]))
	storing_all_xnas$T_level <- factor(storing_all_xnas$T_level, levels=unique(last_two))

	print("plotting")

	# re-set palette
	plotting_palette <- c("antiquewhite", "aquamarine", "chocolate", "cornflowerblue", "brown", "burlywood", "darkorchid", "forestgreen", "gold", "deeppink", "lawngreen", "lightgoldenrodyellow", "lightpink", "midnightblue", "deepskyblue", "gray50", "gray0")
	# set palettes names
	plotting_palette <- plotting_palette[1:length(levels(storing_all_xnas$T_level))]
	names(plotting_palette) <- levels(storing_all_xnas$T_level)

	# get a better plot ready
	# reverse levels because of coord_flip
	storing_all_xnas$T_level <- factor(storing_all_xnas$T_level, levels=rev(levels(storing_all_xnas$T_level)))

	# set DNA RNA levels
	storing_all_xnas$DNA_RNA <- factor(storing_all_xnas$DNA_RNA, levels=c("DNA", "RNA"))

	# get palette
	xna_palette <- brewer.pal(n=4, name = "Set2")[c(1,4)]
	names(xna_palette) <- c("DNA", "RNA")

	# get the plot ready
	better_plot <- storing_all_xnas %>%
				count(T_level, DNA_RNA, SE_up, SE_down, wt=Relative, name="Relative") %>%
				ggplot(aes(fill=DNA_RNA, y=Relative, x=T_level)) +
				geom_col(position = position_dodge(width = 0.9), width=0.8, color="black", lwd=1) +
				scale_fill_manual(name="DNA/RNA", values=xna_palette) +
				geom_errorbar(aes(T_level, ymin = SE_down, ymax = SE_up), position = position_dodge(width=0.9), width=0.3, lwd=2) +
				ggplot_theme(leg_pos="bottom", fnt_size=90, ang_le=0) +
				guides(fill=guide_legend(title="", override.aes=list(shape=21, size=30), ncol=2), shape=guide_legend(override.aes=list(size=30))) +
				xlab(taxa_level) +
				ylab("Relative abundance") +
				ggtitle("") +
				coord_flip() +
				theme(panel.spacing = unit(3, "lines"))

	# export treatment
	export_figs_tabs(paste0(save_taxa_figs, "BIGGER_better_plot_taxa_DNA_RNA_best_", n_best[taxa_level], "_", taxa_level, "_", exp_name, "_", clust_method, "_", taxa_algo, "_MEANS"), better_plot, as_rds=F, width=168*4, height=168*4)

	# export table used for plotting
	write.table(storing_all_xnas, paste0(save_taxa_figs, "DNA_RNA_best_", n_best[taxa_level], "_", taxa_level, "_", exp_name, "_", clust_method, "_", taxa_algo, "_MEANS.csv"), col.names=T, row.names=F, quote=F, sep="\t")
}

