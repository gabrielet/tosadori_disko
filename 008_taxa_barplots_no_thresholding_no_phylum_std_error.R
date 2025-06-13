# source("008_taxa_barplots_no_thresholding_no_phylum_STANDARD_ERROR.R")

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
source("/microbiology/disko2013/code/000_micro_functions_disko2013.R")

######### SELECT EXPERIMENT AND DEFINE PATHS #########

# select the proper organism
orgn <- "bacteria"

# select the clustering method
if (orgn == "bacteria") {
	clust_method <- "ASV"
	# select exp name
	exp_name <- "submission_final_norejects_dada_new_bootstrap"
} else if (orgn == "fungi") {
	clust_method <- "OTU"
	# select exp name
	exp_name <- "submission_final_norejects_dada_final_camelia"
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
save_phylo <- paste0(path_to_taxa, "phyloseq/", taxa_algo, "/")
save_taxa_figs <- paste0(save_img, "taxaFigs/")
save_genera <- paste0(save_taxa_figs, "genera_by_phyla/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_taxa_figs), "dir exists!", dir.create(save_taxa_figs, recursive=TRUE))
ifelse(dir.exists(save_genera), "dir exists!", dir.create(save_genera, recursive=TRUE))

# load non-transformed object
phylo_data <- readRDS(paste0(save_phylo, "phylo_data_", exp_name, ".Rds"))

# get metadata table
phylo_metadata <- sample_data(phylo_data)

# loop through taxa levels, to aggregate appropriately
#taxa_level <- "Phylum"
# get best n taxa for the selected taxa_level
n_best <- c(10, 15, 15)
names(n_best) <- c("Phylum", "Order", "Genus")

######### PLOTTING BEST SCORING TAXA

# create a list of phylo objects containing DNA, RNA, or ALL samples
phylo_xna <- list()
phylo_xna[["DNA"]] <- prune_samples(as.character(phylo_metadata$sampleID[which(phylo_metadata$DNA_RNA == "DNA")]), phylo_data)
phylo_xna[["RNA"]] <- prune_samples(as.character(phylo_metadata$sampleID[which(phylo_metadata$DNA_RNA == "RNA")]), phylo_data)

# loop through taxa_levels
for (taxa_level in c("Phylum", "Order", "Genus")) {

	# loop through sample types
	for (xna in names(phylo_xna)) {

		# remove ASVs which are always zero for the subset under consideration
		phylo_xna[[xna]] <- prune_taxa(taxa_sums(phylo_xna[[xna]]) > 0, phylo_xna[[xna]])

		# aggregate taxa at a specific taxonomy level
		phylo_aggregated <- tax_glom(phylo_xna[[xna]], taxrank=taxa_level)

		# get metadata for phylo_aggregated
		phylo_aggregated_meta <- sample_data(phylo_aggregated)

		# set levels to have controls first
		sample_data(phylo_aggregated)$Treatment <- factor(sample_data(phylo_aggregated)$Treatment, levels=c("C", "W"))

		###################### THIS IS THE PART WHERE THE BEST TAXA ARE RETRIEVED ######################

		# compute sums to obtain the contribution of all sampling sites
		# to a specific taxa, then order from most abundant to less abundant
		tx_sums <- sort(rowSums(otu_table(phylo_aggregated)), decreasing=T)

		# get ordered taxa names, using factors
		taxa_names <- tax_table(phylo_aggregated)[names(tx_sums), taxa_level]

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
			other_taxa <- taxa_names[(n_best[taxa_level]+1):length(taxa_names)]
		}

		################################################################################################

		# define empty to store all timepoints for each xna
		storing_all_timepoints <- vector()

		# loop through timepoints
		for (a_time in c("June", "July", "August")) {

			# subset by timepoint
			phylo_time <- prune_samples(as.character(phylo_aggregated_meta$sampleID[which(phylo_aggregated_meta$TimePoint == a_time)]), phylo_aggregated)

			# remove ASVs which are always zero for the subset under consideration
			phylo_time <- prune_taxa(taxa_sums(phylo_time) > 0, phylo_time)

			# define empty to store all treatments
			storing_all_treatments <- vector()

			print(paste0("Analysing ", xna, " in ", a_time))

			# loop through treatment
			for (a_treat in unique(sample_data(phylo_time)$Treatment)) {

				# define empty lists to store a treatment
				storing_single_treatment <- as.data.frame(matrix(ncol=7))
				# assign colnames
				colnames(storing_single_treatment) <- c("T_level", "Treatment", "TimePoint", "Relative", "SE", "SE_up", "SE_down")

				# select samples belonging to a treatment
				phylo_treat <- prune_samples(as.character(sample_data(phylo_time)$sampleID[which(sample_data(phylo_time)$Treatment == a_treat)]), phylo_time)

				# remove ASVs which are always zero for the subset under consideration
				phylo_treat <- prune_taxa(taxa_sums(phylo_treat) > 0, phylo_treat)

				# compute relative abundance
				phylo_treat <- transform_sample_counts(phylo_treat, function(x) x / sum(x))

				# loop through taxa types
				for (taxa_type in c("best_taxa", "other_taxa", "unidenti_taxa")) {

					print(paste0("Analysing ", taxa_type))

					taxa_names_sub <- ""

					if (taxa_type=="best_taxa") {
						# get best taxa
						treat_sums <- rowMeans(otu_table(phylo_treat)[which(tax_table(phylo_treat)[, taxa_level] %in% get(taxa_type)), ])
						# compute standard errors
						all_se <- apply(otu_table(phylo_treat)[which(tax_table(phylo_treat)[, taxa_level] %in% get(taxa_type)), ], 1, std.error)
						standerr_up <- treat_sums + all_se
						standerr_down <- treat_sums - all_se
						# get taxa names
						taxa_names_sub <- tax_table(phylo_treat)[names(treat_sums), taxa_level]

						# put columns together for phylo_treat samples and compute
						# the percentage of contribution each taxa is providing to
						# the total, to finally compare the taxa. by doing so, we
						# make sure the total will always sum up to 1, i.e. 100%
						storing_single_treatment <- cbind.data.frame(T_level=taxa_names_sub, Treatment=a_treat, TimePoint=a_time, Relative=treat_sums, SE=all_se, SE_up=standerr_up, SE_down=standerr_down)
						rownames(storing_single_treatment) <- NULL
						colnames(storing_single_treatment) <- c("T_level", "Treatment", "TimePoint", "Relative", "SE", "SE_up", "SE_down")

						# sort table by abundance
						storing_single_treatment <- storing_single_treatment[order(storing_single_treatment$Relative, decreasing=T), ]
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
						if (any(tax_table(phylo_treat)[, taxa_level] %in% get(taxa_type))) {
							# get taxa sums
							treat_sums <- sum(rowMeans(otu_table(phylo_treat)[which(tax_table(phylo_treat)[, taxa_level] %in% get(taxa_type)), ]))
							# compute standard errors
							all_se <- std.error(c(otu_table(phylo_treat)[which(tax_table(phylo_treat)[, taxa_level] %in% get(taxa_type)), ]))
							standerr_up <- treat_sums + all_se
							standerr_down <- treat_sums - all_se

							# put columns together for phylo_time samples and compute
							# the percentage of contribution each taxa is providing to
							# the total, to finally compare the taxa. by doing so, we
							# make sure the total will always sum up to 1, i.e. 100%
							storing_single_treatment <- cbind.data.frame(T_level=taxa_names_sub, Treatment=a_treat, TimePoint=a_time, Relative=treat_sums, SE=all_se, SE_up=standerr_up, SE_down=standerr_down)
							rownames(storing_single_treatment) <- NULL
							colnames(storing_single_treatment) <- c("T_level", "Treatment", "TimePoint", "Relative", "SE", "SE_up", "SE_down")
						} else {
							storing_single_treatment <- NULL
						}
					}
					# add to storing_all_treatments
					storing_all_treatments <- rbind.data.frame(storing_all_treatments, storing_single_treatment)
				}
			}

			print("Organising for plots")

			# remove taxa_level label
			which_level <- paste0(tolower(substr(taxa_level, 1, 1)), "__")
			storing_all_treatments$T_level <- gsub(which_level, "", storing_all_treatments$T_level)

			# set levels in order
			last_two <- as.character(c(storing_all_treatments$T_level[!grepl("Other|Unidenti", storing_all_treatments$T_level)], storing_all_treatments$T_level[grepl("Other|Unidenti", storing_all_treatments$T_level)]))
			storing_all_treatments$T_level <- factor(storing_all_treatments$T_level, levels=unique(last_two))

			# re-set labels
			storing_all_treatments$Treatment <- ifelse(storing_all_treatments$Treatment=="W", "Warming", "Control")

			# store a timepoint
			storing_all_timepoints <- rbind.data.frame(storing_all_timepoints, storing_all_treatments)
		}

		print("Final plot")

		# set levels for T_level
		last_two <- as.character(c(storing_all_timepoints$T_level[!grepl("Other|Unidenti", storing_all_timepoints$T_level)], storing_all_timepoints$T_level[grepl("Other|Unidenti", storing_all_timepoints$T_level)]))
		storing_all_timepoints$T_level <- factor(storing_all_timepoints$T_level, levels=unique(last_two))

		# set levels for TimePoint
		storing_all_timepoints$TimePoint <- factor(storing_all_timepoints$TimePoint, levels=c("June", "July", "August"))

		# set levels for Treatment
		storing_all_timepoints$Treatment <- factor(storing_all_timepoints$Treatment, levels=c("Control", "Warming"))

		# re-set palette
		plotting_palette <- c("antiquewhite", "aquamarine", "chocolate", "cornflowerblue", "brown", "burlywood", "darkorchid", "forestgreen", "gold", "deeppink", "lawngreen", "lightgoldenrodyellow", "lightpink", "midnightblue", "deepskyblue", "gray50", "gray0")
		
		# set palettes names
		plotting_palette <- plotting_palette[1:length(levels(storing_all_timepoints$T_level))]
		names(plotting_palette) <- levels(storing_all_timepoints$T_level)

		# get the plot ready
		ribon_plot <- storing_all_timepoints %>%
					count(T_level, Treatment, TimePoint, SE, wt=Relative, name="Relative") %>%
					ggplot() +
					geom_bar(aes(fill=T_level, y=Relative, x=Treatment), position="fill", stat="identity", color="black") +
					scale_fill_manual(name="", values=plotting_palette) +
					ggplot_theme(leg_pos="right", ang_le=45) +
					xlab("") +
					ylab("Relative abundance") +
					facet_wrap(~TimePoint) +
					labs(fill=taxa_level) +
					ggtitle("") +
					theme(panel.spacing = unit(3, "lines"))

		# export treatment
		export_svg(paste0(save_taxa_figs, "taxa_", xna, "_best_", n_best[taxa_level], "_", taxa_level, "_", exp_name, "_", clust_method, "_", taxa_algo, "_NEW_SUMS"), ribon_plot, as_rds=T, width=32, height=24)

		# get a better plot ready
		# reverse levels because of coord_flip
		storing_all_timepoints$T_level <- factor(storing_all_timepoints$T_level, levels=rev(levels(storing_all_timepoints$T_level)))

		# create palette
		treatment_palette <- hue_pal()(2)[c(2, 1)]
		names(treatment_palette) <- c("Control", "Warming")

		# get the plot ready
		better_plot <- storing_all_timepoints %>%
					count(T_level, Treatment, TimePoint, SE_up, SE_down, wt=Relative, name="Relative") %>%
					ggplot(aes(fill=Treatment, y=Relative, x=T_level)) +
					geom_col(position = position_dodge(width = 0.7), width=0.6, color="black", lwd=1) +
					geom_errorbar(aes(T_level, ymin = SE_down, ymax = SE_up), position = position_dodge(width=0.7), width=0.3, lwd=1) +
					scale_fill_manual(name="Treatment", values=treatment_palette) +
					ggplot_theme(leg_pos="bottom", ang_le=0) +
					guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=2), shape=guide_legend(override.aes=list(size=20))) +
					xlab("") +
					ylab("Relative abundance") +
					ggtitle("") +
					facet_wrap(~TimePoint) +
					coord_flip() +
					theme(panel.spacing = unit(3, "lines"))

		# export treatment
		export_svg(paste0(save_taxa_figs, "better_plot_taxa_", xna, "_best_", n_best[taxa_level], "_", taxa_level, "_", exp_name, "_", clust_method, "_", taxa_algo, "_NEW_SUMS"), better_plot, as_rds=T, width=32, height=24)
	}
}

