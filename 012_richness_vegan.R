# source("012_richness_vegan.R")

# load libraries
library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("scales")
library("ellipse")
library("ggpubr")
library("ggsignif")
library("ggrepel")
library("ggtext")
library("gridExtra")
library("phyloseq")
library("vegan")
library("viridis")
library("data.table")
library("nlme")
library("mirlyn")

# import functions
source("/microbiology/disko2013/code/000_micro_functions_disko2013.R")

######### SELECT EXPERIMENT AND DEFINE PATHS #########
# select exp name
exp_name <- "submission_final_norejects_dada_new_bootstrap" ; orgn <- "bacteria" ; clust_method <- "ASV"

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
save_alphadiv <- paste0(save_img, "alpha_diversity/vegan/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_alphadiv), "dir exists!", dir.create(save_alphadiv, recursive=TRUE))

# load un-transformed object
phylo_data <- readRDS(paste0(save_phylo, "phylo_data_", exp_name, ".Rds"))

# get metadata
phylo_metadata <- sample_data(phylo_data)

######### INVESTIGATING RICHNESS USING VEGAN

print("RICHNESS WITH VEGAN")

# reorder Treatment and TimePoints levels to get proper colouring in the plots
sample_data(phylo_data)$Treatment <- factor(sample_data(phylo_data)$Treatment, levels=sort(unique(sample_data(phylo_data)$Treatment), decreasing=TRUE))
sample_data(phylo_data)$TimePoint <- factor(sample_data(phylo_data)$TimePoint, levels=sort(unique(sample_data(phylo_data)$TimePoint), decreasing=TRUE))

# create two phylo objects with only DNA or RNA samples
phylo_xna <- list()
phylo_xna[["DNA"]] <- prune_samples(as.character(phylo_metadata$sampleID[which(phylo_metadata$DNA_RNA == "DNA")]), phylo_data)
phylo_xna[["RNA"]] <- prune_samples(as.character(phylo_metadata$sampleID[which(phylo_metadata$DNA_RNA == "RNA")]), phylo_data)

# this is where the info will be stored
rarefaction_curve_data_summary_merged <- list()
richness_data <- list()

# rarefy or not?
rarefy <- TRUE

# empty to store everything
all_index_dna_rna <- vector()

# loop through DNA and RNA
for (xna in names(phylo_xna)) {

	# xna <- "DNA"

	# remove ASVs which are always zero for the subset under consideration
	phylo_xna[[xna]] <- prune_taxa(taxa_sums(phylo_xna[[xna]]) > 0, phylo_xna[[xna]])

	######### INVESTIGATING RICHNESS USING PHYLOSEQ

	# rarefaction curves represents the number of species as a function of the number of samples
	# the n. of samples is variable, and is used to compute the n. of species at different thresholds

	for (an_index in c("Observed", "Shannon", "Simpson")) {

		if (rarefy == TRUE) {
			print("rarefaction TRUE")
			# rarefy samples
			lib_size <- min(colSums(otu_table(phylo_xna[[xna]])))
			rarefied_xna <- mirl(phylo_xna[[xna]], libsize=lib_size, rep=100, set.seed=131, mc.cores=6L) ; rarefied <- "rarefied_"

			# compute index, using phyloseq
			richness_list <- lapply(rarefied_xna, estimate_richness, measures=an_index)
			# dataframe with mean richness for each sample
			richness_index <- rowMeans(do.call("cbind", richness_list))
			# add info to the just computed indices
			richn_tmp <- cbind.data.frame(sampleID=names(richness_index), richness_index)
		} else {
			print("rarefaction FALSE")
			# compute index, using phyloseq
			richness_index <- estimate_richness(phylo_xna[[xna]], measures=an_index) ; rarefied <- ""
			# add info to the just computed indices
			richn_tmp <- cbind.data.frame(sampleID=rownames(richness_index), richness_index)
		}

		richn_tmp$sampleID <- gsub("^X", "", richn_tmp$sampleID)
		# finalise the table construction
		richness_data[[xna]] <- merge(richn_tmp, as.matrix(sample_data(phylo_xna[[xna]])[, c("sampleID", "Treatment", "TimePoint", "CollectionSite")]), by="sampleID")

		# define palettes for plotting
		timepoint_palette <- hue_pal()(3)
		names(timepoint_palette) <- c("June", "July", "August")

		# create palette for each datapoint
		treatment_palette <- hue_pal()(2)[c(2, 1)]
		names(treatment_palette) <- c("Control", "Warming")

		# set levels
		richness_data[[xna]]$TimePoint <- factor(richness_data[[xna]]$TimePoint, levels=c("June", "July", "August"))
		richness_data[[xna]]$Treatment <- ifelse(richness_data[[xna]]$Treatment=="C", "Control", "Warming")
		richness_data[[xna]]$Treatment <- factor(richness_data[[xna]]$Treatment, levels=c("Control", "Warming"))

		# define shapes for xna. both shapes must be in the same order
		treatment_palette <- hue_pal()(2)[c(2, 1)]
		names(treatment_palette) <- c("Control", "Warming")

		# creating colours palette
		timepoint_palette <- hue_pal()(3)
		names(timepoint_palette) <- c("June", "July", "August")


		# set colnames
		colnames(richness_data[[xna]]) <- c("sampleID", "col_index", "Treatment", "TimePoint", "CollectionSite")

		# set reference levels for testing
		richness_data[[xna]]$TimePoint <- relevel(richness_data[[xna]]$TimePoint, ref="July")
		richness_data[[xna]]$Treatment <- relevel(richness_data[[xna]]$Treatment, ref="Control")

		# test significance
		an_index_model <- lme(col_index~Treatment+TimePoint, random=~1|CollectionSite, data=richness_data[[xna]])

		sink(paste0(save_alphadiv, rarefied, "vegan_lme_testing_", an_index, "_for_", xna, ".txt"))

			print(summary(an_index_model))

		sink()

		# set factors
		richness_data[[xna]]$TimePoint <- factor(richness_data[[xna]]$TimePoint, levels=c("June", "July", "August"))

		# set lims
		if (an_index == "Shannon") {
			ax_lims <- c(4.8, 6.2)
		} else if (an_index == "Simpson"){
			ax_lims <- c(0.985, 1)
		} else {
			ax_lims <- c(350, 750)
		}

		# plot by timepoint
		if (an_index == "Shannon" | an_index == "Simpson") {
			sp_time <- ggplot(richness_data[[xna]], aes(x=Treatment, y=col_index, fill=Treatment)) +
					geom_boxplot(key_glyph="point") +
					geom_point(size=2, alpha=1) +
					scale_fill_manual(values=treatment_palette) +
					ggtitle("") +
					xlab("") +
					ylim(ax_lims[1], ax_lims[2]) +
					ylab(paste0(an_index, "'s diversity")) +
					facet_wrap(~TimePoint) +
					guides(fill=guide_legend(override.aes=list(size=20, shape=22)), ncol=2) +
					ggplot_theme(leg_pos="bottom", ang_le=45) +
					theme(axis.text.x=element_blank())
		} else {
			sp_time <- ggplot(richness_data[[xna]], aes(x=Treatment, y=col_index, fill=Treatment)) +
					geom_boxplot(key_glyph="point") +
					geom_point(size=2, alpha=1) +
					scale_fill_manual(values=treatment_palette) +
					ggtitle("") +
					xlab("") +
					ylab(paste0(an_index, "'s diversity")) +
					facet_wrap(~TimePoint) +
					guides(fill=guide_legend(override.aes=list(size=20, shape=22)), ncol=2) +
					ggplot_theme(leg_pos="bottom", ang_le=45) +
					theme(axis.text.x=element_blank())
		}

		# plot boxplots
		export_svg(paste0(save_alphadiv, rarefied, an_index, "_boxplots_by_time_point_", xna, "_", exp_name), sp_time, as_rds=T)

		all_index_dna_rna <- rbind.data.frame(all_index_dna_rna, cbind.data.frame(richness_data[[xna]], xna, an_index))
	}
}

# final plot
for (an_index in c("Observed", "Shannon", "Simpson")) {

	# subset for pvals
	pvals_subset <- all_index_dna_rna[all_index_dna_rna$an_index==an_index, ]

	dna_comps <- list(c("DNA", "RNA"))
	treat_comps <- list(c("Control", "Warming"))

	tmp_time_xna <- list()
	tmp_time <- list()

	for (a_time in c("June", "July", "August")) {

		# the plot
		tmp_time_xna[[a_time]] <- ggplot(pvals_subset[pvals_subset$TimePoint==a_time, ], aes(x=xna, y=col_index, fill=Treatment)) +
					geom_boxplot(key_glyph="point") +
					scale_fill_manual(values=treatment_palette) +
					ggtitle(a_time) +
					xlab("") +
					ylab(paste0(an_index)) +
					guides(fill=guide_legend(override.aes=list(size=20, shape=22))) +
					geom_pwc(aes(group=xna), method = "wilcox.test", label = "p.signif", size=1, label.size=15, color="black", linetype=1, hide.ns=T) +
					ggplot_theme(leg_pos="bottom", ang_le=45)
	}

	for (xna in c("DNA", "RNA")) {

		# the plot
		tmp_time[[xna]] <-  ggplot(pvals_subset, aes(x=Treatment, y=col_index, fill=Treatment)) +
					geom_boxplot(key_glyph="point") +
					scale_fill_manual(values=treatment_palette) +
					ggtitle("") +
					xlab("") +
					ylab(paste0(xna, " ", an_index)) +
					guides(fill=guide_legend(override.aes=list(size=20, shape=22))) +
					geom_pwc(aes(group=Treatment), method = "wilcox.test", label = "p.signif", size=1, label.size=15, color="black", linetype=1, hide.ns=T) +
					ggplot_theme(leg_pos="bottom", ang_le=45) +
					theme(strip.background = element_rect(colour="white", fill="white")) +
					facet_wrap(~TimePoint)
	}

	# get legend
	the_legend <- cowplot::plot_grid(ggpubr::get_legend(tmp_time_xna[["June"]]))
	# plot
	sp_time_xna <- cowplot::plot_grid(
				cowplot::plot_grid(
					cowplot::plot_grid(tmp_time_xna[["June"]] + theme(legend.position = "none"),
							tmp_time_xna[["July"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
							tmp_time_xna[["August"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
							ncol=3, rel_widths=c(1.8, 1, 1)
							),
					cowplot::plot_grid(tmp_time[["DNA"]] + theme(legend.position = "none", axis.text.x=element_blank()),
							tmp_time[["RNA"]] + theme(legend.position = "none"),
							nrow=2, rel_heights=c(1,1.3)
							),
					ncol=2, rel_widths=c(1,1)
					),
			the_legend, nrow=2, rel_heights=c(1, 0.1)
			)
	# export
	export_svg(paste0(save_alphadiv, "DNA_RNA_", rarefied, an_index, "_boxplots_by_time_point_", exp_name), sp_time_xna, width=26, height=26, as_rds=T)

}

