# source("007_basic_statistics_enviro.R")

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
library("scales")
library("BBmisc")
library("GGally")
library("rstatix")
library("nlme")
library("multcomp")
library("emmeans")
library("ggcorrplot")
library("plotrix")
library("cowplot")
library("viridis")

# import functions
source("/mnt/cinqueg/gabriele/work/microbiology/disko2013/code/000_micro_functions_disko2013.R")

######### SELECT EXPERIMENT AND DEFINE PATHS #########
# select exp name
exp_name <- "submission_final_norejects_dada_new_bootstrap" ; orgn <- "bacteria" ; clust_method <- "ASV"

# select the taxa assigment method
taxa_algo <- "NBC"

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, "  EXPERIMENT NAME ", exp_name))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/mnt/cinqueg/gabriele/work/microbiology/disko2013/"
# this is the general path to the experiment
path_to_exp <- paste0(root_path, orgn_dir, "experiments/", exp_name, "/")
# this is the path where counts results will be stored
path_to_taxa <- paste0(path_to_exp, clust_method, "/")
save_img <- paste0(path_to_taxa, "figures/figures_improved/", taxa_algo, "/")
save_phylo <- paste0(path_to_taxa, "phyloseq/", taxa_algo, "/")
save_box <- paste0(save_img, "boxplots/")
save_box_comparison_env <- paste0(save_box, "comparison_env/")
save_box_comparison_qpcr <- paste0(save_box, "comparison_qpcr/")
save_box_temp <- paste0(save_box, "by_timepoint/")
save_scatter <- paste0(save_img, "scatterplots/")
save_corr <- paste0(save_img, "correlograms/")
# this is the path where counts results will be stored
save_final_figs <- paste0(save_img, "final_figs/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_box), "dir exists!", dir.create(save_box, recursive=TRUE))
ifelse(dir.exists(save_box_comparison_env), "dir exists!", dir.create(save_box_comparison_env, recursive=TRUE))
ifelse(dir.exists(save_box_comparison_qpcr), "dir exists!", dir.create(save_box_comparison_qpcr, recursive=TRUE))
ifelse(dir.exists(save_box_temp), "dir exists!", dir.create(save_box_temp, recursive=TRUE))
ifelse(dir.exists(save_scatter), "dir exists!", dir.create(save_scatter, recursive=TRUE))
ifelse(dir.exists(save_corr), "dir exists!", dir.create(save_corr, recursive=TRUE))
ifelse(!dir.exists(save_final_figs), dir.create(save_final_figs, recursive=TRUE), "dir exists!")

# load non-transformed object
phylo_data <- readRDS(paste0(save_phylo, "phylo_data_", exp_name, ".Rds"))

# load metadata table
phylo_metadata <- sample_data(phylo_data)

# get DNA only bacause metadata are repeted between RNA and DNA
phylo_xna <- prune_samples(as.character(phylo_metadata$sampleID[which(phylo_metadata$DNA_RNA == "DNA")]), phylo_data)

# get metadata table as dataframe
meta_phylo_xna <- as(sample_data(phylo_xna), "data.frame")

# set levels and references
meta_phylo_xna$TimePoint <- factor(meta_phylo_xna$TimePoint, levels=c("June", "July", "August"))
meta_phylo_xna$Treatment <- ifelse(meta_phylo_xna$Treatment=="C", "Control", "Warming")
meta_phylo_xna$Treatment <- factor(meta_phylo_xna$Treatment, levels=c("Control", "Warming"))
meta_phylo_xna$CollectionSite <- factor(meta_phylo_xna$CollectionSite, levels=unique(meta_phylo_xna$CollectionSite))

# rename cols for plotting
colnames(meta_phylo_xna) <- c("sampleID", "DNA_RNA", "Treatment", "TimePoint", "CollectionSite", "SM", "SM_perc", "pH", "NH4", "NH4_SOM", "NO3", "NO3_SOM", "PO4", "PO4_SOM", "N", "C", "mcrA", "mxaF", "nifH", "nirS", "nosZ", "TOC", "TOC_SOM", "microC", "microC_SOM", "DON", "DON_SOM", "TN", "TN_SOM", "microN", "microN_SOM", "microCN", "d15N", "d13C", "SOM", "LOI", "Deciduous_shrub", "Evergreen_shrub", "Lichen", "Forb", "Litter", "plant_diversity", "plant_eveness", "soil_fresh", "soil_dry", "soil_water", "CN")

# load centered log ratio qpcr object
phylo_qpcr <- readRDS(paste0(save_phylo, "phylo_qpcr_", exp_name, ".Rds"))

# get qpcr metadata
meta_phylo_qpcr <- sample_data(phylo_qpcr)

# rename cols for plotting
colnames(meta_phylo_qpcr) <- c("sampleID", "DNA_RNA", "Treatment", "TimePoint", "CollectionSite", "mcrA", "mxaF", "nifH", "nirS", "nosZ", "soil_fresh", "soil_dry")

################################################# PERFORM SIMPLE TEST ON THE ABUNDANCES

# get sample names
control_samples <- meta_phylo_xna[meta_phylo_xna$Treatment=="Control", ]$sampleID
warming_samples <- meta_phylo_xna[meta_phylo_xna$Treatment=="Warming", ]$sampleID

# get counts for each group of samples
control_counts <- c(otu_table(phylo_xna)[, which(colnames(otu_table(phylo_xna)) %in% control_samples)])
warming_counts <- c(otu_table(phylo_xna)[, which(colnames(otu_table(phylo_xna)) %in% warming_samples)])

# perform test to compare counts
t.test(control_counts, warming_counts, paired=T)

################################################# PLOT SCATTERPLOTS

# select all params
physico_par <- c("pH", "SM", "C", "N", "CN", "SOM", "LOI", "TOC", "TN", "DON", "NH4", "NO3", "PO4", "d15N", "d13C", "microC", "microN", "microCN")
qpcr_par <- c("nirS", "nosZ", "nifH", "mcrA", "mxaF")
# put params together
all_available_params <- c(physico_par, qpcr_par)

##################################################################### MAKE TABLE FOR PAPER
# export table
selected_meta <- t(meta_phylo_xna[, all_available_params])

# list with new columns
new_cols <- list()

for (a_time in c("June", "July", "August")) {
	for (a_treat in c("C", "W")) {

#		a_time <- "June" ; a_treat <- "C"
		# subset by time
		a_time_cols <- selected_meta[, grep(a_time, colnames(selected_meta))]
		# subset time by treatment
		a_treat_cols <- a_time_cols[, grep(paste0("T", a_treat), colnames(a_time_cols))]
		# modify label
		treat_lbl <- ifelse(a_treat=="C", "Control", "Warming")
		# compute mean and std.error
		new_cols[[paste0(a_time, "x", treat_lbl)]] <- paste0("$", round(apply(a_treat_cols, 1, mean, na.rm=T), digits=2), " \\pm ", round(apply(a_treat_cols, 1, std.error, na.rm=T), digits=2), "$")
	}
}

# put all together
final_meta <- cbind.data.frame(parameter=rownames(selected_meta), do.call(cbind.data.frame, new_cols))

# export one for overview
write.table(final_meta, paste0(save_scatter, "final_metadata_table.txt"), row.names=F, col.names=T, quote=F, sep=" & ")
# and one with all metadata
write.table(selected_meta, paste0(save_scatter, "final_metadata_table.tsv"), row.names=F, col.names=T, quote=F, sep="\t")
#####################################################################

# create a matrix to store the info
# one col for each parameter
# one row for each combination of C-W and June-July-August
paper_meta <- as.data.frame(matrix(ncol=length(all_available_params), nrow=6))
rownames(paper_meta) <- c("Control-June", "Control-July", "Control-August", "Warming-June", "Warming-July", "Warming-August")
colnames(paper_meta) <- all_available_params

# sort metadata for paper
for (timepoint in unique(meta_phylo_xna$TimePoint)) {
	for (treatment in unique(meta_phylo_xna$Treatment)) {
		# current analysis
		current <- paste0(treatment, "-", timepoint)
		phylo_sub <- meta_phylo_xna[which(meta_phylo_xna$Treatment==treatment & meta_phylo_xna$TimePoint==timepoint), all_available_params]
		paper_meta[current, ] <- colMeans(phylo_sub, na.rm = TRUE)
	}
}

# finalise table and export
paper_meta_final <- cbind.data.frame(rownames(paper_meta), paper_meta)
colnames(paper_meta_final) <- c("Samples", colnames(paper_meta))

write.table(paper_meta_final, paste0(save_scatter, "params_summary.csv"), row.names=F, col.names=T, sep="\t")

print("PLOT SCATTERPLOTS")

# list interesting variables
interesting_params <- list() ; plot_cols <- list()
interesting_params[["qpcr"]] <- qpcr_par ; plot_cols[["qpcr"]] <- 5
# this is final with parameters from RDA
interesting_params[["final"]] <- c("pH", "C", "N", "SOM", "LOI", "TOC", "TN", "DON", "PO4", "d15N", "microC", "microN", "microCN") ; plot_cols[["final"]] <- 8

# select parameters
plot_pars <- c("qpcr", "final")

# standardise data, to avoid issues with different magnitudes
# and measure units between environmental parameters
normalised_env <- apply(meta_phylo_xna[, c(6:ncol(meta_phylo_xna))], 2, normalize, method="range", range=c(0, 1))
colnames(normalised_env) <- paste0("norm_", colnames(normalised_env))
standard_phylo <- data.frame(sampleID=meta_phylo_xna$sampleID, TimePoint=meta_phylo_xna$TimePoint, Treatment=meta_phylo_xna$Treatment, CollectionSite=meta_phylo_xna$CollectionSite, normalised_env, meta_phylo_xna[, c(6:ncol(meta_phylo_xna))])

# plot scatterplot using the subset of parameters defined
# https://fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/
for (intp in plot_pars) {

#	intp <- "final"

	# printing out what's on
	print(paste0("plotting ", paste(interesting_params[[intp]], collapse=", ")))

	# remove nas
	no_nas <- standard_phylo[complete.cases(standard_phylo), ]

	scatter_time <- ggplot_pairs(no_nas, interesting_params[[intp]], "TimePoint", n_cols=plot_cols[[intp]], normalised=T)

	export_figs_tabs(paste0(save_scatter, "standard_phylo_by_timepoint_", intp), scatter_time, width=168*9, height=168*12, base=FALSE)

	scatter_treat <- ggplot_pairs(no_nas, interesting_params[[intp]], "Treatment", n_cols=plot_cols[[intp]], normalised=T)

	export_figs_tabs(paste0(save_scatter, "standard_phylo_by_treatment_", intp), scatter_treat, width=168*9, height=168*12, base=FALSE)
}

################################################# PLOTTING CORRELOGRAMS

# correlograms are redundant with scatterplots, but provide
# a different visual style. this one is more easy to read
# than the plot above.

print("PLOT CORRELOGRAMS")

# test correlations
pvals_corr <- vector()

# plot correlations
corr_plot <- ggcorr(meta_phylo_xna[, all_available_params], label=F, size=10, hjust=0.78)

# for each parameter
for (a_par in 1:nrow(corr_plot$data)) {
	# get params names
	par_i <- corr_plot$data$x[a_par]
	par_j <- as.character(corr_plot$data$y[a_par])
	# compute the correlation and get the p-value
	pvals_corr <- append(pvals_corr, cor.test(meta_phylo_xna[, par_i], meta_phylo_xna[, par_j], alternative="two.sided", method = "pearson")$p.value)
}
# add pvalues
corr_plot$data$pvalues <- format(round(pvals_corr, digits=2), nsmall=2)

corr_plot$data$pvalues[corr_plot$data$pvalues > 0.05] <- ""

final_corr <- corr_plot + geom_text(aes(label=corr_plot$data$pvalues), size=10) + ggplot_theme(leg_pos="right") + guides(fill = guide_colorbar(barheight = 20, barwidth = 2)) + theme(legend.text=element_text(size=30))

# export
export_figs_tabs(paste0(save_corr, "new_correlogram_pearson_", exp_name, "_easy_to_read"), final_corr, base=FALSE)

################################################# TESTING SOME VARIABLES WITH LINEAR MIXED EFFECT MODELS AND THE MULTCOMP PACKAGE

print("TESTING SOME VARIABLES NO INTERACTIONS")

# set seed
set.seed(131)

# test linear
pvals_lin <- vector()
pvals_lin_to_export <- vector()

# for each enviro variable
for (an_env in all_available_params){

#	an_env <- all_available_params[1]

	# subset by param
	subset <- meta_phylo_xna[, c("Treatment", "TimePoint", "CollectionSite", an_env)]
	colnames(subset) <- c("Treatment", "TimePoint", "CollectionSite", "an_env")

	# set reference levels
	subset$TimePoint <- relevel(subset$TimePoint, ref="July")
	subset$Treatment <- relevel(subset$Treatment, ref="Control")

	# perform the test
	if (orgn == "bacteria") {
		an_env_model <- lme(an_env~Treatment+TimePoint, random=~1|CollectionSite, data=subset, na.action=na.omit)
	} else {
		an_env_model <- lme(an_env~Treatment, random=~1|CollectionSite, data=subset, na.action=na.omit)
	}

	# this is for plotting
	pvals_lin <- rbind.data.frame(pvals_lin, cbind.data.frame(an_env, rbind(round(summary(an_env_model)$tTable[, 5], digits=4))))
	pvals_lin_to_export <- rbind.data.frame(pvals_lin_to_export, cbind.data.frame(an_env, param=rownames(summary(an_env_model)$tTable), rbind(round(summary(an_env_model)$tTable, digits=4))))

}
# remove intercept
pvals_lin$parameter <- gsub("[0-9]", "", rownames(pvals_lin))
pvals_lin <- pvals_lin[!grepl("Intercept", pvals_lin$parameter), ]

# export
write.table(pvals_lin_to_export, paste0(save_box, "all_enviro_stats_no_interaction_refs_C_July_NEW_CORRECTED_BUT_OLD.csv"), row.names=F, quote=F, sep="\t")

################################################# PLOTTING BOXPLOTS WARMING VERSUS CONTROLS

print("PLOTTING BOXPLOTS FOR ENVIRO")
print("this is using all parameters")

# put together pvalues for the boxplots
final_all_pvals <- cbind.data.frame(pvals_lin, treat_sym=ifelse(0 <= pvals_lin$TreatmentWarming & pvals_lin$TreatmentWarming <=0.001, "***", ifelse(0.001 < pvals_lin$TreatmentWarming & pvals_lin$TreatmentWarming <= 0.01, "**", ifelse(0.01< pvals_lin$TreatmentWarming & pvals_lin$TreatmentWarming <= 0.05, "*", ifelse(0.05 < pvals_lin$TreatmentWarming & pvals_lin$TreatmentWarming <= 0.1, "°", "")))), june_sym=ifelse(0 <= pvals_lin$TimePointJune & pvals_lin$TimePointJune <=0.001, "***", ifelse(0.001 < pvals_lin$TimePointJune & pvals_lin$TimePointJune <= 0.01, "**", ifelse(0.01< pvals_lin$TimePointJune & pvals_lin$TimePointJune <= 0.05, "*", ifelse(0.05 < pvals_lin$TimePointJune & pvals_lin$TimePointJune <= 0.1, "°", "")))), aug_sym=ifelse(0 <= pvals_lin$TimePointAugust & pvals_lin$TimePointAugust <=0.001, "***", ifelse(0.001 < pvals_lin$TimePointAugust & pvals_lin$TimePointAugust <= 0.01, "**", ifelse(0.01< pvals_lin$TimePointAugust & pvals_lin$TimePointAugust <= 0.05, "*", ifelse(0.05 < pvals_lin$TimePointAugust & pvals_lin$TimePointAugust <= 0.1, "°", "")))))

# empty var for putting all the data together, for plotting
long_tab_env <- vector()

# reshape table for plotting with facets
for (a_var in all_available_params){

#	a_var <- "microN"

	# get a_var sub table
	var_tab <- cbind.data.frame(standard_phylo[, c(c("sampleID", "TimePoint", "Treatment", "CollectionSite", paste0("norm_", a_var), a_var))], variable=rep(a_var, nrow(standard_phylo)), symbol=rep(final_all_pvals$treat_sym[final_all_pvals$an_env==a_var], nrow(standard_phylo)))
	# rename cols to make them standard for all the variable sub tables
	colnames(var_tab) <- c("sampleID", "TimePoint", "Treatment", "CollectionSite", "norm_value", "value", "variable", "symbol")

	# finally, put all together
	long_tab_env <- rbind.data.frame(long_tab_env, var_tab)
}

# merge variable and symbol
long_tab_env$stat_lab <- paste0(long_tab_env$variable, " ", long_tab_env$symbol)

# print with paired version, Warming versus Control
# add plants to get boxplots
pb <- print_boxplots_paired_bacteria(long_tab_env, all_available_params, exp_name, save_box_comparison_env)

# create palette
treatment_palette <- hue_pal()(2)[c(2, 1)]
names(treatment_palette) <- c("Control", "Warming")

# set levels for plotting
long_tab_env$variable <- factor(long_tab_env$variable, levels=all_available_params)

# empty list for plotting
multiplot <- list()

# reshape pvals table
final_all_pvals_reshaped <- cbind.data.frame(var=rep(final_all_pvals$an_env, 2), TimePoint=c(rep("June", nrow(final_all_pvals)), rep("August", nrow(final_all_pvals))), p.signif=c(final_all_pvals$june_sym, final_all_pvals$aug_sym), test=c(final_all_pvals$TimePointJune, final_all_pvals$TimePointAugust), group1=c(rep("June", nrow(final_all_pvals)), rep("August", nrow(final_all_pvals))), group2="July")

# remove all rows where the symbol is not present
only_valid_pvals <- final_all_pvals_reshaped[final_all_pvals_reshaped$p.signif!="", ]

# multiplot
for (an_env in all_available_params) {

	# subset for a parameter
	subset <- long_tab_env[long_tab_env$variable==an_env, ]
	# subset for pvals
	pvals_subset <- only_valid_pvals[only_valid_pvals$var==an_env, ]
	# plot normalised
	dodge <- position_dodge(width = 0.8)
	tmp_multiplot <- ggplot(subset[!is.na(subset$norm_value), ], aes(x=TimePoint, y=norm_value, fill=Treatment)) +
			geom_violin(width=1, position=dodge, linewidth=1) +
			geom_boxplot(width=0.2, color="black", position=dodge, linewidth=1, show.legend=FALSE) +
			geom_point(size=2, position=dodge, show.legend=F) +
			scale_fill_manual(name="Treatment", values=treatment_palette) +
			xlab("") +
			ylim(0, 1.3) +
			ylab("Normalised unit") +
			ggtitle(subset$stat_lab) +
			guides(fill=guide_legend(override.aes=list(shape=22, size=20), ncol=1)) +
			ggplot_theme(leg_pos="bottom", ang_le=60, fnt_size=70)
	if (nrow(pvals_subset) > 0) {
		multiplot[[an_env]] <- tmp_multiplot + stat_pvalue_manual(pvals_subset, label="p.signif", y.position=1.1, inherit.aes=F, size=30, bracket.size=2.5, color="black", linetype=1, step.increase=0.1)
	} else {
		multiplot[[an_env]] <- tmp_multiplot
	}
}

# get legend
the_legend <- cowplot::plot_grid(ggpubr::get_legend(multiplot[["pH"]]))

# enviro
all_together <- cowplot::plot_grid(
			cowplot::plot_grid(
					multiplot[["pH"]] + theme(legend.position = "none", axis.text.x=element_blank()),
					multiplot[["SM"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["C"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["N"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["CN"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					the_legend,
					ncol=6, rel_widths=c(1.3,1,1,1,1,1)
				),
			cowplot::plot_grid(
					multiplot[["SOM"]] + theme(legend.position = "none", axis.text.x=element_blank()),
					multiplot[["LOI"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["TOC"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["TN"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["DON"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["NH4"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["NO3"]] + theme(legend.position = "none", axis.text.x=element_blank()),
					multiplot[["PO4"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["d15N"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["d13C"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["microC"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["microN"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["microCN"]] + theme(legend.position = "none"),
					multiplot[["nirS"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["nosZ"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["nifH"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["mcrA"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["mxaF"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					nrow=3, ncol=6, rel_widths=c(1.3,1,1,1,1,1), rel_heights=c(1,1,1.2)
				),
				nrow=2, rel_heights=c(1,3.3)
		)

# export
export_figs_tabs(paste0(save_box_comparison_env, "ULTRA_NEW_BUT_OLD_all_together_", exp_name), all_together, height=168*9, width=168*9)

# enviro
all_together_no_genes <- cowplot::plot_grid(
					multiplot[["pH"]] + theme(legend.position = "none", axis.text.x=element_blank()),
					multiplot[["SM"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["C"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["N"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["CN"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["SOM"]] + theme(legend.position = "none", axis.text.x=element_blank()),
					multiplot[["LOI"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["TOC"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["TN"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["DON"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["NH4"]] + theme(legend.position = "none", axis.text.x=element_blank()),
					multiplot[["NO3"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["PO4"]] + theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["d15N"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["d13C"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["microC"]] + theme(legend.position = "none"),
					multiplot[["microN"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					multiplot[["microCN"]] + theme(legend.position = "none", axis.text.y=element_blank(), axis.title.y=element_blank()),
					the_legend, align="vh", axis = "tblr"
				)

# export
export_figs_tabs(paste0(save_box_comparison_env, "no_gene_all_together_", exp_name), all_together_no_genes, height=168*3, width=168*3)

# plot with paired, Warming or Control at all timepoints
for (a_treat in unique(long_tab_env$Treatment)) {

#	a_treat <- "Control" ; an_env <- "pH"

	print(paste0("Analysing ", a_treat))

	# subset by treatment
	env_by_treat <- long_tab_env[which(long_tab_env$Treatment==a_treat), ]
	# plot normalised
	dodge <- position_dodge(width = 0.8)
	plt_all <- ggplot(env_by_treat[which(env_by_treat$variable %in% all_available_params), ], aes(x=TimePoint, y=norm_value, fill=TimePoint)) +
		geom_violin(width=1, position=dodge, linewidth=1) +
		geom_boxplot(width=0.2, color="black", position=dodge, linewidth=1, show.legend=FALSE) +
		geom_point(size=2, position=dodge, show.legend=FALSE) +
		xlab("") +
		ylim(c(0, 1.3)) +
		ylab("Normalised unit") +
		ggtitle("") +
		guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=2), shape=guide_legend(override.aes=list(size=20))) +
		facet_wrap(~variable) +
		ggplot_theme(leg_pos="none", ang_le=45)

	# export
	export_figs_tabs(paste0(save_box_comparison_env, "all_enviro_comparison_by_time_", a_treat, "_samples_", exp_name), plt_all)
}

################################################# PLOTTING BOXPLOTS FOR qPCR GENES

print("PLOT BOXPLOTS FOR QPCR")

# get qpcr parameters only
# normalise them for better plotting
normalised_qpcr <- apply(meta_phylo_qpcr[, c("nirS", "nosZ", "nifH", "mcrA", "mxaF")], 2, normalize, method="range", range=c(0, 1))
colnames(normalised_qpcr) <- paste0("norm_", colnames(normalised_qpcr))
standard_qpcr <- data.frame(sampleID=meta_phylo_qpcr$sampleID, TimePoint=meta_phylo_qpcr$TimePoint, Treatment=meta_phylo_qpcr$Treatment, CollectionSite=meta_phylo_qpcr$CollectionSite, normalised_qpcr, meta_phylo_qpcr[, c("nirS", "nosZ", "nifH", "mcrA", "mxaF")])

# set levels for timepoint and treatment
standard_qpcr$TimePoint <- factor(standard_qpcr$TimePoint, levels=c("June", "July", "August"))
standard_qpcr$Treatment <- ifelse(standard_qpcr$Treatment=="C", "Control", "Warming")
standard_qpcr$Treatment <- factor(standard_qpcr$Treatment, levels=c("Control", "Warming"))

# empty var for putting all the data together, for plotting
long_tab_qpcr <- vector()

# reshape table for plotting with facets
for (a_gene in qpcr_par){

#	a_gene <- "mxaF"

	# get gene sub table
	gene_tab <- cbind.data.frame(standard_qpcr[, c(c("sampleID", "TimePoint", "Treatment", "CollectionSite", paste0("norm_", a_gene), a_gene))], gene=rep(a_gene, nrow(standard_qpcr)))
	# rename cols to make them standard for all the genes sub tables
	colnames(gene_tab) <- c("sampleID", "TimePoint", "Treatment", "CollectionSite", "norm_value", "value", "variable")

	# merge together
	long_tab_qpcr <- rbind.data.frame(long_tab_qpcr, gene_tab)
}

# print with paired version, Warming versus Control
pb <- print_boxplots_paired(long_tab_qpcr, qpcr_par, exp_name, save_box_comparison_qpcr)

# create palette
treatment_palette <- hue_pal()(2)[c(2, 1)]
names(treatment_palette) <- c("Control", "Warming")

# plot normalised
dodge <- position_dodge(width = 0.8)
trt_pl_qpcr <- ggplot(long_tab_qpcr, aes(x=TimePoint, y=norm_value, fill=Treatment)) +
	geom_violin(width=1, position=dodge, linewidth=1) +
	geom_boxplot(width=0.2, color="black", position=dodge, linewidth=1, show.legend=FALSE) +
	geom_point(size=2, position=dodge, show.legend=FALSE) +
	scale_fill_manual(name="Treatment", values=treatment_palette) +
	xlab("") +
	ylab("Normalised unit") +
	ggtitle("") +
	guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=2)) +
	facet_wrap(~variable, ncol=5) +
	ggplot_theme(leg_pos="bottom", ang_le=45)

# export
export_figs_tabs(paste0(save_box_comparison_qpcr, "all_genes_comparison_by_treatment_samples_", exp_name), trt_pl_qpcr, width=168*2, height=168*1.5)

# plot with paired, Warming or Control at all timepoints
for (a_treat in unique(long_tab_qpcr$Treatment)) {

#	a_treat <- "Control" ; a_gene <- "mcrA"

	print(paste0("Analysing ", a_treat))

	# subset by treatment
	qpcr_by_treat <- long_tab_qpcr[which(long_tab_qpcr$Treatment==a_treat), ]

	# plot normalised
	dodge <- position_dodge(width = 0.8)
	plt_all <- ggplot(qpcr_by_treat, aes(x=TimePoint, y=norm_value, fill=TimePoint)) +
		geom_violin(width=1, position=dodge, linewidth=1) +
		geom_boxplot(width=0.2, color="black", position=dodge, linewidth=1, show.legend=FALSE) +
		geom_point(size=2, position=dodge, show.legend=FALSE) +
		xlab("") +
		ylab("Normalised unit") +
		guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=2), shape=guide_legend(override.aes=list(size=20))) +
		facet_wrap(~variable, ncol=3) +
		ggplot_theme(leg_pos="none", ang_le=45)

	# export
	export_figs_tabs(paste0(save_box_comparison_qpcr, "all_genes_comparison_by_time_", a_treat, "_samples_", exp_name), plt_all)

}

################################################# SURFACE TEMP DATA
if (TRUE) {
	print("PLOTTING SURFACE TEMPERATURE DATA")

	# load surface temperature data
	surf_temp_data <- read.csv("/mnt/cinqueg/gabriele/work/microbiology/disko2013/metadata/csv_original/surface_temp_analysed_timeframe.csv", header=T, sep=",")

	# set date as date, in proper format
	surf_temp_data$day <- as.Date(surf_temp_data$day, format="%m-%d-%Y")

	# discard july 2012
	surf_temp_data <- subset(surf_temp_data, day>="2012-08-01")

	info_mean_surf <- vector()

	# loop through the days
	for (day in unique(surf_temp_data$day)) {
		# loop through treatments
		for (treat in unique(surf_temp_data$treatment)) {
			# compute average temperature
			mean_temp <- mean(surf_temp_data$t_surf[which(surf_temp_data$day==day & surf_temp_data$treatment==treat)][!is.na(surf_temp_data$t_surf[which(surf_temp_data$day==day & surf_temp_data$treatment==treat)])])
			# get info for a single day
			info_mean_surf <- rbind.data.frame(info_mean_surf, cbind.data.frame(day=unique(surf_temp_data$day[which(surf_temp_data$day==day & surf_temp_data$treatment==treat)]), treatment=treat, mean_temp))
		}
	}

	# set levels
	info_mean_surf$treatment <- ifelse(info_mean_surf$treatment=="C", "Control", "Warming")
	info_mean_surf$treatment <- factor(info_mean_surf$treatment, levels=c("Control", "Warming"))

	# create palette for each datapoint
	treatment_palette <- hue_pal()(2)[c(2, 1)]
	names(treatment_palette) <- c("Control", "Warming")

	# plot temperatures timeline
	timeline_surf <- ggplot(info_mean_surf, aes(x=day, y=mean_temp, color=treatment)) +
			geom_line(lwd=2.5, key_glyph="point") +
			scale_color_manual(values=treatment_palette) +
			scale_x_date("", date_breaks = "1 months") +
			scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 5)) +
			ylab("Surface Temperature, °C") +
			guides(color=guide_legend(override.aes=list(size=20, shape=15), title=""), ncol=2) +
			ggplot_theme(leg_pos="right", ang_le=45)

	# export
	export_figs_tabs(paste0(save_box_temp, "surface_yearly_temperatures_time_line_", exp_name), timeline_surf)

	# compute stats for temperatures
	by_month_surf <- as.data.frame(surf_temp_data %>% group_by(month = lubridate::floor_date(surf_temp_data$day, 'month')))
	by_month_surf$treatment <- ifelse(by_month_surf$treatment=="C", "Control", "Warming")
	by_month_surf$treatment <- factor(by_month_surf$treatment, levels=c("Control", "Warming"))

	# add column for aggregating temps
	by_month_surf$month_by_treat <- paste(by_month_surf$treatment, by_month_surf$month, sep="_")

	# get month only
	by_month_surf$string_month <- zoo::as.yearmon(format(by_month_surf$day, "%Y-%m"))

	# compute mean and std.error
	all_means <- vector()
	all_ses <- vector()
	all_stdvs <- vector()
	months <- vector()
	for (aggreg in unique(by_month_surf$month_by_treat)) {
		# get subset
		subset <- by_month_surf[by_month_surf$month_by_treat == aggreg, ]
		months <- append(months, aggreg)
		all_means <- append(all_means, mean(subset$t_surf))
		all_ses <- append(all_ses, std.error(subset$t_surf))
		all_stdvs <- append(all_stdvs, sd(subset$t_surf))
	}
	final_data_surf <- cbind.data.frame(months, avg=all_means, stdev=all_stdvs, se=all_ses)

	# create table to export aggregated temperatures
	write.table(final_data_surf, paste0(save_box_temp, "aggregated_surface_temps.csv"), col.names=T, row.names=F, quote=F, sep="\t")

	# compute stats to add to temperature boxplots
	stat_test <- by_month_surf %>%
		group_by(string_month) %>%
		pairwise_t_test(t_surf ~ treatment, ref.group="Control") %>%
		adjust_pvalue(method = "fdr") %>%
		add_significance()

	stat_test <- stat_test %>% add_xy_position(x="treatment")
	stat_test$p.adj <- format(round(stat_test$p.adj, digits=3), nsmall = 3)

	# create boxplot
	dodge <- position_dodge(width = 0.8)
	bxp_surf <- ggplot(by_month_surf, aes(x=treatment, y=t_surf)) +
			geom_boxplot(aes(fill=treatment), color="black", position=dodge, linewidth=1, key_glyph = "point") +
			geom_point(size=2, position=dodge) +
			scale_fill_manual(name="", values=treatment_palette) +
			xlab("") +
			ylim(c(-20, 25)) +
			ylab("Surface Temperature, °C") +
			facet_wrap(~string_month, ncol=9) +
			guides(fill=guide_legend(override.aes=list(shape=22, size=20), ncol=2), shape=guide_legend(override.aes=list(size=20))) +
			ggplot_theme(leg_pos="bottom", ang_le=45) +
			theme(axis.text.x=element_blank())

	# add stats to faceted boxplots
	temp_box_surf <- bxp_surf + stat_pvalue_manual(stat_test, label="p.adj", hide.ns = TRUE, size=10)

	# export
	export_figs_tabs(paste0(save_box_temp, "surface_yearly_temperatures_box_", exp_name), temp_box_surf, width=168*2, height=168*1.5)

	# define summer
	summer_season <- c("2013-06-01", "2013-07-01", "2013-08-01")

	# define empty
	surf_differences <- vector()

	# loop through summer
	for (a_month in summer_season) {

		# subset by month and treatment
		# and compute average temp
		summer_surf_warming <- mean(by_month_surf$t_surf[by_month_surf$treatment=="Warming" & by_month_surf$string_month==a_month][!is.na(by_month_surf$t_surf[by_month_surf$treatment=="Warming" & by_month_surf$string_month==a_month])])
		summer_surf_control <- mean(by_month_surf$t_surf[by_month_surf$treatment=="Control" & by_month_surf$string_month==a_month][!is.na(by_month_surf$t_surf[by_month_surf$treatment=="Control" & by_month_surf$string_month==a_month])])
		surf_differences <- rbind.data.frame(surf_differences, cbind.data.frame(a_month, diff=summer_surf_warming-summer_surf_control))

	}

################################################# 5 CM DEEP TEMP DATA

	print("PLOTTING FIVE DEEP TEMPERATURE DATA")

	# load fivedace temperature data
	five_deep_temp <- read.csv("/mnt/cinqueg/gabriele/work/microbiology/disko2013/metadata/csv_original/five_cm_deep_temp_analysed_timeframe.csv", header=T, sep=",")

	# set date as date, in proper format
	five_deep_temp$day <- as.Date(five_deep_temp$day, format="%m-%d-%Y")

	# discard july 2012
	five_deep_temp <- subset(five_deep_temp, day>="2012-08-01")

	info_mean_fived <- vector()

	# loop through the days
	for (day in unique(five_deep_temp$day)) {
		# loop through treatments
		for (treat in unique(five_deep_temp$treatment)) {
			# compute average temperature
			mean_temp <- mean(five_deep_temp$five_deep[which(five_deep_temp$day==day & five_deep_temp$treatment==treat)][!is.na(five_deep_temp$five_deep[which(five_deep_temp$day==day & five_deep_temp$treatment==treat)])])
			# get info for a single day
			info_mean_fived <- rbind.data.frame(info_mean_fived, cbind.data.frame(day=unique(five_deep_temp$day[which(five_deep_temp$day==day & five_deep_temp$treatment==treat)]), treatment=treat, mean_temp))
		}
	}

	# set levels
	info_mean_fived$treatment <- ifelse(info_mean_fived$treatment=="C", "Control", "Warming")
	info_mean_fived$treatment <- factor(info_mean_fived$treatment, levels=c("Control", "Warming"))

	# create palette for each datapoint
	treatment_palette <- hue_pal()(2)[c(2, 1)]
	names(treatment_palette) <- c("Control", "Warming")

	# plot temperatures timeline
	dodge <- position_dodge(width = 0.8)
	timeline_five <- ggplot(info_mean_fived, aes(x=day, y=mean_temp, color=treatment)) +
			geom_line(lwd=2.5, key_glyph="point") +
			scale_color_manual(values=treatment_palette) +
			scale_x_date("", date_breaks = "1 months") +
			scale_y_continuous(limits = c(-15, 20), breaks = seq(-10, 20, by = 5)) +
			ylab("Topsoil Temperature, °C") +
			guides(color=guide_legend(override.aes=list(size=20, shape=15), title=""), ncol=2) +
			ggplot_theme(leg_pos="right", ang_le=45)

	# export
	export_figs_tabs(paste0(save_box_temp, "five_deep_yearly_temperatures_time_line_", exp_name), timeline_five)

	# compute stats for temperatures
	by_month_fived <- as.data.frame(five_deep_temp %>% group_by(month = lubridate::floor_date(five_deep_temp$day, 'month')))
	by_month_fived$treatment <- ifelse(by_month_fived$treatment=="C", "Control", "Warming")
	by_month_fived$treatment <- factor(by_month_fived$treatment, levels=c("Control", "Warming"))

	# add column for aggregating temps
	by_month_fived$month_by_treat <- paste(by_month_fived$treatment, by_month_fived$month, sep="_")

	# get month only
	by_month_fived$string_month <- zoo::as.yearmon(format(by_month_fived$day, "%Y-%m"))

	# compute mean and std.error
	all_means <- vector()
	all_ses <- vector()
	all_stdvs <- vector()
	months <- vector()
	for (aggreg in unique(by_month_fived$month_by_treat)) {
		# get subset
		subset <- by_month_fived[by_month_fived$month_by_treat == aggreg, ]
		months <- append(months, aggreg)
		all_means <- append(all_means, mean(subset$five_deep))
		all_ses <- append(all_ses, std.error(subset$five_deep))
		all_stdvs <- append(all_stdvs, sd(subset$five_deep))
	}
	final_data_five_d <- cbind.data.frame(months, avg=all_means, stdev=all_stdvs, se=all_ses)

	# create table to export aggregated temperatures
	write.table(final_data_five_d, paste0(save_box_temp, "aggregated_five_deep_temps.csv"), col.names=T, row.names=F, quote=F, sep="\t")

	# compute stats to add to temperature boxplots
	stat_test <- by_month_fived %>%
		group_by(string_month) %>%
		pairwise_t_test(five_deep ~ treatment, ref.group="Control") %>%
		adjust_pvalue(method = "fdr") %>%
		add_significance()

	stat_test <- stat_test %>% add_xy_position(x="treatment")
	stat_test$p.adj <- format(round(stat_test$p.adj, digits=3), nsmall = 3)

	# create boxplot
	bxp_fived <- ggplot(by_month_fived, aes(x=treatment, y=five_deep)) +
			geom_boxplot(aes(fill=treatment), color="black", position=dodge, linewidth=1, key_glyph = "point") +
			geom_point(size=2, position=dodge) +
			scale_fill_manual(name="", values=treatment_palette) +
			ylim(c(-20, 25)) +
			xlab("") +
			ylab("Topsoil Temperature, °C") +
			facet_wrap(~string_month, ncol=9) +
			guides(fill=guide_legend(override.aes=list(shape=22, size=20), ncol=2), shape=guide_legend(override.aes=list(size=20))) +
			ggplot_theme(leg_pos="bottom", ang_le=45) +
			theme(axis.text.x=element_blank())

	# add stats to faceted boxplots
	temp_box_five <- bxp_fived + stat_pvalue_manual(stat_test, label="p.adj", hide.ns = TRUE, size=10)

	# export
	export_figs_tabs(paste0(save_box_temp, "five_deep_yearly_temperatures_box_", exp_name), temp_box_five, width=168*2, height=168*1.5)

	# define summer
	summer_season <- c("2013-06-01", "2013-07-01", "2013-08-01")

	# define empty
	fived_differences <- vector()

	# loop through summer
	for (a_month in summer_season) {

		# subset by month and treatment
		# and compute average temp
		summer_fived_warming <- mean(by_month_fived$five_deep[by_month_fived$treatment=="Warming" & by_month_fived$string_month==a_month][!is.na(by_month_fived$five_deep[by_month_fived$treatment=="Warming" & by_month_fived$string_month==a_month])])
		summer_fived_control <- mean(by_month_fived$five_deep[by_month_fived$treatment=="Control" & by_month_fived$string_month==a_month][!is.na(by_month_fived$five_deep[by_month_fived$treatment=="Control" & by_month_fived$string_month==a_month])])
		fived_differences <- rbind.data.frame(fived_differences, cbind.data.frame(a_month, diff=summer_fived_warming-summer_fived_control))

	}

	###################################################################################################

	# finally, create images for publication

	supplementary_fig_four_a <- plot_grid(temp_box_surf, temp_box_five, labels = c('A', 'B'), label_size = 50, ncol=1)

	# export it
	export_figs_tabs(paste0(save_final_figs, "supplementary_fig_four_a"), supplementary_fig_four_a, base=FALSE, width=168*4, height=168*4)

	###################################################################################################

	supplementary_fig_four_b <- plot_grid(timeline_surf, timeline_five, labels = c('A', 'B'), label_size = 50, ncol=1)

	# export it
	export_figs_tabs(paste0(save_final_figs, "supplementary_fig_four_b"), supplementary_fig_four_b, base=FALSE, width=168*4, height=168*4)

}
cat("DON'T WORRY ABOUT THIS WARNING: Removed 48 rows containing non-finite values (`stat_boxplot()`).\nIT IS COMING FROM THE NAs IN THE TEMP DATASET\n")

