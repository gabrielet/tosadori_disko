# source("010_ordinations_partial_rda_clr.R")

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
library("ggvegan")
library("RColorBrewer")

# import functions
source("/microbiology/disko2013/code/000_micro_functions_disko2013.R")

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
save_orderings_partial <- paste0(save_orderings, "partial/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_phylo), dir.create(save_phylo, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_orderings_partial), "dir exists!", dir.create(save_orderings_partial, recursive=TRUE))

# load centered-transformed object
centered_phylo_data <- readRDS(paste0(save_phylo, "centered_phylo_data_", exp_name, ".Rds")) ; guilded <- ""

# set taxa level of interest
t_level <- "ASV"

# get metadata
centered_metadata <- sample_data(centered_phylo_data)
# rename cols for plotting
colnames(centered_metadata) <- c("sampleID", "DNA_RNA", "Treatment", "TimePoint", "CollectionSite", "SM", "SM_perc", "pH", "NH4", "NH4_SOM", "NO3", "NO3_SOM", "PO4", "PO4_SOM", "N", "C", "mcrA", "mxaF", "nifH", "nirS", "nosZ", "TOC", "TOC_SOM", "microC", "microC_SOM", "DON", "DON_SOM", "TN", "TN_SOM", "microN", "microN_SOM", "microCN", "d15N", "d13C", "SOM", "LOI", "Deciduous_shrub", "Evergreen_shrub", "Lichen", "Forb", "Litter", "plant_diversity", "plant_eveness", "soil_fresh", "soil_dry", "soil_water", "CN")

# load centered-transformed qpcr object
centered_phylo_qpcr <- readRDS(paste0(save_phylo, "centered_phylo_qpcr_", exp_name, ".Rds"))

# get qpcr metadata
centered_metadata_qpcr <- sample_data(centered_phylo_qpcr)
# rename cols for plotting
colnames(centered_metadata_qpcr) <- c("sampleID", "DNA_RNA", "Treatment", "TimePoint", "CollectionSite", "mcrA", "mxaF", "nifH", "nirS", "nosZ", "soil_fresh", "soil_dry")

## and the features selections
#forward_selection <- list()

######################################################################## PERFORMING PARTIAL RDA ON DNA and RNA SEPARATELY

print("PERFORMING PARTIAL RDA ON DNA AND RNA, SEPARATELY")

# list environmental pararameters that might require imputation
selected_vars_qpcr <- colnames(centered_metadata)[c(6:ncol(centered_metadata))]

# define list for agglomerated, using tax_glom(), phylo objects
phylo_glom_xna <- list()

# define lists to store the results from the ordinations
r_square_partial <- list()
signif_test_all_vars_partial <- list()
mod_selected <- list()
var_fitting_enviro <- list()
imputed_meta_glom <- list()
stats_to_plot <- list()
variance_expl <- list()

# define formula for computing the ordination
formula_rda <- "phylo_glom_xna[[xna]]~Treatment+TimePoint+Condition(CollectionSite)" ; withint <- ""

# define variables for the fitting
if (orgn == "bacteria") {
	selected_vars <- c("pH", "SM", "C", "N", "CN", "SOM", "LOI", "TOC", "TN", "DON", "NH4", "NO3", "PO4", "d15N", "d13C", "microC", "microN", "microCN")
} else if (orgn == "fungi") {
	selected_vars <- c("pH", "SM", "C", "N", "CN", "SOM", "LOI", "TOC", "TN", "DON", "NH4", "NO3", "PO4", "d15N", "d13C", "microC", "microN", "microCN", "Deciduous_shrub", "Evergreen_shrub", "Lichen", "Forb", "Litter", "plant_diversity", "plant_eveness")
}

# loop through the xnas
for (xna in c("DNA", "RNA")) {

	print(paste0("analysing ", xna))

	# subset cols
	centered_metadata_subset <- centered_metadata[, c("sampleID", "DNA_RNA", "Treatment", "TimePoint", "CollectionSite", selected_vars)]

	# impute metadata for xna
	imputed_phylo <- impute_missing_xna(centered_metadata_subset[which(centered_metadata_subset$DNA_RNA==xna)], centered_phylo_data, selected_vars)

	# get imputed metadata
	imputed_metadata <- sample_data(imputed_phylo)

	# get DNA or RNA only, from the imputed_phylo object4
	phylo_xna <- prune_samples(as.character(imputed_metadata$sampleID[which(imputed_metadata$DNA_RNA == xna)]), imputed_phylo)

	# remove ASVs which are always zero for the subset under consideration
	# to do so, with clr transformation, find how many zeros are found for each row. if the
	# number of zeros per each row is equal to the length of the row, remove the taxa. here,
	# the rows which length is different are selected because of the PRUNE way of thinking
	phylo_xna <- prune_taxa(apply(tax_table(phylo_xna), 1, function(x) length(which(x==0))!=length(x)), phylo_xna)

	# at this point it is required to decide whether the analysis will be performed
	# at ASV level or at a specific taxonomic level. If ASV level, go ahead. If ASV
	# level, then the following needs to be un-commented and the following line needs
	# to be commented in order to not let it work

	# if not aggregating at a specific taxonomic level
	if(t_level=="ASV" | t_level=="OTU") {

		# create object with same name as if it were agglomerated
		# to keep it consistent with the else option, below
		phylo_glom_xna[[xna]] <- phylo_xna

		# if not guilded use names at a specific taxonomic level
		if (guilded == "") {
			# set new names for taxa, using the level of interest
			taxa_names(phylo_glom_xna[[xna]]) <- make.unique(tax_table(phylo_glom_xna[[xna]])[, "Phylum"])
		} else {
			# otherwise set names using guild info
			taxa_names(phylo_glom_xna[[xna]]) <- make.unique(tax_table(phylo_glom_xna[[xna]])[, 2])
		}

		# set max distance for taxa to show up in the plot
		arr_length <- 0.25
	} else {
		# merge by taxa level of interest
		phylo_glom_xna[[xna]] <- tax_glom(phylo_xna, taxrank=t_level)

		# if not guilded use names at a specific taxonomic level
		if (guilded == "") {
			# set new names for taxa, using the level of interest
			taxa_names(phylo_glom_xna[[xna]]) <- make.unique(tax_table(phylo_glom_xna[[xna]])[, t_level])
		} else {
			# otherwise set names using guild info
			taxa_names(phylo_glom_xna[[xna]]) <- make.unique(tax_table(phylo_glom_xna[[xna]])[, 2])
		}

		# set max distance for taxa to show up in the plot
		arr_length <- 0.5
	}

	# get metadata for aggregated-imputed object
	imputed_meta_glom[[xna]] <- as(sample_data(phylo_glom_xna[[xna]]), "data.frame")

	# set levels to make Controls the reference level
	imputed_meta_glom[[xna]]$Treatment <- factor(imputed_meta_glom[[xna]]$Treatment, levels=unique(imputed_meta_glom[[xna]]$Treatment))
	imputed_meta_glom[[xna]]$Treatment <- relevel(imputed_meta_glom[[xna]]$Treatment, ref="Control")
	imputed_meta_glom[[xna]]$TimePoint <- factor(imputed_meta_glom[[xna]]$TimePoint, levels=unique(imputed_meta_glom[[xna]]$TimePoint))
	imputed_meta_glom[[xna]]$TimePoint <- relevel(imputed_meta_glom[[xna]]$TimePoint, ref="June")

	# loop through scaling parameters
	for (scale_parameter in c(1, 2)) {

		set.seed(131)

		# compute partial RDA using formula
		invisible(capture.output(mod_selected[[xna]] <- ordinate(phylo_glom_xna[[xna]], "RDA", formula=as.formula(formula_rda))))
		# same as doing
#		otu_mat <- as(otu_table(phylo_glom_xna[[xna]]), "matrix")
#		rda_data <- data.frame(sample_data(phylo_glom_xna[[xna]]))
#		rda(t(otu_mat)~Treatment+TimePoint+Condition(CollectionSite), data=rda_data)
		# https://www.rdocumentation.org/packages/vegan/versions/2.6-10/topics/varpart
#		r_squared <- varpart(t(otu_mat), ~Treatment, ~TimePoint, partial=~CollectionSite, data=rda_data)
#		r_squared$part$indfrac

		# variance explained by axes
		variance_expl[[xna]] <- as.data.frame(summary(mod_selected[[xna]])$cont$importance)

		# compute fitting
		var_fitting_enviro[[xna]] <- envfit(mod_selected[[xna]], imputed_meta_glom[[xna]][, selected_vars], perm=1000, na.rm=TRUE, scaling=scale_parameter)

		# get adjusted R2
		r_square_partial[[xna]] <- RsquareAdj(mod_selected[[xna]])$adj.r.squared

		# get p-value
		signif_test_all_vars_partial[[xna]] <- anova.cca(mod_selected[[xna]], step=1000, by="terms")

		# variance explained to R2
		print(signif_test_all_vars_partial[[xna]]$Variance / sum(signif_test_all_vars_partial[[xna]]$Variance))

		# get stats to plot
		stats_to_plot[[xna]] <- signif_test_all_vars_partial[[xna]][c(1,2), c(3,4)]
		colnames(stats_to_plot[[xna]]) <- c("F", "pval")
		# create ggplot plot
		ggplot_scaling <- plot_ordination_ggplot(data=mod_selected[[xna]], var_fitting=var_fitting_enviro[[xna]], all_stats=stats_to_plot[[xna]], v_exp=variance_expl[[xna]], scale_param=scale_parameter)
		# export ggplot
		export_svg(paste0(save_orderings_partial, "ggplot2_partial_RDA_", xna, "_", exp_name, "_", t_level, "_scaling_", scale_parameter, "_params_clr", withint, guilded), ggplot_scaling, base=FALSE, as_rds=T)
	}
}

######################################################################## PERFORMING PARTIAL RDA ON QPCR GENES

if (orgn == "bacteria") {
	print("PERFORMING PARTIAL RDA on qPCR")

	# list qpcr pararameters that might require imputation
	selected_vars_qpcr <- c("nirS", "nosZ", "nifH", "mcrA", "mxaF")

	# make a copy of centered_metadata qpcr to perform imputation
	imputed_phylo_qpcr <- centered_phylo_qpcr#impute_missing(centered_metadata_qpcr, centered_phylo_qpcr, selected_vars_qpcr)

	# build imputed phylo qpcr object
	imputed_metadata_qpcr <- centered_metadata_qpcr#sample_data(imputed_phylo_qpcr)

	# at this point it is required to decide whether the analysis will be performed
	# at ASV level or at a specific taxonomic level. If ASV level, go ahead. If ASV
	# level, then the following needs to be un-commented and the following line needs
	# to be commented in order to not let it work

	# if not aggregating at a specific taxonomic level
	if(t_level=="ASV" | t_level=="OTU") {
		phylo_glom_qpcr <- imputed_phylo_qpcr

		# set new names for taxa, using the level of interest
		taxa_names(phylo_glom_qpcr) <- make.unique(tax_table(phylo_glom_qpcr)[, "Phylum"])

		# set max distance for taxa to show up in the plot
		arr_length <- 0.25
	} else {
		# merge by taxa level of interest
		phylo_glom_qpcr <- tax_glom(imputed_phylo_qpcr, taxrank=t_level)

		# set new names for taxa, using the level of interest
		taxa_names(phylo_glom_qpcr) <- make.unique(tax_table(phylo_glom_qpcr)[, t_level])

		# set max distance for taxa to show up in the plot
		arr_length <- 0.5
	}

	# get metadata for aggregated-imputed objects
	imputed_enviro_qpcr <- as(merge(as(imputed_meta_glom[["DNA"]][, c("sampleID", "Treatment", "TimePoint", "CollectionSite", selected_vars)], "data.frame"), as(imputed_metadata_qpcr, "data.frame")[, c("sampleID", selected_vars_qpcr)], by="sampleID"), "data.frame")

	# make sure the variables in the formula are factors since
	# statistical models require it, to work properly
	imputed_enviro_qpcr$Treatment <- ifelse(imputed_enviro_qpcr$Treatment=="C", "Control", "Warming")
	imputed_enviro_qpcr$Treatment <- factor(imputed_enviro_qpcr$Treatment, levels=c("Control", "Warming"))
	imputed_enviro_qpcr$TimePoint <- factor(imputed_enviro_qpcr$TimePoint, levels=c("June", "July", "August"))
	imputed_enviro_qpcr$CollectionSite <- factor(imputed_enviro_qpcr$CollectionSite, levels=unique(imputed_enviro_qpcr$CollectionSite))

	# set levels to make Controls the reference level
	imputed_enviro_qpcr$Treatment <- relevel(imputed_enviro_qpcr$Treatment, ref="Control")
	imputed_enviro_qpcr$TimePoint <- relevel(imputed_enviro_qpcr$TimePoint, ref="June")

	# NOTE: the result of this ordination should be the same, if the formula
	# is the same, as the one that was obtained for DNA since, we are using
	# the same abundances. What will change is the fitting which will be done
	# using qPCR variables instead of environmental variables

	# loop through scaling parameters
	for (scale_parameter in c(1)) { #, 2)) {

		set.seed(131)

		# compute fitting on qPCR variables as Jana asked for the paper
		# here, use only DNA since qPCR info are available for DNA only
		# this step is adding qPCR variables to the already fit environmental
		# varibales, to obtain a single plot with all the variables
		var_fitting_qpcr <- envfit(mod_selected[["DNA"]], imputed_enviro_qpcr[, c(selected_vars, selected_vars_qpcr)], perm=1000, na.rm=TRUE, scaling=scale_parameter)

		# get stats to plot
		stats_to_plot_qpcr <- signif_test_all_vars_partial[["DNA"]][c(1,2), c(3,4)]
		colnames(stats_to_plot_qpcr) <- c("F", "pval")

		# create ggplot plot
		ggplot_scaling_qpcr <- plot_ordination_ggplot(data=mod_selected[["DNA"]], var_fitting=var_fitting_qpcr, all_stats=stats_to_plot_qpcr, v_exp=variance_expl[["DNA"]], scale_param=scale_parameter)

		# export ggplot
		export_svg(paste0(save_orderings_partial, "ggplot2_partial_RDA_on_qPCR_genes_and_enviro_", exp_name, "_", t_level, "_scaling_", scale_parameter, "_clr", withint, guilded), ggplot_scaling_qpcr, base=FALSE, as_rds=T)

	}

	# export info in text file
	sink(paste0(save_orderings_partial, "var_fitting_pRDA_qpcr_", t_level, ".txt"))
		print(var_fitting_qpcr)
	sink()
}
######################## FINALLY, EXPORT EVERYTHING

## export forward feature selection
#saveRDS(forward_selection, paste0(save_orderings_partial, "centered_forward_selection_", exp_name, ".Rds"))

# export info in text files
sink(paste0(save_orderings_partial, "signif_pRDA_", t_level, ".txt"))
	print(signif_test_all_vars_partial)
sink()

sink(paste0(save_orderings_partial, "r_squares_pRDA_", t_level, ".txt"))
	print(r_square_partial)
sink()

sink(paste0(save_orderings_partial, "var_fitting_pRDA_", t_level, ".txt"))
	print(var_fitting_enviro)
sink()

# export adjusted R2
saveRDS(r_square_partial, paste0(save_orderings_partial, "r_square_partial_", exp_name, "_", t_level, "_clr", withint, guilded, ".Rds"))

# export signif_test_all_vars_partial
saveRDS(signif_test_all_vars_partial, paste0(save_orderings_partial, "signif_test_all_vars_partial_", exp_name, "_", t_level, "_clr", withint, guilded, ".Rds"))

######################## DNA VERSUS RNA

# set taxa level of interest
t_level <- "ASV"

formula_rda <- "phylo_glom~DNA_RNA+TimePoint+Treatment+Condition(CollectionSite)" ; withint <- ""

# impute metadata for xna
imputed_phylo <- impute_missing(centered_metadata, centered_phylo_data, selected_vars)

# get imputed metadata
imputed_metadata <- sample_data(imputed_phylo)

# set levels to make Controls the reference level
imputed_metadata$Treatment <- ifelse(imputed_metadata$Treatment=="C", "Control", "Warming")
imputed_metadata$Treatment <- factor(imputed_metadata$Treatment, levels=c("Control", "Warming"))
imputed_metadata$TimePoint <- factor(imputed_metadata$TimePoint, levels=c("June", "July", "August"))
imputed_metadata$DNA_RNA <- factor(imputed_metadata$DNA_RNA, levels=c("DNA", "RNA"))
imputed_metadata$Treatment <- relevel(imputed_metadata$Treatment, ref="Control")
imputed_metadata$TimePoint <- relevel(imputed_metadata$TimePoint, ref="June")
imputed_metadata$DNA_RNA <- relevel(imputed_metadata$DNA_RNA, ref="DNA")

# at this point it is required to decide whether the analysis will be performed
# at ASV level or at a specific taxonomic level. If ASV level, go ahead. If ASV
# level, then the following needs to be un-commented and the following line needs
# to be commented in order to not let it work

# if not aggregating at a specific taxonomic level
if(t_level=="ASV" | t_level=="OTU") {

	# create object with same name as if it were agglomerated
	# to keep it consistent with the else option, below
	phylo_glom <- imputed_phylo

	# if not guilded use names at a specific taxonomic level
	if (guilded == "") {
		# set new names for taxa, using the level of interest
		taxa_names(phylo_glom) <- make.unique(tax_table(phylo_glom)[, "Phylum"])
	} else {
		# otherwise set names using guild info
		taxa_names(phylo_glom) <- make.unique(tax_table(phylo_glom)[, 2])
	}

	# set max distance for taxa to show up in the plot
	arr_length <- 0.25
} else {
	# merge by taxa level of interest
	phylo_glom <- tax_glom(imputed_phylo, taxrank=t_level)

	# if not guilded use names at a specific taxonomic level
	if (guilded == "") {
		# set new names for taxa, using the level of interest
		taxa_names(phylo_glom) <- make.unique(tax_table(phylo_glom)[, t_level])
	} else {
		# otherwise set names using guild info
		taxa_names(phylo_glom) <- make.unique(tax_table(phylo_glom)[, 2])
	}

	# set max distance for taxa to show up in the plot
	arr_length <- 0.5
}

# get metadata for aggregated-imputed object
imputed_meta_glom_xna <- sample_data(phylo_glom)

# compute partial RDA using phyloseq function ordinate()
# https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
# the problem concerns the how-to select the most significant variables in the model

# loop through scaling parameters
for (scale_parameter in c(1)) {

	set.seed(131)

	# compute partial RDA using formula
	invisible(capture.output(mod_selected_xna <- ordinate(phylo_glom, "RDA", formula=as.formula(formula_rda))))

	# variance explained by axes
	variance_expl_xna <- as.data.frame(summary(mod_selected_xna)$cont$importance)

	# compute fitting
	var_fitting_enviro_xna <- envfit(mod_selected_xna, imputed_meta_glom_xna[, selected_vars], perm=1000, na.rm=TRUE, scaling=scale_parameter)

	# get adjusted R2
	r_square_partial_xna <- RsquareAdj(mod_selected_xna)$adj.r.squared

	# get p-value
	signif_test_all_vars_partial_xna <- anova.cca(mod_selected_xna, step=1000, by="terms")

	# get stats to plot
	stats_to_plot_xna <- signif_test_all_vars_partial_xna[c(1,2), c(3,4)]
	colnames(stats_to_plot_xna) <- c("F", "pval")

	# create ggplot plot
	ggplot_scaling <- plot_ordination_ggplot(data=mod_selected_xna, var_fitting=var_fitting_enviro_xna, all_stats=stats_to_plot_xna, v_exp=variance_expl_xna, scale_param=scale_parameter, pattern=T)

	# export ggplot
	export_svg(paste0(save_orderings_partial, "ggplot2_partial_RDA_DNA_RNA_", exp_name, "_", t_level, "_scaling_", scale_parameter, "_params_clr", withint, guilded), ggplot_scaling, base=FALSE, as_rds=T)

	# there is also the phyloseq version of the same plots here:
	# https://github.com/joey711/phyloseq/issues/274#issuecomment-30553161
}

