# this function was created to sort the issue
# with the merge_samples function from phyloseq
# as mentioned in several threads
# https://github.com/joey711/phyloseq/issues/493
# https://github.com/joey711/phyloseq/issues/243
# https://github.com/joey711/phyloseq/issues/608
# https://github.com/joey711/phyloseq/issues/386

assign_factors_back <- function(map_it, ref_mapping) {

	# create mapping table
	back_to_names <- cbind.data.frame(name=unique(ref_mapping), value=as.numeric(unique(ref_mapping)))

	# create null vector
	str_vector <- NULL
	# loop through the numerics in sample_data
	for (vrb in map_it) {
		# find the position of each numeric
		pos_itions <- which(vrb == back_to_names$value)
		# and append its corresponding character vector
		str_vector <- append(str_vector, as.character(back_to_names$name[pos_itions]))
	}
	# return a characters vector
	return(str_vector)
}

############################################################################### PRINTING BOXPLOTS FOR PAIRED DATA

# this funtion is used to plot environmental variables
# as boxplots using wilcoxon test since we can't assume data has normal distribution
# http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
print_boxplots_paired <- function(meta_data, env_vars, exp_nm, save_path){

	# loop through environmental variables
	for (a_var in env_vars) {

		print(paste0("analysing ", a_var))

		# plot, using t-test, assuming environmental variables have normal distribution
		# see histograms generated above
		# http://www.sthda.com/english/wiki/t-test

		# get data for a single gene
		meta_sub <- meta_data[which(meta_data$variable==a_var), ]

		# check is all the data are paired and discard those that are not
		# create a vector combining season and site
		meta_sub$time_by_site <- paste0(meta_sub$TimePoint, meta_sub$CollectionSite)

		# create palette
		treatment_palette <- hue_pal()(2)[c(2, 1)]
		names(treatment_palette) <- c("Control", "Warming")

		# get the plot ready
		dodge <- position_dodge(width=0.8)
		pl <- ggplot(meta_sub, aes(x=TimePoint, y=value, fill=Treatment)) +
			geom_violin(width=1, position=dodge, linewidth=1) +
			geom_boxplot(width=0.2, color="black", position=dodge, linewidth=1, show.legend=FALSE) +
			geom_point(size=2, alpha=1, position=dodge, show.legend=FALSE) +
			scale_fill_manual(name="Treatment", values=treatment_palette) +
			xlab("") +
			ylab("Standardised unit") +
			ggtitle(a_var) +
			guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=2), shape=guide_legend(override.aes=list(size=20))) +
			ggplot_theme(leg_pos="bottom", ang_le=45)
		# and export
		export_svg(paste0(save_path, "analysing_", a_var, "_", exp_nm), pl)
	}
	return(0)
}

############################################################################### PRINTING BOXPLOTS WITHIN

# this funtion is used to plot environmental variables
# as boxplots
print_boxplots_within <- function(meta_data, env_vars, exp_nm, save_path){

	# define comparisons
	comparisons <- list(c("June", "July"), c("July", "August"))

	# loop through environmental variables
	for (envv in env_vars) {

		# create subset of meta_data using one variable at a time
		meta_sub <- meta_data[which(meta_data$env_var==envv), ]

		# create palette
		treatment_palette <- hue_pal()(2)[c(2, 1)]
		names(treatment_palette) <- c("Control", "Warming")

		# plot, using using wilcoxon test since we can't assume data has normal distribution
		# http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
		dodge <- position_dodge(width=0.8)
		pl <- ggplot(meta_sub, aes(x=TimePoint, y=value, fill=Treatment)) +
			geom_violin(width=1, position=dodge, linewidth=1) +
			geom_boxplot(width=0.2, color="black", position=dodge, linewidth=1, show.legend=FALSE) +
			geom_point(size=2, alpha=1, position=dodge, show.legend=FALSE) +
			scale_fill_manual(name="Treatment", values=treatment_palette) +
			geom_point(position=position_jitterdodge()) +
			ylab("Standardised unit") +
			xlab("") +
			guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=2), shape=guide_legend(override.aes=list(size=20))) +
			ggplot_theme(leg_pos="bottom", ang_le=45)

		# and export
		export_svg(paste0(save_path, "analysing_", envv, "_", exp_nm), pl)

		print(paste0(envv, " done"))
	}
	return(0)
}

############################################################################### RAREFACTION CURVES

#rarefaction curves https://github.com/joey711/phyloseq/issues/143
calculate_rarefaction_curves <- function(psdata, measures, depths) {

	set.seed(131)

	estimate_rarified_richness <- function(psdata, measures, depth) {

		if(max(sample_sums(psdata)) < depth) return()

		psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)

		rarified_psdata <- rarefy_even_depth(psdata, depth, verbose=FALSE)

		alpha_diversity <- estimate_richness(rarified_psdata, measures=measures)

		# as.matrix forces the use of melt.array, which includes the Sample names (rownames)
		molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames=c("SampleID", "Measure"), value.name="Alpha_diversity")

		return(molten_alpha_diversity)
	}

	# enable automatic addition of the Depth to the output by ldply
	names(depths) <- depths
	# with dply, the parameters with the . are specific parameters of the function
	rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata=psdata, measures=measures, .id="Depth", .progress=ifelse(interactive(), "text", "none"))

	# convert Depth from factor to numeric
	rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]

	return(rarefaction_curve_data)
}

############ PLOTTING FOR PARTIAL RDA

# define plotting rda with base r
plot_ordination <- function(data, metadata, var_fitting, scale_param, arrlen) {

	# adding species to triplot as mentioned here
	# https://github.com/vegandevs/vegan/issues/341

	# define shapes for timepoint legend
	shape_legend <- c(16, 18, 17)
	# define shapes for timepoint. both shapes must be in the same order
	shape <- c(21, 23, 24)
	names(shape) <- unique(metadata$TimePoint)

	# define colours for treatment
	treatment_palette <- hue_pal()(2)[c(2, 1)]
	names(treatment_palette) <- c("Control", "Warming")

	# set levels as factor
	metadata$Treatment <- factor(metadata$Treatment, levels=c("Control", "Warming"))
	metadata$TimePoint <- factor(metadata$TimePoint, levels=c("June", "July", "August"))

	# get multiplying factors for the arrows
	arrow_factor <- ordiArrowMul(var_fitting)

	# adding species to triplot as mentioned here
	# https://github.com/vegandevs/vegan/issues/341

	# plot
	par(mar=c(5, 5, 4, 2))
	pl <- plot(data, type="none", cex.lab=2.5, cex.axis=2.5, cex.main=3.5, font=2, scaling=scale_param, main="")
#	splen <- sqrt(rowSums(pl$species^2)) # to add species
#	with(metadata, text(pl, "species", select=splen > arrlen, arrow=TRUE, length=0.05)) # to add species
	with(metadata, plot(var_fitting, add=T, col="black", cex=3, font=2, p.max=0.05, arrow.mul=arrow_factor))
	with(metadata, points(data, "sites", pch=shape[TimePoint], bg=treatment_palette[Treatment], col="black", cex=8, lwd=5, scaling=scale_param))
	with(metadata, legend("topright", legend=levels(Treatment), fill=treatment_palette, title="Treatment", cex=4))
	with(metadata, legend("topleft", legend=levels(TimePoint), pch=shape_legend, title="Time point", cex=4))
#	with(metadata, text(data, cex=0.9, scaling=scale_param)) # to add sample names
}

# define plotting rda with ggplot2
plot_ordination_ggplot <- function(data, var_fitting, all_stats, v_exp, scale_param, pattern=FALSE) {

	# create data frame for ggplot plotting
	# get scores for sites
	ggplot_rda <- as.data.frame(vegan::scores(data, display="sites", scaling=scale_param))
	ggplot_rda$label <- rownames(ggplot_rda)

	# collect info about treatment
	treats <- c(nrow(ggplot_rda))
	treats[grep("TW", ggplot_rda$label)] <- "Warming"
	treats[grep("TC", ggplot_rda$label)] <- "Control"

	# collect info about treatment
	xna <- c(nrow(ggplot_rda))
	xna[grep("AD006", ggplot_rda$label)] <- "DNA"
	xna[grep("AD012", ggplot_rda$label)] <- "RNA"

	# collect info about timepoint
	times <- c(nrow(ggplot_rda))
	times[grep("June", ggplot_rda$label)] <- "June"
	times[grep("July", ggplot_rda$label)] <- "July"
	times[grep("August", ggplot_rda$label)] <- "August"

	# add the info just collect as columns to ggplot_rda
	ggplot_rda$Treatment <- treats
	ggplot_rda$TimePoint <- times
	ggplot_rda$DNA_RNA <- xna

	# set factors
	ggplot_rda$Treatment <- factor(ggplot_rda$Treatment, levels=c("Control", "Warming"))
	ggplot_rda$TimePoint <- factor(ggplot_rda$TimePoint, levels=c("June", "July", "August"))
	ggplot_rda$DNA_RNA <- factor(ggplot_rda$DNA_RNA, levels=c("DNA", "RNA"))

	# get multiplying factors for the arrows
	arrow_factor <- ordiArrowMul(var_fitting)

	# get arrows and scale them
	arrows <- as.data.frame(vegan::scores(var_fitting, display="vectors")) * arrow_factor
	arrows$names <- rownames(arrows)

	# select arrows which have p-value <= 0.05
	signif_arrows <- names(which(var_fitting$vectors$pvals<=0.05))
	# subset
	print_arrows <- arrows[signif_arrows, ]

	# define shapes for timepoint. both shapes must be in the same order
	time_shape <- c(21, 23, 24)
	names(time_shape) <- unique(ggplot_rda$TimePoint)

	# define colours for treatment
	treatment_palette <- hue_pal()(2)[c(2, 1)]
	names(treatment_palette) <- c("Control", "Warming")

	# define palette for xna
	xna_palette <- brewer.pal(n=4, name = "Set2")[c(1,4)]
	names(xna_palette) <- c("DNA", "RNA")

	if (pattern == TRUE) {

		# ggplotting
		rda_with_ggplot <- ggplot(ggplot_rda) +
				geom_point(mapping=aes(x=RDA1, y=RDA2, fill=Treatment, shape=TimePoint, color=DNA_RNA), size=20, stroke=7) +
				coord_fixed() +
				scale_shape_manual(values=time_shape) +
				scale_color_manual(values=xna_palette) +
				scale_fill_manual(values=treatment_palette) +
				geom_segment(data=print_arrows, aes(x=0, xend=RDA1, y=0, yend=RDA2), arrow=arrow(length=unit(0.25, "cm")), colour="black", lwd=1) +
				geom_text_repel(data=print_arrows, aes(x=RDA1, y=RDA2, label=names), col="black", fontface="bold", size=15, max.overlaps=20) +
				guides(color=guide_legend(override.aes=list(shape=22, size=20)), fill=guide_legend(override.aes=list(shape=22, size=20)), shape=guide_legend(override.aes=list(size=20))) +
				labs(x = paste("RDA1 (",round(v_exp[2,"RDA1"]*100,digits = 2),"%)",sep = ""), y = paste("RDA2 (",round(v_exp[2,"RDA2"]*100,digits = 2),"%)",sep = "")) +
				ggplot_theme("bottom") +
				theme(aspect.ratio=1.5)
	} else {

		efs <- all_stats$F
		pvals <- all_stats$pval

		# ggplotting
		rda_with_ggplot <- ggplot(ggplot_rda) +
				geom_point(mapping=aes(x=RDA1, y=RDA2, fill=Treatment, shape=TimePoint), color="black", size=20, stroke=3) +
				coord_fixed() +
				xlim(c(-2, 2)) +
				annotate("text", x=-1.7, y=-1.3, label= "Treatment", size=15, fontface="bold") +
				annotate("text", x=-1.7, y=-1.4, label= "Timepoint", size=15, fontface="bold") +
				annotate("text", x=-0.9, y=-1.2, label="F", size=15, fontface="bold") +
				annotate("text", x=-0.9, y=-1.3, label= round(efs[1], digits=2), size=15, fontface="bold") +
				annotate("text", x=-0.9, y=-1.4, label= round(efs[2], digits=2), size=15, fontface="bold") +
				annotate("text", x=-0.2, y=-1.2, label="p-value", size=15, fontface="bold") +
				annotate("text", x=-0.2, y=-1.3, label= round(pvals[1], digits=3), size=15, fontface="bold") +
				annotate("text", x=-0.2, y=-1.4, label= round(pvals[2], digits=3), size=15, fontface="bold") +
				scale_shape_manual(values=time_shape) +
				scale_fill_manual(values=treatment_palette) +
				geom_segment(data=print_arrows, aes(x=0, xend=RDA1, y=0, yend=RDA2), arrow=arrow(length=unit(0.25, "cm")), colour="black", lwd=1) +
				geom_text_repel(data=print_arrows, aes(x=RDA1, y=RDA2, label=names), col="black", fontface="bold", size=15, max.overlaps=20) +
				guides(fill=guide_legend(override.aes=list(shape=22, size=20)), shape=guide_legend(override.aes=list(size=20))) +
				labs(x = paste("RDA1 (",round(v_exp[2,"RDA1"]*100,digits = 2),"%)",sep = ""), y = paste("RDA2 (",round(v_exp[2,"RDA2"]*100,digits = 2),"%)",sep = "")) +
				ggplot_theme("bottom") +
				theme(aspect.ratio=1.5)
	}

	return(rda_with_ggplot)
}

# create consistent theme for all ggplots
ggplot_theme <- function(leg_pos="right", ang_le=0, fnt_sz=50) {

	if (ang_le == 45) {
		# set theme
		theme_bw() +
		theme(text=element_text(size=fnt_sz, color="black"), legend.position=leg_pos, legend.text=element_text(size=fnt_sz), legend.key.height=unit(1,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=ang_le, hjust=1, vjust=1, color="black", size=fnt_sz), axis.text.y=element_text(color="black", size=fnt_sz), strip.text=element_text(size=fnt_sz), panel.spacing = unit(3, "lines"))
	} else if (ang_le == 90) {
		# set theme
		theme_bw() +
		theme(text=element_text(size=fnt_sz, color="black"), legend.position=leg_pos, legend.text=element_text(size=fnt_sz), legend.key.height=unit(1,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=ang_le, hjust=0.5, vjust=1, color="black", size=fnt_sz), axis.text.y=element_text(color="black", size=fnt_sz), strip.text=element_text(size=fnt_sz), panel.spacing = unit(3, "lines"))
	} else {
		# set theme
		theme_bw() +
		theme(text=element_text(size=fnt_sz, color="black"), legend.position=leg_pos, legend.text=element_text(size=fnt_sz), legend.key.height=unit(1,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(color="black", size=fnt_sz), axis.text.y=element_text(color="black", size=fnt_sz), strip.text=element_text(size=fnt_sz), panel.spacing = unit(3, "lines"))
	}
}

# export images in case they are needed
export_svg <- function(filename, a_plot, width=22, height=22, base=F, as_rds=F) {

	# if plot is made with ggplot2
	if (base==F) {
		# plot
		svg(paste0(filename, ".svg"), width=width, height=height)
			plot(a_plot)
		dev.off()

	# else with base R
	} else {
		# export plot as svg
		svg(paste0(filename, ".svg"), width=width, height=height)
			a_plot
		dev.off()
	}
	if (as_rds ==T) {
		saveRDS(a_plot, paste0(filename, ".Rds"))
	}
}

# plot that mimic the pairs function of base R
ggplot_pairs <- function(data, vars_to_plot, colour_by, n_cols, normalised=T) {

	# normalised?
	if (normalised ==T ) {
		vars_to_plot <- paste0("norm_", vars_to_plot)
	}

	# melt the table
	melted <- melt(data[, c(colour_by, vars_to_plot)])
	# remove the norm from variable
	melted$variable <- gsub("norm_", "", melted$variable)

	# define empty vector
	final_all <- vector()

	# for each interesting variable
	for (a_var in vars_to_plot) {

		# get it's actual values
		temp <- data[, a_var]

		# replicate them for each variable
		temp_col <- rep(temp, length(unique(melted$variable)))

		# rbind them
		final_all <- rbind.data.frame(final_all, cbind(melted, var=rep(gsub("norm_", "", a_var), length(temp_col)), temp_col))
	}

	# remove duplicate comparisons
	final_no_dup <- final_all[-which(final_all$variable==final_all$var), ]

	# assign levels for plotting
	final_no_dup$variable <- factor(final_no_dup$variable, levels=unique(final_no_dup$variable))
	final_no_dup$var <- factor(final_no_dup$var, levels=unique(final_no_dup$var))

	if (colour_by == "Treatment") {
		# create palette for each datapoint
		treatment_palette <- hue_pal()(2)[c(2, 1)]
		names(treatment_palette) <- c("Control", "Warming")

		# plot
		my_ggpairs <- ggplot(final_no_dup, aes(x = value, y = temp_col, fill=get(colour_by))) +
			geom_point(shape=21, size=12, color="black", stroke=2) +
			scale_fill_manual(name="Treatment", values=treatment_palette) +
			facet_wrap(var ~ variable, ncol=n_cols, labeller = label_wrap_gen(multi_line=FALSE)) +
			xlab("") +
			ylab("") +
			guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=3), shape=guide_legend(override.aes=list(size=20))) +
			theme_bw() +
			theme(text=element_text(size=50, colour="black"), legend.position="bottom", legend.text=element_text(size=80), legend.key.height=unit(0.7,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, colour="black", size=50), axis.text.y=element_text(colour="black", size=50), strip.text=element_text(size=60)) +
			theme(legend.title=element_blank(), panel.spacing = unit(3, "lines"))

	} else {
		# create palette for each datapoint
		timepoint_palette <- magma(3)
		names(timepoint_palette) <- c("June", "July", "August")
		# plot
		my_ggpairs <- ggplot(final_no_dup, aes(x = value, y = temp_col, fill=get(colour_by))) +
			geom_point(shape=21, size=12, color="black", stroke=2) +
			scale_fill_manual(name="TimePoint", values=timepoint_palette) +
			facet_wrap(var ~ variable, ncol=n_cols, labeller = label_wrap_gen(multi_line=FALSE)) +
			xlab("") +
			ylab("") +
			guides(fill=guide_legend(override.aes=list(shape=21, size=20), ncol=3), shape=guide_legend(override.aes=list(size=20))) +
			theme_bw() +
			theme(text=element_text(size=50, colour="black"), legend.position="bottom", legend.text=element_text(size=80), legend.key.height=unit(0.7,'cm'), plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, colour="black", size=50), axis.text.y=element_text(colour="black", size=50), strip.text=element_text(size=60)) +
			theme(legend.title=element_blank(), panel.spacing = unit(3, "lines"))
		}
	return(my_ggpairs)
}

# imputation for partial RDA
impute_missing <- function(imputed_metadata, phylo_obj, sel_vars){

	# for the remaining columns (if any!), fill the NA with
	# the median value computed using the available values
	# that are collected from the biological replicates
	for (clm_index in seq(1, length(sel_vars), by=1)) {

		# get the name of the column of interest
		clm <- sel_vars[clm_index]
		# check if there is a NA in the column
		any_true <- any(is.na(imputed_metadata[, clm])==TRUE)
		# if yes, compute average/median and substitute
		if (any_true == TRUE) {
			print(paste0("imputing ", clm))
			# find NAs
			is_na <- which(is.na(imputed_metadata[, clm]))
			# for each NA that is found, compute the average using
			# samples from same timepoint and treatment
			for (sn in is_na) {
				# find sample treatment
				treat <- imputed_metadata$Treatment[sn]
				# find sample timepoint
				timep <- imputed_metadata$TimePoint[sn]
				# find sample xna
				xna <- imputed_metadata$DNA_RNA[sn]
				# get all the datapoints with that treatment and that timepoint
				interesting_data <- which(imputed_metadata$TimePoint==timep & imputed_metadata$Treatment==treat & imputed_metadata$DNA_RNA==xna)
				all_data <- imputed_metadata[interesting_data, clm]
				# compute median with the available values
				average <- median(unlist(all_data[!is.na(all_data)]))
				# substitute the NAs into the original column
				imputed_metadata[sn, clm] <- average
			}
		}
	}

	# make sure the variables that will be used in the formula are
	# factors since statistical models require it, to work properly
	imputed_metadata$Treatment <- ifelse(imputed_metadata$Treatment=="C", "Control", "Warming")
	imputed_metadata$Treatment <- factor(imputed_metadata$Treatment, levels=c("Control", "Warming"))
	imputed_metadata$TimePoint <- factor(imputed_metadata$TimePoint, levels=c("June", "July", "August"))
	imputed_metadata$CollectionSite <- factor(imputed_metadata$CollectionSite, levels=unique(imputed_metadata$CollectionSite))
	imputed_metadata$DNA_RNA <- factor(imputed_metadata$DNA_RNA, levels=c("DNA", "RNA"))

	# build imputed phylo object
	imputed_phylo <- phyloseq(otu_table(phylo_obj), imputed_metadata, tax_table(phylo_obj))

	return(imputed_phylo)
}

# imputation for partial RDA
impute_missing_xna <- function(imputed_metadata, phylo_obj, sel_vars){

	# for the remaining columns (if any!), fill the NA with
	# the median value computed using the available values
	# that are collected from the biological replicates
	for (clm_index in seq(1, length(sel_vars), by=1)) {

		# get the name of the column of interest
		clm <- sel_vars[clm_index]
		# check if there is a NA in the column
		any_true <- any(is.na(imputed_metadata[, clm])==TRUE)
		# if yes, compute average/median and substitute
		if (any_true == TRUE) {
			print(paste0("imputing ", clm))
			# find NAs
			is_na <- which(is.na(imputed_metadata[, clm]))
			# for each NA that is found, compute the average using
			# samples from same timepoint and treatment
			for (sn in is_na) {
				# find sample treatment
				treat <- imputed_metadata$Treatment[sn]
				# find sample timepoint
				timep <- imputed_metadata$TimePoint[sn]
				# get all the datapoints with that treatment and that timepoint
				interesting_data <- which(imputed_metadata$TimePoint==timep & imputed_metadata$Treatment==treat)
				all_data <- imputed_metadata[interesting_data, clm]
				# compute median with the available values
				average <- median(unlist(all_data[!is.na(all_data)]))
				# substitute the NAs into the original column
				imputed_metadata[sn, clm] <- average
			}
		}
	}

	# make sure the variables that will be used in the formula are
	# factors since statistical models require it, to work properly
	imputed_metadata$Treatment <- ifelse(imputed_metadata$Treatment=="C", "Control", "Warming")
	imputed_metadata$Treatment <- factor(imputed_metadata$Treatment, levels=c("Control", "Warming"))
	imputed_metadata$TimePoint <- factor(imputed_metadata$TimePoint, levels=c("June", "July", "August"))
	imputed_metadata$CollectionSite <- factor(imputed_metadata$CollectionSite, levels=unique(imputed_metadata$CollectionSite))

	# build imputed phylo object
	imputed_phylo <- phyloseq(otu_table(phylo_obj), imputed_metadata, tax_table(phylo_obj))

	return(imputed_phylo)
}

first_capital <- function(strng) {
	substr(strng, 1, 1) <- toupper(substr(strng, 1, 1))
	return(strng)
}
