# source("005_dada_analysis.R")

# The data we are using are already demultiplexed and barcodes were already removed. So, the FASTQ files are ready to perform the DADA2 pipeline

# load libraries
library("dada2")
library("ggplot2")
library("phyloseq")
library("DECIPHER")

# import functions
source("/mnt/cinqueg/gabriele/work/microbiology/disko2013/code/000_micro_functions_disko2013.R")

##################### select the experiment #####################
exp_name <- "submission_final_norejects_dada_new_bootstrap" ; orgn <- "bacteria"

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, " EXPERIMENT NAME ", exp_name))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/mnt/cinqueg/gabriele/work/microbiology/disko2013/"
# this is the general path to the experiment
path_to_exp <- paste0(root_path, orgn_dir, "experiments/", exp_name, "/")
# these are the paths where the files from the dada2 analysis go
path_to_dada_files <- paste0(path_to_exp, "dada_files/")
save_img <- paste0(path_to_dada_files, "figures/")
save_data <- paste0(path_to_dada_files, "tables/")

# this is where the cutadapted reads are
cuta_files <- paste0(root_path, orgn_dir, "experiments/cutadapted/dada/")

# if directories do not exist, create it
ifelse(dir.exists(save_img), TRUE, dir.create(save_img, recursive=T))
ifelse(dir.exists(save_data), TRUE, dir.create(save_data, recursive=T))

# get forward and reverse reads
fn_fwds <- sort(list.files(cuta_files, pattern="_R1_", full.names = TRUE))
fn_revs <- sort(list.files(cuta_files, pattern="_R2_", full.names = TRUE))

# extract sample names
sample_names <- vector()
# for each line in fn_fwds
for (ln in c(1:length(fn_fwds))) {
	# split the string and get the first occurrence. also remove the reference to the read direction, i.e. R1 in the case of fn_fwds
	sample_names <- append(sample_names, strsplit(strsplit(basename(fn_fwds[ln]), "\\.")[[1]][1], "_R1")[[1]][1])
}

############ now use a DADA2 function to evaluate the quality of the reads
############ forward and reverse, first two samples
########### plotQualityProfile(fn_fwds[1:2])
########### plotQualityProfile(fn_revs[1:2])

############ Export a svg file for each sample, to get a complete view of the data.
############ forward reads all
###########for (smpl in c(1:length(fn_fwds))){
###########   p <- plotQualityProfile(fn_fwds[smpl])
###########   export_svg(paste0(save_img, "seqs_quality_profile_sample_", smpl, "_forward.svg"), p)
###########}

############ reverse reads all
###########for (smpl in c(1:length(fn_fwds))){
###########   p <- plotQualityProfile(fn_fwds[smpl])
###########   export_svg(paste0(save_img, "seqs_quality_profile_sample_", smpl, "_reverse.svg"), p)
###########}

# NOTE this is an optional step, just to organise the figures in two separated directories.
# To do so, enter the save_img directory from the terminal:
# cd /home/gabriele/microbiology/disko2013/metadata/analysis/cutadapted/figures/
# and copy-paste these commands:
# mkdir forward; mkdir reverse; mv *forward.svg forward; mv *reverse.svg reverse/; exit;

# At this point, having observed the quality profiles of all the samples, we need to define a threshold, if needed, to trim the sequences in order to keep only high-quality, reliable nucleotides. To do this, we need to remember how [Phred scores](https://gatk.broadinstitute.org/hc/en-us/articles/360035531872-Phred-scaled-quality-scores) work. In general, scores > 30 should be acceptable, for amplicon data. The higher the better.

# Create a directory for filtered FASTQ files.

# define a new directory for trimmed and filtered cutadapted files
filtered_path <- paste0(cuta_files, "filtered/")
# create the directory
ifelse(dir.exists(filtered_path), TRUE, dir.create(filtered_path, recursive=T))

# Then prepare the files and perform the trimming step.

# get the files ready for the trimming and filtering
filt_fwds <- paste0(filtered_path, paste0(sample_names, "_F_filt.fastq"))
filt_revs <- paste0(filtered_path, paste0(sample_names, "_R_filt.fastq"))

# trim at the position where the quality of the sequences drops
# IN THIS SPECIFIC CASE, WE HAVE
# - R1 trimmed at 220
# - R2 trimmed ad 220, too since there is no visible drop in quality
# - no need to remove phiX sequences

if (orgn == "bacteria") {
	# 16S calls for truncation
	out <- filterAndTrim(fwd=fn_fwds, filt=filt_fwds, rev=fn_revs, filt.rev=filt_revs, truncLen=c(160, 220), maxN=0, maxEE=c(4, 2), truncQ=2, rm.phix=FALSE, compress=TRUE, multithread=TRUE)
} else if (orgn == "fungi") {
	# no truncation because ITS has variable length
	out <- filterAndTrim(fwd=fn_fwds, filt=filt_fwds, rev=fn_revs, filt.rev=filt_revs, truncQ=2, minLen=50, maxN=0, maxEE=c(4, 4), rm.phix=FALSE, compress=TRUE, multithread=TRUE)
}

# Now, learn the error rates.

# learning the error rates, forward
err_fwd <- learnErrors(filt_fwds, multithread=TRUE)
# and reverse
err_rev <- learnErrors(filt_revs, multithread=TRUE)

# and plot the errors, to verify if everything is ok. forward reads
plot_fwd <- plotErrors(err_fwd, nominalQ=TRUE)
# export plot
export_svg(paste0(save_img, "plot_errors_forward.svg"), plot_fwd)

# same for reverse
plot_rev <- plotErrors(err_rev, nominalQ=TRUE)
# export plot
export_svg(paste0(save_img, "plot_errors_reverse.svg"), plot_rev)

# if a warning is thrown, such as:
# Warning: Transformation introduced infinite values in continuous y-axis
# no need to worry, as reported here https://github.com/benjjneb/dada2/issues/742

# Dereplicate, merge paired ends an perform sample inference.

# performing dereplication, forward
#################################derep_fwds <- derepFastq(filt_fwds, verbose=FALSE)
################################## and reverse
#################################derep_revs <- derepFastq(filt_revs, verbose=FALSE)

################################## name the derep-class objects by the sample names
#################################names(derep_fwds) <- sample_names
#################################names(derep_revs) <- sample_names

# perform sample inference, forward
dada_fwds <- dada(filt_fwds, err=err_fwd, multithread=TRUE, pool="pseudo")
# and reverse
dada_revs <- dada(filt_revs, err=err_rev, multithread=TRUE, pool="pseudo")

# and merge paired end reads, if there are R1 and R2
# running with parameters that should be comparable with the ones from SEED2
# i.e. min overlap = 40 and max perc of differences = 15
# for bacteria
mergers <- mergePairs(dada_fwds, filt_fwds, dada_revs, filt_revs, verbose=FALSE, returnRejects=FALSE, minOverlap=12, maxMismatch=0)
# Inspect the merger data.frame from the first sample
# head(mergers[1])

# create the sequence table
seq_tab <- makeSequenceTable(mergers)

# export seq_tab, before removing chimeras
saveRDS(seq_tab, paste0(save_data, "seq_tab.Rds"))

# The sequences being tabled vary in length.
# dim(seq_tab)

# Inspect distribution of sequence lengths
# table(nchar(getSequences(seq_tab)))

# remove chimeras
seq_tab_no_chim <- removeBimeraDenovo(seq_tab, method="consensus", multithread=TRUE, verbose=FALSE)

# get the table that summarises how many reads were discarded for each filtering step
get_names <- function(x) sum(getUniques(x))
filter_info <- cbind.data.frame(out, sapply(dada_fwds, get_names), sapply(dada_revs, get_names), sapply(mergers, get_names), rowSums(seq_tab_no_chim))
colnames(filter_info) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(filter_info) <- sample_names

# n. of input sequences
print(paste0("n. of input sequences ", filter_info$input))

# percentage of filtered
print(paste0("percentage of filtered ", round(filter_info$merged/filter_info$input, digits=2)))

# percentage of chimeras
print(paste0("percentage of non chimeras ", round(sum(seq_tab_no_chim)/sum(seq_tab), digits=2)))

# export filtering information
write.table(filter_info, paste0(save_data, "dada_filtering_steps.csv"), col.names = T, row.names = T, sep = "\t")

# export the final dada2 table clean from chimeras
saveRDS(seq_tab_no_chim, paste0(save_data, "counts_table_ASV.Rds"))

#set seed for reproducibility
set.seed(131)

# select min bootstrap to assign a taxa
boot <- 80

# run assignment for ASVs
taxa_asv <- assignTaxonomy(seq_tab_no_chim, "/mnt/cinqueg/gabriele/work/microbiology/databases/silva_nr99_v138.1_train_set.fa.gz", minBoot=boot, outputBootstraps=TRUE, tryRC=TRUE, multithread=TRUE)

# select the clustering method
clust_method <- "ASV"

# this is the path where counts tables and taxonomy go
path_to_taxa <- paste0(path_to_exp, clust_method, "/taxonomy/")
# create the folder
ifelse(dir.exists(path_to_taxa), TRUE, dir.create(path_to_taxa, recursive=T))

# export taxonomy table
saveRDS(taxa_asv[[1]], paste0(path_to_taxa, "taxonomy_table_boot_", boot, "_NBC_", clust_method, ".Rds"))
# and bootstrap scores
saveRDS(taxa_asv[[2]], paste0(path_to_taxa, "taxonomy_table_boot_", boot, "_NBC_", clust_method, "_scores.Rds"))

######################### get table ready for latex

# load data
dada_info <- read.csv(paste0(save_data, "dada_filtering_steps.csv"), header=T, sep="\t")
match_table <- read.csv(paste0(save_data, "matching_names_table.csv"), header=T, sep="\t")

# create rna names
all_xnas <- rbind.data.frame(match_table, cbind.data.frame(sampleID=gsub("AD006", "AD012", match_table$sampleID), sample_info=match_table$sample_info, xna=rep("RNA", nrow(match_table))))

# create a sampleID column
dada_info$sampleID <- unlist(lapply(strsplit(rownames(dada_info), "_",), "[[", 1))

# merge by sampleID
final_latex <- merge(dada_info, all_xnas, by="sampleID")

# split col into cols
final_latex$treatment <- unlist(lapply(strsplit(final_latex$sample_info, "_"), "[[", 1))
final_latex$treatment <- ifelse(final_latex$treatment=="C", "Control", "Warming")
final_latex$site <- unlist(lapply(strsplit(final_latex$sample_info, "_"), "[[", 2))
final_latex$timepoint <- unlist(lapply(strsplit(final_latex$sample_info, "_"), "[[", 3))

# sort by site and month and drop colums
final_latex_export <- final_latex[with(final_latex, order(treatment, site, -xtfrm(timepoint))), c("treatment", "timepoint", "site", "input", "filtered", "denoisedF", "denoisedR", "merged")]

colnames(final_latex_export) <- c("Treatment", "Timepoint", "Site", "Input", "Filtered", "Denoised forward", "Denoised reverse", "Merged")

# export
write.table(final_latex_export[which(final_latex$xna=="DNA"), ], paste0(save_data, "dada_info_dna.tex"), col.names=T, row.names=F, sep=" & ", quote=F)

write.table(final_latex_export[which(final_latex$xna=="RNA"), ], paste0(save_data, "dada_info_rna.tex"), col.names=T, row.names=F, sep=" & ", quote=F)

