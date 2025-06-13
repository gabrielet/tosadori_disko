#The data we are using are already demultiplexed and barcodes were already removed. So, the FASTQ files are ready to perform the DADA2 pipeline

# following the pipeline suggested here, for inferring ASVs
# https://benjjneb.github.io/dada2/tutorial.html

# source("005_dada_analysis_fungi.R")

# load libraries
library("dada2")
library("ggplot2")
library("phyloseq")
library("DECIPHER")

# define function for exporting sequences as fasta files
# to blast them in a format qiime likes, or other purposes
# https://docs.qiime2.org/2021.11/tutorials/importing/#sequences-without-quality-information-i-e-fasta
# The ID in each header must follow the format <sample-id>_<seq-id>.
# <sample-id> is the identifier of the sample the sequence belongs to, and
# <seq-id> is an identifier for the sequence within its sample.

make_fastas <- function(sequences) {

	# create fake fasta headers
	snames <- paste(">", sequences, sep="")
	fastas <- vector()

	# loop through all the names
	for (i in c(1:length(snames))) {

		# and put together the table
		fastas <- append(fastas, rbind(snames[i], sequences[i]))
	}

	# results
	return(fastas)
}

# import functions
source("/microbiology/disko2013/code/000_micro_functions_disko2013.R")

##################### select the experiment #####################
exp_name <- "september_2023"

# select the proper organism
orgn <- "fungi"

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn, " EXPERIMENT NAME ", exp_name))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/microbiology/diskodry2021/"
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

if (orgn == "bacteria") {
	# 16S calls for truncation
	out <- filterAndTrim(fwd=fn_fwds, filt=filt_fwds, rev=fn_revs, filt.rev=filt_revs, truncLen=c(220, 220), maxN=0, maxEE=c(4, 4), truncQ=2, rm.phix=FALSE, compress=TRUE, multithread=TRUE)
} else if (orgn == "fungi") {
	# no truncation because ITS has variable length
	out <- filterAndTrim(fwd=fn_fwds, filt=filt_fwds, rev=fn_revs, filt.rev=filt_revs, truncQ=2, minLen=50, maxN=0, maxEE=c(4, 4), rm.phix=FALSE, compress=TRUE, multithread=TRUE)
}

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

# select min bootstrap to assign a taxa. see:
# https://benjjneb.github.io/dada2/assign.html
# At this more stringent threshold, sequence 5 and 7 are not classified
# at the genus level. Instead an NA is returned, indicating that less
# than minBoot=80 of the 100 bootstraps returned the same genus for 
# those sequences.
boot <- 60

# create a DNAStringSet from the ASVs
xna_set <- DNAStringSet(getSequences(seq_tab_no_chim))

# load trainingSet. it is already a variable
load("/microbiology/narsarsuaq/original_metadata/SILVA_SSU_r138_2019.RData")

# run idtaxa
ids <- IdTaxa(xna_set, trainingSet, strand="both", threshold=boot, processors=NULL, verbose=FALSE)

# ranks of interest
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxa_asv <- t(sapply(ids, function(x) {
	m <- match(ranks, x$rank)
	tx_asv <- x$taxon[m]
	return(tx_asv)
}))

# assign col and rownames
colnames(taxa_asv) <- ranks;
rownames(taxa_asv) <- getSequences(seq_tab_no_chim)

# export taxonomy table
saveRDS(taxa_asv, paste0(path_to_taxa, "taxonomy_table_boot_", boot, "_idtaxa_", clust_method, ".Rds"))

###########################################################################################################################
# run assignment for ASVs
boot <- 70

taxa_asv <- assignTaxonomy(seq_tab_no_chim, "/microbiology/databases/silva_nr99_v138.1_train_set.fa.gz", minBoot=boot, outputBootstraps=TRUE, tryRC=TRUE, multithread=TRUE)

saveRDS(taxa_asv[[1]], paste0(path_to_taxa, "taxonomy_table_boot_", boot, "_NBC_", clust_method, ".Rds"))
# and bootstrap scores
saveRDS(taxa_asv[[2]], paste0(path_to_taxa, "taxonomy_table_boot_", boot, "_NBC_", clust_method, "_scores.Rds"))

