# source("003_run_cutadapt.R")

# load libraries
library("Biostrings")

##################### select the experiment #####################
#exp_name <- "AD006" # DNA
exp_name <- "AD012" # RNA

# select the proper organism, hence directory
orgn <- "bacteria"

# set proper reverse complement primers
if (orgn == "bacteria") {

	# define primers
	fwd_orig <- "GTGCCAGCMGCCGCGGTAA"
	rev_orig <- "GGACTACHVGGGTWTCTAAT"
	
	# set forward and reverse, reverse complement
	fwd_rev_comp <- as.character(reverseComplement(DNAString(fwd_orig)))
	rev_rev_comp <- as.character(reverseComplement(DNAString(rev_orig)))
} else {

	# define primers
	fwd_orig <- "GTGARTCATCGARTCTTTG"
	rev_orig <- "TCCTCCGCTTATTGATATGC"
	
	# set forward and reverse, reverse complement
	fwd_rev_comp <- as.character(reverseComplement(DNAString(fwd_orig)))
	rev_rev_comp <- as.character(reverseComplement(DNAString(rev_orig)))
}

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/mnt/cinqueg/gabriele/work/microbiology/disko2013/"
path_to_exp <- paste0(root_path, orgn_dir, exp_name, "/")
fastq_path <- paste0(path_to_exp, "shiny_new_demultiplexed_fastq/")
cuta_one <- paste0(path_to_exp, "cutadapted_step_one/")
cuta_two <- paste0(path_to_exp, "cutadapted_step_two/")
cuta_three <- paste0(path_to_exp, "cutadapted_step_three/")
cuta_final <- paste0(path_to_exp, "cutadapted_final/")

# create necessary directories if they do not exist
ifelse(!dir.exists(cuta_one), dir.create(cuta_one, recursive = T), "dir exists!")
ifelse(!dir.exists(cuta_two), dir.create(cuta_two, recursive = T), "dir exists!")
ifelse(!dir.exists(cuta_three), dir.create(cuta_three, recursive = T), "dir exists!")
ifelse(!dir.exists(cuta_final), dir.create(cuta_final, recursive = T), "dir exists!")

# get sample info table and load this info
infoSample <- read.csv(paste0(root_path, orgn_dir, "files_for_analysis/", orgn, "_", exp_name, ".csv"), header=T, sep="\t")

# call cutadapt
cutadapt <- "/home/gabriele/.local/bin/cutadapt"
system2(cutadapt, args = "--version")

# step one: remove forward and reverse from R1 and R2
# using regular 5' adapter (see above) because we expect
# primers to be found at the beginning of the reads: so,
# we want whatever was sequenced AFTER it.
for (smpl in c(1:nrow(infoSample))) {

	# get samples name
	sNameROne <- paste0(fastq_path, infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq")
	sNameRTwo <- paste0(fastq_path, infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq")

	# test if both R1 and R2 files exist
	if (file.exists(sNameROne) & file.exists(sNameRTwo)) {
	
		print(paste0("STEP ONE - remove primers from sample ", smpl))
		
		# create filenames for the cutadapt-step-one outputs
		fwdCut <- paste0(cuta_one, paste0(infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq"))
		revCut <- paste0(cuta_one, paste0(infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq"))

		# primers can be anywhere so don't define specific rules for
		# their position, i.e. by attaching "^" or something
		fwdPrimer <- infoSample$sequence_of_forward_primer[smpl]
		revPrimer <- infoSample$sequence_of_reverse_primer[smpl]
		
		# THIS IS HOW CUTADAPT WORKS
		# see here for references
		# https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-read
		# actual primer removal
		# first, specify which primer to remove from R1 and R2
		# second, specify the output
		# third, specify the input
		# fourth, tell cutadapt not to print anything: it is faster!
		system2(cutadapt, args = c("-g", fwdPrimer, "-n", 2, "-G", revPrimer, "-n", 2, "--minimum-length", 200, "-o", fwdCut, "-p", revCut, sNameROne, sNameRTwo, "--quiet"))
	}
}

# step two: remove remaining primers, if any
# is still inside the sequence. Again, treat them as
# regular 5' adapters
for (smpl in c(1:nrow(infoSample))) {

	# get samples name from step two
	sNameROne <- paste0(cuta_one, infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq")
	sNameRTwo <- paste0(cuta_one, infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq")

	# test if both R1 and R2 files exist
	if (file.exists(sNameROne) & file.exists(sNameRTwo)) {

		print(paste0("STEP TWO - remove primers from sample ", smpl))

		# create filenames for the final-cutadapt-step outputs
		fwdCut <- paste0(cuta_two, paste0(infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq"))
		revCut <- paste0(cuta_two, paste0(infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq"))

		# remove remaining primers, BUT USING THE PRIMERS ONLY
		# WITHOUT THE ADAPTERS
		system2(cutadapt, args = c("-g", fwd_orig, "-n", 2, "-G", rev_orig, "-n", 2, "--minimum-length", 200, "-o", fwdCut, "-p", revCut, sNameROne, sNameRTwo, "--quiet"))
	}
}

# step three: remove remaining primers, if any
# this time in the opposite reads, i.e. forward in R2
# reverse in R1. Again, treat them as regular 5'
for (smpl in c(1:nrow(infoSample))) {

	# get samples name from step two
	sNameROne <- paste0(cuta_two, infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq")
	sNameRTwo <- paste0(cuta_two, infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq")

	# test if both R1 and R2 files exist
	if (file.exists(sNameROne) & file.exists(sNameRTwo)) {

		print(paste0("STEP THREE - remove primers from sample ", smpl))

		# create filenames for the final-cutadapt-step outputs
		fwdCut <- paste0(cuta_three, paste0(infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq"))
		revCut <- paste0(cuta_three, paste0(infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq"))

		# remove remaining primers, BUT USING THE PRIMERS ONLY
		# WITHOUT THE ADAPTERS.
		# here we search reverse primer in R1 and forward
		# primer in the R2
		system2(cutadapt, args = c("-g", rev_orig, "-n", 2, "-G", fwd_orig, "-n", 2, "--minimum-length", 200, "-o", fwdCut, "-p", revCut, sNameROne, sNameRTwo, "--quiet"))
	}
}

# step four: remove reverseComplement primer, if any
# in this case, however, the reverse complement is treated
# as a regular 3' adapter (see above) since it is expected
# to ligate at the end of the amplicon: we want whatever
# was sequenced BEFORE it
for (smpl in c(1:nrow(infoSample))) {

	# get samples name from step two
	sNameROne <- paste0(cuta_three, infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq")
	sNameRTwo <- paste0(cuta_three, infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq")

	# test if both R1 and R2 files exist
	if (file.exists(sNameROne) & file.exists(sNameRTwo)) {

		print(paste0("STEP FOUR - remove reverse complement from sample ", smpl))

		# create filenames for the final-cutadapt-step outputs
		fwdCut <- paste0(cuta_final, paste0(infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq"))
		revCut <- paste0(cuta_final, paste0(infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq"))

		# remove any reverseComplement primer
		# system2(cutadapt, args = c("-a", rev_rev_comp, "-n", 2, "-A", rev_rev_comp, "-n", 2, "-a", fwd_rev_comp, "-n", 2, "-A", fwd_rev_comp, "-n", 2, "--minimum-length", 200, "-o", fwdCut, "-p", revCut, sNameROne, sNameRTwo, "--quiet"))
		
		# remove sequences containing a reverseComplement primer
		# since they might be chimera
		system2(cutadapt, args = c("-a", rev_rev_comp, "-A", rev_rev_comp, "-a", fwd_rev_comp, "-A", fwd_rev_comp, "--discard-trimmed", "-o", fwdCut, "-p", revCut, sNameROne, sNameRTwo, "--quiet"))
	}
}

# print out final warning!
hstgs <- paste(rep("#", nchar(paste("gzip *", cuta_final, "*"))), collapse="")

print(hstgs)
print("AT THIS POINT IT IS NECESSARY TO GZIP ALL THE CUTADAPTED FASTQ FILES BY RUNNING:")
print(paste0("gzip ", cuta_final, "*"))
print(hstgs)

