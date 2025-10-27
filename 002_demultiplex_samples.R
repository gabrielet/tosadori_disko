# source("002_demultiplex_samples.R")

# IMPORTANT NOTE
# this file is used to create all the fastq files for all the samples
# no matter the treatment.
# this should be the starting point of all the analyses

# import library to run python scripts
library("reticulate")

print("VERY IMPORTANT NOTE!")
print("RUN THIS CODE ONLY AFTER REMOVING ALL EXISTING DEMULTIPLEX/CUTADAPTED FILES")

##################### select the experiment #####################

xna <- "AD006" # DNA
#xna <- "AD012" # RNA

# select the proper organism, hence directory
orgn <- "bacteria"

print(paste0("THE ANALYSIS IS PERFORMED ON ", orgn))

# and set the organism directory accordingly
ifelse(orgn=="fungi", orgn_dir <- "analyses_fungi/", orgn_dir <- "analyses_bacteria/")

# set path according to the experiment
root_path <- "/home/gabriele/work/microbiology/disko2013/"
path_to_exp <- paste0(root_path, orgn_dir, xna, "/")
path_to_refs <- paste0(path_to_exp, "references_fastq/")
path_to_csvs <- paste0(root_path, orgn_dir, "files_for_analysis/")
save_in_demux <- paste0(path_to_exp, "shiny_new_demultiplexed_fastq/")

# create necessary directories if they do not exist
ifelse(!dir.exists(path_to_exp), dir.create(path_to_exp, recursive = T), "dir exists!")
ifelse(!dir.exists(path_to_refs), dir.create(path_to_refs, recursive = T), "dir exists!")
ifelse(!dir.exists(save_in_demux), dir.create(save_in_demux, recursive = T), "dir exists!")

# source python parsing script, i.e. load the python function to parse the fastq files
source_python(paste0(root_path, "code/000_parsing_FASTQs.py"))

# perform the parsing using the loaded python function
pRes <- parsing_fastqs(xna, orgn, path_to_refs)

# perform actual demultiplexing with R

print("demultiplexing")

# get sample info table and load this info. in this table
# we have info regarding the forward and reverse primers that
# were used to identify a specific sample
infoSample <- read.csv(paste0(path_to_csvs, tolower(orgn), "_", xna, ".csv"), header=T, sep="\t")

# for each sample, i.e. each row of infoSample
# what we do is, for each sample get the file of its corresponding forward
# and reverse primers and find the headers that match. matching headers
# represents the R1 and R2 reads which belongs to that specific samples,
# since they matched its primers.
for (smpl in c(1:nrow(infoSample))) {

#	smpl <-1

	# each sample is identified by a couple of primers (fwd, rev). what we need to do
	# is to get reads that starts with only these primers
	
	# load file containing all sequences starting with a sample-specific forward primer only
	fwD <- read.csv(paste0(path_to_refs, infoSample$sequence_of_forward_primer[smpl], "_fwd_", xna, ".csv"), header=F, sep="\t")
	# load file containing all sequences starting with a sample-specific reverse primer only
	reV <- read.csv(paste0(path_to_refs, infoSample$sequence_of_reverse_primer[smpl], "_rev_", xna, ".csv"), header=F, sep="\t")

	# set colnames
	colnames(fwD) = c("header", "direction", "sequence", "plus", "quality")
	colnames(reV) = c("header", "direction", "sequence", "plus", "quality")

	# headers are used to identify a sample. so, for all the sequences starting with
	# a sample-specific couple of primers (fwd, rev), we need to find shared headers
	matchingHeads <- intersect(reV$header, fwD$header)

######	debug only
######	write.table(matchingHeads, "/home/gabriele/work/microbiology/disko2013/analyses_bacteria/AD006/matching_heads.csv", quote=F, row.names=F, col.names=F

	# now get headers-corresponding lines from the two primers files, no matter the
	# direction of the reads. the important thing is that each line is a read
	# identified by one of the shared headers
	
	# get forward primer reads
	fwdFoundAllDir <- fwD[fwD$header %in% matchingHeads, ]
	# get reverse primer reads
	revFoundAllDir <- reV[reV$header %in% matchingHeads, ]

	# now get reads starting with a specific primer
	if (xna=="AD006") {
	
		# the problem here is that the file contaning the reads starting with a specific primer
		# can containg both forward and reverse reads. This is an artifact of the sequenceing
		# what we need to do, then, is to get all the forward and reverse reads that starts with
		# a specific primer:
		# forward primer with forward and reverse reads (i.e. R1 forward primer and R2 forward primer, i.e. forward primer from 1:N:0:1 and 2:N:0:1)
		# reverse primer with forward and reverse reads (i.e. R1 reverse primer and R2 reverse primer, i.e. reverse primer from 1:N:0:1 and 2:N:0:1)
		
		# rbind forward reads (1:N:0:1) starting with forward primer and reverse reads (2:N:0:1) starting with forward primer
		fwdFoundOneWay <- rbind.data.frame(fwdFoundAllDir[which(fwdFoundAllDir$direction=="1:N:0:1"), ], fwdFoundAllDir[which(fwdFoundAllDir$direction=="2:N:0:1"), ])
		# rbind forward reads (1:N:0:1) starting with reverse primer and reverse reads (2:N:0:1) starting with reverse primer
		revFoundOneWay <- rbind.data.frame(revFoundAllDir[which(revFoundAllDir$direction=="1:N:0:1"), ], revFoundAllDir[which(revFoundAllDir$direction=="2:N:0:1"), ])
	} else if (xna=="AD012") {
		# same as above!
		fwdFoundOneWay <- rbind.data.frame(fwdFoundAllDir[which(fwdFoundAllDir$direction=="1:N:0:2"), ], fwdFoundAllDir[which(fwdFoundAllDir$direction=="2:N:0:2"), ])
		revFoundOneWay <- rbind.data.frame(revFoundAllDir[which(revFoundAllDir$direction=="1:N:0:2"), ], revFoundAllDir[which(revFoundAllDir$direction=="2:N:0:2"), ])
	} else {
		print("experiment not found!")
	}

	# at this sort lines by header, in order to have two files which lines are matching
	# this is required by qiime
	fwdFound <- fwdFoundOneWay[order(fwdFoundOneWay$header), ]
	revFound <- revFoundOneWay[order(revFoundOneWay$header), ]

	# now merge header and direction together to prepare for exporting
	fwdMerged <- cbind.data.frame(header=paste(fwdFound$header, fwdFound$direction, sep=" "), fwdFound[, c(3, 4, 5)])
	revMerged <- cbind.data.frame(header=paste(revFound$header, revFound$direction, sep=" "), revFound[, c(3, 4, 5)])

	# get final sample name, in casava 1.8 format as specified here:
	# https://docs.qiime2.org/2021.4/tutorials/importing/#casava-1-8-paired-end-demultiplexed-fastq
	# i.e. sample:name_CD18_L000_R[12]_001.fastq.gz
	sNameROne <- paste0(infoSample$sampleID[smpl], "_CD18_L000_R1_001.fastq")
	sNameRTwo <- paste0(infoSample$sampleID[smpl], "_CD18_L000_R2_001.fastq")

	# now export the two tables as R1 and R2 for a sample!
	write.table(as.vector(t(fwdMerged)), paste0(save_in_demux, sNameROne), row.names=F, col.names=F, quote=F, sep="\t")
	write.table(as.vector(t(revMerged)), paste0(save_in_demux, sNameRTwo), row.names=F, col.names=F, quote=F, sep="\t")
}

# NEXT STEP WILL BE THE REMOVAL OF PRIMERS FROM R1 AND R2 USING CUTADAPT

# FOR PAIRED-ENDS REMEMBER THAT R1 MATCH THE REVERSE COMPLEMENT OF R2

# ALSO, PRIMERS ARE FOUND MULTIPLE TIMES INSIDE A READ (which is weird?)

