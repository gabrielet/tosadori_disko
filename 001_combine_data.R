# load libraries
library("Biostrings")
library("openxlsx")

# set raw_paths
orgn <- "fungi" # can be bacteria
raw_path <- "/microbiology/disko2013/"
raw_path_original <- paste0(raw_path, "metadata/csv_original/")
save_in <- paste0(raw_path, "analyses_", orgn, "/files_for_analysis/")

# create analyses directories if it does not exist
ifelse(dir.exists(save_in), TRUE, dir.create(save_in, recursive=T))

# select one exp between AD006 and AD0012
for (exp_ad in c("AD006", "AD012")) {

	############################
	# OBTAIN CLEAN METADATA TABLES THAT WILL BE LATER MERGED WITH SAMPLE_TO_TAG TABLES

	# load data
	if (exp_ad == "AD006") {
		meta_data <- read.csv(paste0(raw_path_original, "metaDNA.csv"), sep="\t", header=T)
		# remove unwanted characters from colnames
		colnames(meta_data) <- gsub("\\.", "_", colnames(meta_data))

		# remove unwanted cols
		meta_data$DosT <- NULL
		meta_data$DosD <- NULL
		
		# rename collection site columns
		colnames(meta_data)[which(colnames(meta_data)=="Fence")] <- "CollectionSite"

		# set clearer collection sites
		meta_data$CollectionSite <- paste0("s", meta_data$CollectionSite)
		
		# remove existing sampleID columns and other useless cols
		meta_data$X_SampleID <- NULL
		meta_data$SampleID <- NULL
		meta_data$BarcodeSequence <- NULL
		meta_data$LinkerPrimerSequence <- NULL
		meta_data$ReversePrimer <- NULL
	} else if (exp_ad == "AD012") {
		meta_data <- read.csv(paste0(raw_path_original, "metaRNA.csv"), sep="\t", header=T)
		# remove unwanted characters from colnames
		colnames(meta_data) <- gsub("\\.", "_", colnames(meta_data))
		meta_data$new_sample_ID <- NULL
		
		# rename collection site columns
		colnames(meta_data)[which(colnames(meta_data)=="Description")] <- "CollectionSite"
		
		# set clearer collection sites
		meta_data$CollectionSite <- paste0("s", gsub("R", "", meta_data$CollectionSite))
		
		# remove existing sampleID columns and other useless cols
		meta_data$Sample_ID <- NULL
		meta_data$BarcodeSequence <- NULL
		meta_data$LinkerPrimerSequence <- NULL
		meta_data$ReversePrimer <- NULL
		# and rename weird column names
		names(meta_data)[names(meta_data) == 'Chao1___2420seq_over_sample'] <- 'Chao1_2420_seq_over_sample'
	} else {
		print("the experiment you selected does not exist!")
	}

	# rename time point columns
	colnames(meta_data)[which(colnames(meta_data)=="DOS")] <- "TimePoint"

	# set clearer time points
	time_p <- c("June", "July", "August")

	meta_data$TimePoint[which(meta_data$TimePoint=="T1")] <- time_p[1]
	meta_data$TimePoint[which(meta_data$TimePoint=="T2")] <- time_p[2]
	meta_data$TimePoint[which(meta_data$TimePoint=="T3")] <- time_p[3]

	# create a new sampleID column
	merged_ids <- paste(meta_data$Treatment, meta_data$CollectionSite, meta_data$TimePoint, sep="_")

	# finally, build the tables to be exported
	final_clean <- cbind.data.frame(merged_ids, meta_data)

	# substitute zeroes with NAs since these 0 represent missing data
	# find zeroes
	zeroes <- which(final_clean==0, arr.ind=T)
	# substitute all of them
	for (i in 1:nrow(zeroes)) {
		final_clean[zeroes[i, 1], zeroes[i, 2]] <- NA
	}

	############################
	# OBTAIN CLEAN SAMPLE_TO_TAG TABLE

	# load original table
	original_treat <- read.csv(paste0(raw_path_original, "samples_treatments.csv"), header=T, sep="\t")

	# at this point get necessary info from the tag_to_sample mapping table
	map_tab <- read.csv(paste0(raw_path_original, "sample_to_tag_mapping.csv"), sep="\t", header=T)

	# remove unwanted character from colnames
	colnames(map_tab) <- gsub("\\.", "_", colnames(map_tab))

	# subset the table for bacteria
	ifelse(orgn=="bacteria", orgn <- "Bacteria", orgn <- "Fungi")
	map_tab_bacteria <- map_tab[map_tab$Bacteria_Fungi==orgn, ]

	# rename sample_name column to avoid qiime conflicts
	names(map_tab_bacteria)[names(map_tab_bacteria)=="sample_name"] <- "sample"

	# list possible treatments
	treats <- c("S", "W", "R")

	# zeroes to false, ones to true
	original_treat[, c(3,4,5)] <- ifelse(original_treat[, c(3,4,5)]==0, FALSE, TRUE)
	clean_treats <- vector()
	for (ln in c(1:nrow(original_treat))) {
		# get controls in form of string
		string_t <- paste(treats[unlist(original_treat[ln, c(3,4,5)])], collapse="")
		# then put everything together. if string_t is empty all the conditions are
		# FALSE then it's a control WITHOUT snow addition
		if (nchar(string_t) == 0) {
			clean_treats <- append(clean_treats, "C")
		}
		# else it is a control WITH snow addition
		else if (string_t == "S"){
			clean_treats <- append(clean_treats, "SC")
		}
		# else, it is treated somehow
		else {
			clean_treats <- append(clean_treats, string_t)
		}
	}

	# merge info with sample_to_tag_mapping.csv, before cleaning the table itself
	bact_subset <- map_tab_bacteria[map_tab_bacteria$ligation_tag==exp_ad, ]

	# prepare to add the treatment
	bact_subset <- cbind.data.frame(bact_subset, treatment=rep(0, nrow(bact_subset)))

	# for each time point
	for (tp in c(1, 2, 3)) {
		# get position for a specific time point
		pos <- which(bact_subset$time==tp)
		# set month
		bact_subset$time[pos] <- time_p[tp]
		# append experimental information to the mapping table
		bact_subset$treatment[pos] <- clean_treats
	}

	# ok, the remove missing samples, i.e. "no sample"
	bacteria_clean <- bact_subset[!bact_subset$sample=="no sample", ]

	# if notes is empty, remove the whole column
	if (all(bacteria_clean$notes == "")) {
		bacteria_clean$notes <- NULL
	}

	# re-arrange tables
	export_this <- cbind.data.frame(bacteria_clean[, c(1,2)], treatment=bacteria_clean[, ncol(bacteria_clean)], bacteria_clean[, c(3,4,5,6,7,8,10)])

	# at this point, add the combined sample names that will be used as:
	# filenames for cutadapted reads
	# sampleIDs for qiime and other referencing, without any symbol
	cuta_subset <- apply(export_this[, c(1,3,2,6,7)] , 1, paste, collapse="")
	export_this <- cbind.data.frame(sampleID=cuta_subset, export_this)

	############################
	# FINALISE THE MERGING BETWEEN SAMPLE_TO_TAG AND METADATA TABLES

	# create the merged_ids column to cross-reference this table with the metadata table
	# to achieve this goal we need to find the sampling site using the sample name

	# create 8 sub-sites for each site and paste them with its sample name
	sites_by_samples <- cbind.data.frame(site=paste0("s", rep(c(1:6), each=8)), sample=paste0(c(1:48), "T"))

	# for each sample, retrieve its site
	matched_site <- vector()
	# loop through the samples
	for (smpl in c(1:nrow(export_this))) {
		# loop through the sites
		for (st in c(1:nrow(sites_by_samples))) {
			# is sample and sample match
			if (export_this$sample[smpl] == sites_by_samples$sample[st]) {
				# get the site
				matched_site <- append(matched_site, sites_by_samples$site[st])
			}
		}
	}

	# create the merged_ids column
	merged_ids_subset <- paste(export_this$treatment, matched_site, export_this$time, sep="_")

	# adding a column specifying both time and treatment
	treat_by_time <- paste0(export_this$treatment, "_", export_this$time)

	# put everything together. the first columns HAVE TO be sampleID
	# as required by qiime
	export_this <- cbind.data.frame(sampleID=export_this$sampleID, merged_ids=merged_ids_subset, TreatBYTime=treat_by_time, export_this[, c(2:ncol(export_this))])

	# export the tables that are used to perform the demultiplexing
	write.table(export_this, paste0(save_in, orgn, "_", exp_ad, ".csv"), row.names=F, quote=F, sep="\t")

	# cross reference metadata and sample_to_tag. some columns are duplicated but useful to check if they match
	finalmeta_data_subset <- merge(export_this, final_clean, by="merged_ids")

	# put everything together. the first columns HAVE TO be sampleID
	# as required by qiime
	finalmeta_data_subset <- cbind.data.frame(sampleID=finalmeta_data_subset$sampleID, finalmeta_data_subset[, -2])

	# remove unwanted columns from final metadata table
	finalmeta_data_subset$time <- NULL
	finalmeta_data_subset$treatment <- NULL
	finalmeta_data_subset$name_of_forward_primer <- NULL
	finalmeta_data_subset$sequence_of_forward_primer <- NULL
	finalmeta_data_subset$sequence_of_reverse_primer <- NULL

	# export the merged tables
	write.table(finalmeta_data_subset, paste0(save_in, "metadata_disko_", exp_ad, ".csv"), row.names=F, quote=F, sep="\t")
	############################
	# OBTAIN CLEAN-FROM-WOBBLE-BASES PRIMERS TABLES

	# create un-wobbled tables by removing wobble bases and substituting them with a regular expression

	# import un-wobbling, manually curated table
	wobble <- read.csv(paste0(raw_path_original, "sub_wobble_bases.csv"), sep="\t", header=T)

	# get some variables ready
	unique_fwd <- unique(bacteria_clean$sequence_of_forward_primer)
	unique_rev <- unique(bacteria_clean$sequence_of_reverse_primer)

	forward_unwobbled <- cbind.data.frame(wobbled = unique_fwd, unwobbled = unique_fwd)
	reverse_unwobbled <- cbind.data.frame(wobbled = unique_rev, unwobbled = unique_rev)

	# and parse each wobble basis with a regexp
	for (linE in 1:nrow(wobble)) {
		forward_unwobbled$unwobbled <- gsub(wobble$wobble[linE], wobble$bio[linE], forward_unwobbled$unwobbled)
		reverse_unwobbled$unwobbled <- gsub(wobble$wobble[linE], wobble$bio[linE], reverse_unwobbled$unwobbled)
	}

	# prepare for reverseComplementing
	rev_comp_forward_unwobbled <- vector()
	rev_comp_reverse_unwobbled <- vector()
	reverse_forward_unwobbled <- vector()
	reverse_reverse_unwobbled <- vector()
	rev_comp_forward_wobbled <- vector()
	rev_comp_reverse_wobbled <- vector()
	reverse_forward_wobbled <- vector()
	reverse_reverse_wobbled <- vector()

	# find reverse and reverseComplement of each forward primer
	for (fwd in c(1:length(forward_unwobbled$unwobbled))) {
		rev_comp_forward_unwobbled <- append(rev_comp_forward_unwobbled, as.character(reverseComplement(DNAString(forward_unwobbled$wobbled[fwd]))))
		reverse_forward_unwobbled <- append(reverse_forward_unwobbled, as.character(reverse(DNAString(forward_unwobbled$wobbled[fwd]))))
		rev_comp_forward_wobbled <- append(rev_comp_forward_wobbled, as.character(reverseComplement(DNAString(forward_unwobbled$wobbled[fwd]))))
		reverse_forward_wobbled <- append(reverse_forward_wobbled, as.character(reverse(DNAString(forward_unwobbled$wobbled[fwd]))))
	}

	# find reverse and reverseComplement of each reverse primer
	for (fwd in c(1:length(reverse_unwobbled$unwobbled))) {
		rev_comp_reverse_unwobbled <- append(rev_comp_reverse_unwobbled, as.character(reverseComplement(DNAString(reverse_unwobbled$wobbled[fwd]))))
		reverse_reverse_unwobbled <- append(reverse_reverse_unwobbled, as.character(reverse(DNAString(reverse_unwobbled$wobbled[fwd]))))
		rev_comp_reverse_wobbled <- append(rev_comp_reverse_wobbled, as.character(reverseComplement(DNAString(reverse_unwobbled$wobbled[fwd]))))
		reverse_reverse_wobbled <- append(reverse_reverse_wobbled, as.character(reverse(DNAString(reverse_unwobbled$wobbled[fwd]))))
	}

	# finalise the tables
	forward_unwobbled_final <- cbind.data.frame(forward_unwobbled, revUnwobbled=reverse_forward_unwobbled, revCompUnwobbled=rev_comp_forward_wobbled, rev=reverse_forward_wobbled, revComp=rev_comp_forward_wobbled)
	reverse_unwobbled_final <- cbind.data.frame(reverse_unwobbled, revUnwobbled=reverse_reverse_unwobbled, revCompUnwobbled=rev_comp_reverse_wobbled, rev=reverse_reverse_wobbled, revComp=rev_comp_reverse_wobbled)

	# and substitute the wobble bases in the reverseComplement-ed primers
	# and parse each wobble basis with a regexp
	for (linE in 1:nrow(wobble)) {
		forward_unwobbled_final$revCompUnwobbled <- gsub(wobble$wobble[linE], wobble$bio[linE], forward_unwobbled_final$revCompUnwobbled)
		reverse_unwobbled_final$revCompUnwobbled <- gsub(wobble$wobble[linE], wobble$bio[linE], reverse_unwobbled_final$revCompUnwobbled)
		forward_unwobbled_final$revUnwobbled <- gsub(wobble$wobble[linE], wobble$bio[linE], forward_unwobbled_final$revUnwobbled)
		reverse_unwobbled_final$revUnwobbled <- gsub(wobble$wobble[linE], wobble$bio[linE], reverse_unwobbled_final$revUnwobbled)
	}

	# export tables
	write.table(forward_unwobbled_final, paste0(save_in, "forward_unwobbled.csv"), sep="\t", row.names = F, quote = F)
	write.table(reverse_unwobbled_final, paste0(save_in, "reverse_unwobbled.csv"), sep="\t", row.names = F, quote = F)

	# these tables can be used to create the multi-primer strings to remove
	# multiple primers and their reverseComplement that are found all across
	# across the reads. Also, they are used at the end of the cutadapt trim
	# step to test whether primers are still found in the trimmed sequences
}
