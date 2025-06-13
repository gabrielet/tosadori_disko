################################################################### IMPORTANT NOTE ##################################################################
# this pipeline is run on a subset of the cutadapted samples. Indeed, only Warming and Control samples were considered for this analysis. If it is
# necessary to add more samples, this script should be run after adding the samples of interest to the directory that is used as input in the first
# qiime step. See below for the first command involving the /microbiology/disko2013/analyses_bacteria/experiments/november/cutadapted_selected
# directory.
# In order to create such directory, the following should be done:

# create the local cutadapted directory, i.e. the one for the experiment that needs to be performed
# for example, we want to perform the experiment Warming versus Control, DNA+RNA.
# so, we create the merged exp directory and we copy Warming and Control samples for both DNA and RNA

echo "copying cutadapted files"

# set some variables
orig_dir_RNA="/microbiology/disko2013/analyses_bacteria/AD012/cutadapted_final/"
orig_dir_DNA="/microbiology/disko2013/analyses_bacteria/AD006/cutadapted_final/"
dest_dir="/microbiology/disko2013/analyses_bacteria/experiments/cutadapted/dada/"

# create destination directory
mkdir -p $dest_dir

# compile scripts to copy samples that are necessary for the analysis
# modify the regexps as needed

# grep C and W only, for R1
ls $orig_dir_RNA | grep T["C|W"]["A|J"] | sed "s|^|cp \\$orig_dir_RNA|" | sed "s|$| \\$dest_dir|" > rna_files.sh

# grep C and W only, for R2
ls $orig_dir_DNA | grep T["C|W"]["A|J"] | sed "s|^|cp \\$orig_dir_DNA|" | sed "s|$| \\$dest_dir|" > dna_files.sh

# run the scripts
chmod +x rna_files.sh dna_files.sh
./rna_files.sh
./dna_files.sh

# remove files that are not needed anymore
rm rna_files.sh dna_files.sh

# DONE!
# AT THIS POINT THE EXPERIMENT'S DIRECTORY CONTAINS ALL THE FASTQs THAT ARE REQUIRED TO RUN A SPECIFIC EXPERIMENT
