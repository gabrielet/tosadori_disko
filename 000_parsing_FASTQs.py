def parsing_fastqs(xna, orgn, save_in):

	# debug only
	#xna = "AD006"
	#xna = "AD012"
	
	# importing required libraries
	import re
	import pandas as pd
	import sys

	# path to where the parsed files will be stored
	root_path="/home/gabriele/work/microbiology/disko2013/"
	wobble_path=root_path+"analyses_"+orgn+"/files_for_analysis/"
	fastq_path=root_path+"fastq/"

	# set paths accordingly to the chosen experiment name
	if xna=="AD006":
		print("detected experiment", xna)
		# AD006 R1
		rone_file="PRI-LANY-Jana-AD006_S1_L001_R1_001.fastq"
		# AD006 R2
		rtwo_file="PRI-LANY-Jana-AD006_S1_L001_R2_001.fastq"
	elif xna=="AD012":
		print("detected experiment", xna)
		# AD012 R1
		rone_file="PRI-LANY-Jana-AD012_S2_L001_R1_001.fastq"
		# AD012 R2
		rtwo_file="PRI-LANY-Jana-AD012_S2_L001_R2_001.fastq"
	else:
		sys.exit("experiment not found!")

	# read file
	forward_table = pd.read_csv(wobble_path+"forwardUnwobbled.csv", sep="\t")
	reverse_table = pd.read_csv(wobble_path+"reverseUnwobbled.csv", sep="\t")

	# get second column and create a whole regex
	forward_primers = "(^"+")|(^".join(forward_table[forward_table.columns[1]])+")"
	reverse_primers = "(^"+")|(^".join(reverse_table[reverse_table.columns[1]])+")"

	print("parsing R1")

	# read through R1
	with open(fastq_path+rone_file) as rone:
		for current_line in rone:
				
				# if the script does not read a header, throw a fatal error!!!!!!!
				if not "@" in current_line:
					sys.exit("header not found, huge problem in R1!")
					
				# search for regex in the next line, so that the script works into
				# the actual sequence
				next_line = next(rone)
				
				# search for pattern, i.e. forward and reverse primer
				# we are not using re.match but re.search simply because, as reported here:
				# https://docs.python.org/3/library/re.html#search-vs-match
				# re.match() checks for a match only at the beginning of the string,
				# while re.search() checks for a match anywhere in the string
				
				# re.match is also faster and primers are, supposedly, found at
				# the beginning of a read
				
				matched_fwd = re.match(forward_primers, next_line)
				matched_rev = re.match(reverse_primers, next_line)

				# if primer matches current line, then we can print the current_line,
				# i.e. the header
				if matched_fwd:
					# get position of matching regex
					post = list(map(bool, matched_fwd.groups())).index(True)
					
					# create references file name
					reference_fname = save_in+forward_table[forward_table.columns[0]][post]+"_fwd_"+xna+".csv"
					# open the reference file in append + read mode ("a+")
					with open(reference_fname, "a+") as sampleN:
						# then, append header and just found matching read
						to_write = current_line.strip().split(" ")[0] + "\t" + current_line.strip().split(" ")[1] + "\t" + next_line.strip() + "\t" + next(rone).strip() + "\t" + next(rone)
						sampleN.write(to_write)
					
				# if match is not found, skip next two lines, i.e. + and quality score
				elif matched_rev:
					# get position of matching regex
					post = list(map(bool, matched_rev.groups())).index(True)
					
					# create references file name
					reference_fname = save_in+reverse_table[reverse_table.columns[0]][post]+"_rev_"+xna+".csv"
					
					# open the reference file in append + read mode ("a+")
					with open(reference_fname, "a+") as sampleN:
						# then, append header and just found matching read
						to_write = current_line.strip().split(" ")[0] + "\t" + current_line.strip().split(" ")[1] + "\t" + next_line.strip() + "\t" + next(rone).strip() + "\t" + next(rone)
						sampleN.write(to_write)

				# if match is not found, skip next two lines, i.e. + and quality score
				else:
					# skip the + and the quality score
					to_skip = next(rone) + next(rone)

	print("parsing R2")

	# read through R1
	with open(fastq_path+rtwo_file) as rtwo:
		for current_line in rtwo:

				# if the script does not read a header, throw a fatal error!!!!!!!
				if not "@" in current_line:
					sys.exit("header not found, huge problem in R2!")

				# search for regex in the next line, so that the script works into
				# the actual sequence
				next_line = next(rtwo)
				
				# search for pattern, i.e. forward and reverse primer
				# we are not using re.match but re.search simply because, as reported here:
				# https://docs.python.org/3/library/re.html#search-vs-match
				# re.match() checks for a match only at the beginning of the string,
				# while re.search() checks for a match anywhere in the string
				
				# re.match is also faster and primers are, supposedly, found at
				# the beginning of a read
				
				matched_fwd = re.match(forward_primers, next_line)
				matched_rev = re.match(reverse_primers, next_line)

				# if primer matches current line, then we can print the current_line,
				# i.e. the header
				if matched_fwd:
					# get position of matching regex
					post = list(map(bool, matched_fwd.groups())).index(True)

					# create references file name
					reference_fname = save_in+forward_table[forward_table.columns[0]][post]+"_fwd_"+xna+".csv"
					# open the reference file in append + read mode ("a+")
					with open(reference_fname, "a+") as sampleN:
						# then, append header and just found matching read
						to_write = current_line.strip().split(" ")[0] + "\t" + current_line.strip().split(" ")[1] + "\t" + next_line.strip() + "\t" + next(rtwo).strip() + "\t" + next(rtwo)
						sampleN.write(to_write)

				# if match is not found, skip next two lines, i.e. + and quality score
				elif matched_rev:
					# get position of matching regex
					post = list(map(bool, matched_rev.groups())).index(True)

					# create references file name
					reference_fname = save_in+reverse_table[reverse_table.columns[0]][post]+"_rev_"+xna+".csv"
					# open the reference file in append + read mode ("a+")
					with open(reference_fname, "a+") as sampleN:
						# then, append header and just found matching read
						to_write = current_line.strip().split(" ")[0] + "\t" + current_line.strip().split(" ")[1] + "\t" + next_line.strip() + "\t" + next(rtwo).strip() + "\t" + next(rtwo)
						sampleN.write(to_write)
						
				# if match is not found, skip next two lines, i.e. + and quality score
				else:
					# skip the + and the quality score
					to_skip = next(rtwo) + next(rtwo)

	# at this point the original fastq files are parsed into several csv tables,
	# one for each primer either forward or reverse
