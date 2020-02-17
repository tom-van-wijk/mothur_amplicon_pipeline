#!/usr/bin/env python


# Name:		amplicon_pipeline.py
# Author:	Tom van Wijk
# Date:		16-01-2018
# Licence:	GNU General Public License v3.0 (copy provided in directory)

# For detailed information and instruction on how to install and use this software
# please view the included "README.md" file in this repository


# import python libraries
from argparse import ArgumentParser
from time import gmtime, strftime
import os
import sys
import logging
import logging.handlers
import random
import re


# Function to parse the command-line arguments
# Returns a namespace with argument keys and values
def parse_arguments(args, log):
	log.info("Parsing command line arguments...")
	parser = ArgumentParser(prog="mothur_amplicon_pipeline.py")
	parser.add_argument("-i", "--indir", dest = "input_dir",
		action = "store", default = None, type = str,
		help = "Location of input directory (required)",
		required = True)
	parser.add_argument("-o", "--outdir", dest = "output_dir",
		action = "store", default = 'inputdir', type = str,
		help = "Location of output directory (default=inputdir)")
	parser.add_argument("-t", "--threads", dest = "threads",
		action = "store", default = 8, type = int,
		help = "Number of threads to be used (default=8)")
	parser.add_argument("-a", "--amplicon", dest = "amplicon",
		action = "store", default = '16sv4', type = str,
		help = "Name of the amplicon (default=16sv4)")
	parser.add_argument("-x", "--savetemp", dest = "save_temp",
		action = "store", default = "false", type = str,
		help = "Option to save temporary files (default=false)")
	return parser.parse_args()


# Function creates logger with handlers for both logfile and console output
# Returns logger
def create_logger(logid):
	# create logger
	log = logging.getLogger()
	log.setLevel(logging.INFO)
	# create file handler
	fh = logging.FileHandler(str(logid)+'_mothur_amplicon_pipeline.log')
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(logging.Formatter('%(message)s'))
	log.addHandler(fh)
	# create console handler
	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(logging.Formatter('%(message)s'))
	log.addHandler(ch)
	return log


# Create output directories and generate a project name based on parameter and timestamp
# Returns the filename and output directory
def create_output_dirs(input_dir, output_dir, user_input, log):
	# check if output_dir parameter was given, if so: create output directory, if not:
	# create subdirectory in input_dir for output
	if output_dir == "inputdir":
		output_dir = input_dir
	else:
		os.system("mkdir -p "+output_dir)
	if user_input != 'timestamp':
		project_name = re.sub(r'[^0-9a-zA-Z_-]+', '', user_input)+"_"+strftime("%Y%m%d%H%M%S", gmtime())
	else:
		project_name = strftime("%Y%m%d%H%M%S", gmtime())
	log.info("Project name: "+project_name)
	output_directory = output_dir+"/"+project_name
	return project_name, output_directory


# Function to validate and prepare the following input data for mothur
# This function validates the sequence data files, samplesheet (generates samplesheet
# is nessesary) and filteres all characters from filenames that might cause problems in mothur.
# Returns the name of the final samplesheet, terminates when data is not sufficient or valid.
def organise_input_data(input_dir, output_dir, amplicon, log):
	log.info("Validating input data...")
	# extract compressed fastq files and move regular fastq files and samplesheet
	# (if extists) to from <input_dir> to <output_dir>/temp
	for file in list_directory(input_dir, 'files', 1):
		if file.endswith('.fastq.gz'):
			os.system('gunzip -c '+input_dir+"/"+file+" > "+output_dir+"/"+file.replace('.gz',''))
		elif file.endswith(('.fastq','.samplesheet')):
			os.system('cp '+input_dir+"/"+file+" "+output_dir+"/"+file)
	# create filelist of files in <output_dir>
	filelist = list_directory(output_dir, 'files', 1)
	# rename files to avoid problematic characters for mothur
	for file in filelist:
		if file != re.sub(r'[^\.0-9a-zA-Z_]+', '_', file):
			os.system('mv '+output_dir+'/'+file+" "+output_dir+'/'+re.sub(r'[^\.0-9a-zA-Z_]+', '_', file))
	# check if valid samplesheet if present
	samplesheet_found = find_samplesheet(input_dir, output_dir, log)
	if samplesheet_found == False:
		# If no valid samplesheet was found, once will be automatically generated
		create_sample_sheet(output_dir, log)


# Find samplesheet in input directory, if found, validates samplesheet.
# Returns booloan (True if validated samplesheet is present) and samplesheet name.
def find_samplesheet(input_dir, output_dir, log):
	fastq_files, samplesheets = [], []
	# create list of both fastq files and samplesheets
	for file in list_directory(input_dir, 'files', 1):
		if file.endswith('.samplesheet'):
			samplesheets.append(file)
		elif file.endswith('.fastq'):
			fastq_files.append(file)
	if len(samplesheets) == 0:
		log.info("No samplesheet found, script will generate samplesheet")
		return False
	elif len(samplesheets) > 1:
		log.warning("Multiple samplesheets found, only a single samplesheet is allowed.Script will generate samplesheet")
		return False
	# if a single samplesheet found, validate
	elif len(samplesheets) == 1:
		return validate_samplesheet(samplesheets[0], input_dir, output_dir, fastq_files, log)


# Function to correct and validate samplesheet with data files available in input directory.
# Returns booloan (True if validated samplesheet is present) and samplesheet name.
def validate_samplesheet(samplesheet, input_dir, output_dir, fastq_files, log):
	log.info("Samplesheet found: '"+samplesheet+"'.Validating samplesheet...")
	valid_filepairs = []
	# create a 'corrected' samplesheet that includes the fastq name changes made and eliminates faulthy lines
	with open(input_dir+"/"+samplesheet,"r") as inputfile, open(output_dir+"/mothur.samplesheet", "w") as outputfile:
		for line in inputfile:
			line_segments = line.replace('\n','').replace('\r\n','').split("\t")
			if len(line_segments) < 3:
				log.warning("Invalid line in "+samplesheet+": '"+line+"'Missing columns")
			else:
				header = re.sub(r'[^\.0-9a-zA-Z_]+', '_', line_segments[0])
				file1 = re.sub(r'[^\.0-9a-zA-Z_]+', '_', line_segments[1])
				file2 = re.sub(r'[^\.0-9a-zA-Z_]+', '_', line_segments[2])
				# Checks is files in samplesheet are present and correctly provided with "_R1"/"_R2" flags.
				# Only filepairs that meet these requirements are included in the corrected samplesheet.
				if file1 not in fastq_files:
					log.warning("Invalid line in "+samplesheet+": File: '"+file1+"'not found")
				elif file2 not in fastq_files:
					log.warning("Invalid line in "+samplesheet+": File: '"+file2+"'not found")
				else:
					if "_R1" not in file1:
						log.warning("Invalid line in "+samplesheet+": File: '"+file1+"'missing '_R1' flag")
					elif "_R2" not in file2:
						log.warning("Invalid line in "+samplesheet+": File: '"+file2+"'missing '_R2' flag")
					else:
						valid_filepairs.append(file1+"\t"+file2)
						outputfile.write(header+"\t"+file1+"\t"+file2+"\n")
	inputfile.close(), outputfile.close()
	# Return name of the corrected samplesheet is a minimum 1 filepair is valid.
	if len(valid_filepairs) < 1:
		log.warning("No valid filepairs found in "+samplesheet+"Script will try to generate samplesheet")
		os.system("rm "+output_dir+"/mothur.samplesheet")
		return False
	else:
		message = "Correct filepairs found based on samplesheet: '"+samplesheet+"'.Files included:"
		for p in valid_filepairs:
			message+=("\n"+p)
		log.info(message)
		return True


# Function generates samplesheet that can be used by mothur for processing the fastq filepairs.
# Returns the name of the generated samplesheet or terminates when no valid filepairs are found.
def create_sample_sheet(output_dir, log):
	# Create temporary files with a list of the forward and reverse files
	os.system("ls "+output_dir+"/*_R1*.fastq > "+output_dir+"/r1.txt")
	os.system("ls "+output_dir+"/*_R2*.fastq > "+output_dir+"/r2.txt")
	valid_filepairs = []
	with open(output_dir+"/r1.txt","r") as r1file, open(output_dir+"/r2.txt", "r") as r2file, open(output_dir+"/mothur.samplesheet", "w") as samplesheet:
		rownumber = 1
		r1dir = {}
		for r1line in r1file:
			r1dir[rownumber] = r1line.replace("\n","").split("/")[-1]
			rownumber += 1
		rownumber = 1
		for r2line in r2file:
			# If a reverse file can be correctly matched to the forward file, the filepair
			# is added to the samplesheet and list of valid filepairs
			if str(r1dir.get(rownumber)).replace('_R1','_R2') == r2line.replace("\n","").split("/")[-1]:
				samplesheet.write(str(r1dir.get(rownumber).split("_R1")[0])+"\t"+str(r1dir.get(rownumber))+"\t"+r2line.replace("\n","").split("/")[-1]+"\n")
				valid_filepairs.append(str(r1dir.get(rownumber))+"\t"+r2line.replace("\n","").split("/")[-1])
			rownumber += 1
	r1file.close(), r2file.close(), samplesheet.close()
	# Return name of the corrected samplesheet is a minimum 1 filepair is valid
	# Terminates script if no valid filepairs are found
	if len(valid_filepairs) < 1:
		log.critical("No valid filepairs found.")
	else:
		message = "Correct filepairs included in generated samplesheet:"
		for p in valid_filepairs:
			message+=("\n"+p)
		log.info(message)


# Function creates a list of files or directories in <inputdir>
# on the specified directory depth
def list_directory(input_dir, obj_type, depth):
	dir_depth = 1
	for root, dirs, files in os.walk(input_dir):
		if dir_depth == depth:
			if obj_type ==  'files':
				return files
			elif obj_type == 'dirs':
				return dirs
		dir_depth += 1


# Function to trim reads and filter ambiguous/orphan seqeunces from the input data,
# Also creates fastqc and multiqc quality report of trimmed, filtered files.
# input: directory with target files, number of threads and logger
def trimm_and_filter(directory, threads, log):
	with open(directory+"/mothur.samplesheet","r") as samplesheet_file:
		for line in samplesheet_file:
			# Read trimming
			R1_file = line.split("\t")[1]
			R2_file = line.split("\t")[2].replace("\n", "").replace("\r", "")
			log.info("Trimming filepair:\t"+R1_file+",\t"+R2_file)
			os.system("erne-filter --query1 "+directory+"/"+R1_file+" --query2 "+directory+"/"+R2_file
					+" --output-prefix "+directory+"/"+(R1_file.split("_R1")[0]+"_trimmed")
					+" --threads "+str(threads)+" --q 25 | tee "+directory+"/"+(R1_file.split("_R1")[0])+"_erne.log")
			os.system("rm "+directory+"/"+R1_file+" "+directory+"/"+R2_file+" "+directory+"/"
					  +R1_file.split("_R1")[0]+"_trimmed_unpaired.fastq")
			# Filtering ambigious sequences
			log.info("Filtering ambigious sequences from:\t"+R1_file)
			os.system("fastq_screener.py -i "+directory+"/"+(R1_file.split("_R1")[0]+"_trimmed_1.fastq")
					  +" | tee "+directory+"/"+(R1_file.split("_R1")[0])+"_R1_filter.log")
			os.system("rm "+directory+"/"+(R1_file.split("_R1")[0]+"_trimmed_1.fastq"))
			log.info("Filtering ambigious sequences from:\t"+R2_file)
			os.system("fastq_screener.py -i "+directory+"/"+(R2_file.split("_R2")[0]+"_trimmed_2.fastq")
					  +" | tee "+directory+"/"+(R2_file.split("_R2")[0])+"_R2_filter.log")	
			os.system("rm "+directory+"/"+(R2_file.split("_R2")[0]+"_trimmed_2.fastq"))
			# Removing orphaned sequences
			log.info("Removing orphaned sequences from:\t"+R1_file+",\t"+R2_file)
			os.system("fastq_pair_mapper.py "+directory+"/"+(R1_file.split("_R1")[0]+"_trimmed_1_dada2_filtered.fastq ")
					 +directory+"/"+(R2_file.split("_R2")[0]+"_trimmed_2_dada2_filtered.fastq "))
			os.system("rm "+directory+"/"+(R1_file.split("_R1")[0]+"_trimmed_1_dada2_filtered.fastq "))
			os.system("rm "+directory+"/"+(R2_file.split("_R2")[0]+"_trimmed_2_dada2_filtered.fastq "))
			os.system("rm "+directory+"/"+"*singles.fastq")
			# renaming final files to original names
			os.system("mv "+directory+"/"+(R1_file.split("_R1")[0]+"_trimmed_1_dada2_filtered.fastq_pairs_R1.fastq ")
					 +directory+"/"+R1_file)
			os.system("mv "+directory+"/"+(R2_file.split("_R2")[0]+"_trimmed_2_dada2_filtered.fastq_pairs_R2.fastq ")
					 +directory+"/"+R2_file)
			# run fastqc
			os.system("fastqc "+directory+"/"+R1_file+" "+directory+"/"+R2_file)
	samplesheet_file.close()
	# moving logfiles and quality reports to correct directory
	os.system ("mv "+directory+"/*erne.log "+(directory.replace("/temp", "/logfiles/quality_trimming/")))
	os.system ("mv "+directory+"/*filter.log "+(directory.replace("/temp", "/logfiles/filtering/")))
	os.system ("mv "+directory+"/*fastqc.html "+(directory.replace("/temp", "/quality_control/fastqc/")))
	os.system ("mv "+directory+"/*fastqc.zip "+(directory.replace("/temp", "/quality_control/fastqc/")))
	# unzipping fastqc reports
	files = list_directory(directory.replace("/temp", "/quality_control/fastqc/"), "files", 1)
	for file in files:
		if file.endswith(".zip"):
			os.system("unzip "+(directory.replace("/temp", "/quality_control/fastqc/"))+file
					 +" -d "+(directory.replace("/temp", "/quality_control/fastqc/"))+file.replace('.zip', ''))
	# create multiqc report
	os.system("multiqc "+(directory.replace("/temp", "/quality_control/fastqc/ -o "))
			  +(directory.replace("/temp", "/quality_control/multiqc/")))
	


# Funtion that creates, personalises and executes the mothur script
# to further process the amplicon data and classify the contigs
def run_mothur_classify(workdir, amplicon_name, threads, log):
	# These dictionaries contain the specific mothur parameters for each amplicon
	min_seq_size = {"16sv4":230, "16sv34":415}
	max_seq_size = {"16sv4":280, "16sv34":475}
	alignment_reference_database = {"16sv4":"silva_v4.fasta", "16sv34":"silva_v4.fasta"}
	classification_reference_database_fasta = {"16sv4":"trainset16_022016.pds.fasta", "16sv34":"trainset16_022016.pds.fasta"}
	classification_reference_database_tax = {"16sv4":"trainset16_022016.pds.tax", "16sv34":"trainset16_022016.pds.tax"}
	# Get the location of the MOTHUR_AMPLICON_REF path variable
	path = os.environ['MOTHUR_AMPLICON_HOME']
	# Create a personalised mothur batch scipt for this pipeline run
	log.info ("Creating personalised mothur script...")
	os.system("cp "+path+"/mothur_classify_contigs.batch "+workdir+"/mothur_classify_contigs_personalised.batch")
	os.system("sed -i 's/processors_placeholder/"+str(threads)+"/g' "+workdir+"/mothur_classify_contigs_personalised.batch")
	os.system("sed -i 's?directory_placeholder?"+workdir+"?g' "+workdir+"/mothur_classify_contigs_personalised.batch")
	os.system("sed -i 's/min_seq_size_placeholder/"+str(min_seq_size[amplicon_name.lower()])+"/g' "+workdir+"/mothur_classify_contigs_personalised.batch")
	os.system("sed -i 's/max_seq_size_placeholder/"+str(max_seq_size[amplicon_name.lower()])+"/g' "+workdir+"/mothur_classify_contigs_personalised.batch")
	os.system("sed -i 's|alignment_database_placeholder|"+path+"/"+str(alignment_reference_database[amplicon_name.lower()])+"|g' "+workdir+"/mothur_classify_contigs_personalised.batch")
	os.system("sed -i 's|classification_database_fasta_placeholder|"+path+"/"+str(classification_reference_database_fasta[amplicon_name.lower()])+"|g' "+workdir+"/mothur_classify_contigs_personalised.batch")
	os.system("sed -i 's|classification_database_tax_placeholder|"+path+"/"+str(classification_reference_database_tax[amplicon_name.lower()])+"|g' "+workdir+"/mothur_classify_contigs_personalised.batch")
	# logging parameters and preparing/logging mothur command
	log.info ("Selected parameters for amplicon '%s'\nmin_seq_size:\t%s\nmax_seq_size:\t%s\nalignment reference database:\t%s\nclassification reference database:\t%s\nNo.of threads:\t%s"
			% (amplicon_name, min_seq_size[amplicon_name.lower()],max_seq_size[amplicon_name.lower()], alignment_reference_database[amplicon_name.lower()],
			classification_reference_database_fasta[amplicon_name.lower()], threads))
	command = ("mothur "+workdir+"/mothur_classify_contigs_personalised.batch 2>&1 | tee "+workdir.replace('temp', 'logfiles/mothur_classify_amplicons.log'))
	log.info ("Executing mothur script: 'mothur_classify_contigs_personalised.batch' with command:\n"+command)
	# execute the personalised mothur batch script
	os.system(command)


# Function closes logger handlers
def close_logger(log):
	for handler in log.handlers:
		handler.close()
		log.removeFilter(handler)


# MAIN function
def main():
	# create logger
	logid = random.randint(99999, 9999999)
	log = create_logger(logid)
	# parse command line arguments
	args = parse_arguments(sys.argv, log)
	# creating output directory
	if args.output_dir == 'inputdir':
		output_dir = os.path.abspath(args.input_dir+"/mothur_amplicon_pipeline_output")
	else:
		output_dir = os.path.abspath(args.output_dir)
	log.info("Creating output directory: "+output_dir)
	os.system("mkdir -p "+output_dir)
	# create directories for logfiles, temporary data and output in output_dir
	os.system("mkdir -p "+output_dir+"/temp")
	os.system("mkdir -p "+output_dir+"/logfiles/quality_trimming")
	os.system("mkdir -p "+output_dir+"/logfiles/filtering")
	os.system("mkdir -p "+output_dir+"/quality_control/fastqc")
	os.system("mkdir -p "+output_dir+"/quality_control/multiqc")
	# validate and prepare input data: (data files and samplesheet) for mothur script
	# Generates new samplesheet if none are provided or provided file is invalid
	organise_input_data(args.input_dir, output_dir+"/temp", args.amplicon.lower(), log)
	# trimming, filtering reads and create quality rapports
	log.info("Performing read trimming and read-pair filtering...")
	trimm_and_filter(output_dir+"/temp", args.threads, log)
	# Run mothur script to process amplicon data and classify contigs
	run_mothur_classify(output_dir+"/temp", args.amplicon.lower(), args.threads, log)
	# Move relevant mothur output to output/output_mothur directory
	os.system("cp "+output_dir+"/temp/*pick.tax.summary "+output_dir+"/mothur.classified.tax.summary")
	os.system("cp "+output_dir+"/temp/*pick.taxonomy "+output_dir+"/mothur.classified.taxonomy")
	os.system("cp "+output_dir+"/temp/*pick.pick.count_table "+output_dir+"/mothur.classified.count_table")
	os.system("cp "+output_dir+"/temp/*pick.pick.fasta "+output_dir+"/mothur.classified.fasta")
	# Move batchfile and samplesheet to output/log directory
	os.system("mv "+output_dir+"/temp/mothur_classify_contigs_personalised.batch "+output_dir+"/logfiles/")
	os.system("mv "+output_dir+"/temp/mothur.samplesheet "+output_dir+"/logfiles/")
	# Remove the temp directory if temp parameter is not set to anything but false
	if args.save_temp == "false":
                log.info("save_temp parameter is false (default). Removing temporary files and directories...")
		os.system("rm -rf "+output_dir+"/temp")
	else:
		log.info("save_temp parameter is not false (set in input parameter). Keeping temporary files and directories...")
	log.info("Finishing amplicon_pipeline.py...")
	# close logger handlers
	close_logger(log)
	# move logfiles to output/logs directory
	os.system("mv "+str(logid)+"_mothur_amplicon_pipeline.log "+output_dir+"/logfiles/mothur_amplicon_pipeline.log")
	os.system("mv current_files.summary "+output_dir+"/logfiles/")
	os.system("mv mothur.*.logfile "+output_dir+"/logfiles/")


main()
