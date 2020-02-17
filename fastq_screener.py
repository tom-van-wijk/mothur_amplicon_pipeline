#!/usr/bin/env python


# Name:         fastq_screener.py
# Author:       Tom van Wijk
# Date:         19-12-2017
# Licence:      GNU General Public License v3.0 (copy provided in directory)

# Script to chacks if sequences in a .fastq file contain any characters other than A/T/C/G.
# Generates a gererates new .fastq file named <input_file)_dada2-filtered without these sequences.
# This script is created to fullfill the requirements of DADA2, however do keep in mind that
# using this scripts on paired end .fastq files migth create orphan sequences.
# Using dada2, mothur or qiime on the data does require you to sort de sequences and remove orphans.


# import python libraries
from argparse import ArgumentParser
import sys


# Function to parse the command-line arguments
# Returns a namespace with argument keys and values
def parse_arguments(args):
        #log.info("Parsing command line arguments...")
	print "Parsing command line arguments..."
        parser = ArgumentParser(prog="fastq_screener.py")
        parser.add_argument("-i", "--infile", dest = "input_file",
                action = "store", default = None, type = str,
                help = "Location of input file (required)",
                required = True)
	return parser.parse_args()


# Function to process fastq_file
# input: path to input .fastq file and output .fastq file
def process_fastq(input_file, output_file):
	line_counter = 1
	with open(input_file, "r") as input_fastq, open(output_file, "w") as output_fastq:
		for line in input_fastq:
			if line.startswith("@M0"):
				line_counter = 1
				header_line = line
				sequence_validation = False
			else:
				line_counter += 1
                        if line_counter == 2:
				# process sequence line further in process_line
                                sequence_validation = process_line(line)
				sequence_line = line
			elif line_counter == 3:
				delimiter_line = line
			elif line_counter == 4:
				quality_line = line
				if sequence_validation == True:
					output_fastq.write(header_line+sequence_line+delimiter_line+quality_line)
	input_fastq.close()
	output_fastq.close()


# Function to process line of fastq_file
# input: sequence line of fastq as string
# output: boolean which is True if line is valid, otherwise False
def process_line(line):
	for char in line.replace("\n",""):
		if char not in "ATCG":
			print "Invalid character:\t'%s'\tfound in line:\n%s" % (char, line.replace("\n",""))
			return False
	return True


# MAIN function
def main():
	# Parse command line arguments
	args = parse_arguments(sys.argv)
	# Process input_file
	output_file = (args.input_file).replace(".fastq","_dada2_filtered.fastq")
	process_fastq(args.input_file, output_file)


main()
