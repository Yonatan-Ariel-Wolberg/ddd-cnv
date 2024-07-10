#!/usr/bin/python

# CODE THAT CONVERTS CLAMMS OUTPUT BED FILE TO VCF

# Imports os module for file operations
import os

# Imports sys module for system functions
import sys

# Imports re module for regular expressions
import re

# Imports argparse module for argument parsing
import argparse

# Defines function call 'convert_clamms_to_vcf'
def convert_clamms_to_vcf(bed_file, vcf_file):
    
    # Takes BED file as input for reading and VCF file as output for writing
    with open(bed_file, 'r') as infile, open(vcf_file, 'w') as outfile:
        
        # Writes header 
        # Specifies file format as VCF format
        print("Writing VCF header")
        outfile.write("##fileformat=VCFv4.1\n")

        # Specifies source as 'CustomBedToVcfConverter'
        outfile.write("##source=CustomBedToVcfConverter\n")

        # Describes a deletion event
        outfile.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")

        # Describes a duplication event
        outfile.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")

        # Describes an inversion event
        outfile.write("##ALT=<ID=INV,Description=\"Inversion\">\n")

        # Describes a translocation event
        outfile.write("##ALT=<ID=BND,Description=\"Translocation\">\n")

        # Describes an insertion event
        outfile.write("##ALT=<ID=INS,Description=\"Insertion\">\n")

        # Provides the chromosome for the END coordinate in case of a translocation
        outfile.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")

        # Provides the end position for the structural variant
        outfile.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the the structural variant\">\n")

        # Describes an imprecise structural variation
        outfile.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description\"Imprecise structural variation\">\n")

        # Describes a precise structural variation
        outfile.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n")

        # Provides the length of the structural variant
        outfile.write("##INFO=<ID=SVLEN,Number=1,Type=Float,Description=\"Length of the SV\">\n")

        # Provides the method used to identify the SV
        outfile.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Vector of samples supporting the SV\">\n")

        # Provides the type of the SV
        outfile.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")

        # Provides the confidence interval around POS
        outfile.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n")

        # Provides the confidence interval around END
        outfile.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n")

        # Provides the direction of the reads
        outfile.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")

        # Describes the genotype format
        outfile.write("FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

        # Writes CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO fields
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Loops through each line of BED file
        for line in infile:

            # Prints line being processed to console for debugging
            print(f"Processing line: {line}")
            
            # Strips any leading/trailing whitespace from the line and splits it into fields based on tabs
            fields = line.strip().split('\t')
            
            # Assigns each field to a variable
            chrom, start, end, interval, sample, type, cn, num_window, q_some, q_exact, q_left_extend, left_extend_coord, q_right_extend, right_extend_coord, q_left_contract, left_contract_coord, q_right_contract, right_contract_coord = fields
            
            # Converts start and end coordinates to integers
            start = int(start)
            end = int(end)

            # Sets the ID for the VCF entry
            id_ = fields[4]

            # Sets the REF and ALT fields for the VCF entry
            ref = 'N'
            alt = '<' + type + '>'

            # Sets the QUAL field for the VCF entry to Q_SOME
            qual = q_some

            # Sets the FILTER field for the VCF entry to PASS or LowQuality based on Q_EXACT
            filter = 'PASS' if q_exact >= 0 else 'LowQuality'

            # Constructs the INFO field with info about SV type, copy number, sample, number of windows in call, Q_SOME, Q_EXACT, Q_LEFT_EXTEND, LEFT_EXTEND_COORD, Q_RIGHT_EXTEND, RIGHT_EXTEND_COORD, Q_LEFT_CONTRACT, LEFT_CONTRACT_COORD, Q_RIGHT_CONTRACT, RIGHT_CONTRACT_COORD
            info = f'SVTYPE=CNV;CN={cn};SAMPLE={sample};NUM_WINDOW={num_window};Q_SOME={q_some};Q_EXACT={q_exact};Q_LEFT_EXTEND={q_left_extend};LEFT_EXTEND_COORD={left_extend_coord};Q_RIGHT_EXTEND={q_right_extend};RIGHT_EXTEND_COORD={right_extend_coord};Q_LEFT_CONTRACT={q_left_contract};LEFT_CONTRACT_COORD={left_contract_coord};Q_RIGHT_CONTRACT={q_right_contract};RIGHT_CONTRACT_COORD={right_contract_coord}'
            
            # Formats the VCF entry as a tab-separated string
            vcf_line = f'{chrom}\t{start}\t{id_}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n'
            
            # Prints the VCF line being written to the console for debugging
            print(f"Writing line: {vcf_line}")
            
            # Writes formatted VCF line to the output file
            outfile.write(vcf_line)

# Main function execution
# Script only run if executedas the main module
if __name__ == "__main__":
    
    # Creates an argument parser object with a description of the script
    parser = argparse.ArgumentParser(description="Convert CLAMMS output BED file to VCF file.")
    
    # Adds a positional argument for the input BED file
    parser.add_argument("bed_file", help="Input CLAMMS output BED file")
    
    # Adds a positional argument for the output VCF file
    parser.add_argument("vcf_file", help="Output VCF file")
    
    # Parses the command-line arguments and stores them in args.
    args = parser.parse_args()

# Calls the convert_clamms_to_vcf function with the provided BED and VCF file paths
convert_clamms_to_vcf(bed_file, vcf_file)