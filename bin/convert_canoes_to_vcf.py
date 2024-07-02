#!/usr/bin/python

# Import cvs library
import csv

# Define CANOES output file 
canoes_output_file = 'Sample_CNVs_filtered.csv'

# Output VCF file
vcf_output_file = 'canoes-cnv_calls.vcf'

# Function to convert CANOES output to VCF format
def convert_canoes_to_vcf(canoes_file, vcf_file):
    with open(canoes_file, 'r') as csvfile, open(vcf_file, 'w') as vcf:

        # Write a VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write('##INFO=<ID=SVTYPE, Number=1, Type=String, Description="Type of structural variant">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Create a DictReader object
        reader = csv.DictReader(csvfile)
        for row in reader:
            chrom = row['Chromosome']
            start = row['Start']
            end = row['End']
            cn = row['Copy_Number']

            # Define VCF fields
            ID = '.' # You can set a unique identifier here if available
            REF = '<CNV>'
            ALT = '<CNV>' # Adjust according to your CNV representation
            QUAL = '.' # Quality score, typically not applicable directly
            FILTER = '.' # Filter status, can be set based on CNV quality or other criteria
            INFO = 'SVTYPE=CNV;CN=' + cn # Custom INFO field adjust as needed

            # Write VCF line
            vcf.write(f"{chrom}\t{start}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\n")

            # Convert CANOES output to VCF format
            convert_canoes_to_vcf(canoes_output_file, vcf_output_file)

            print(f"VCF file generated: {vcf_output_file}")
