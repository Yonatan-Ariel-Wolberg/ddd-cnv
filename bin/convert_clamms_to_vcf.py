#!/usr/bin/python

# clamms_to_vcf.py

import os
import sys

def convert_clamms_to_vcf(cnv_file, output_dir):
    def parse_clamms_cnv(cnv_file):
        sample_cnvs = {}
        with open(cnv_file, 'r') as infile:
            for line in infile:
                fields = line.strip().split('\t')
                sample = fields[0]
                chrom = fields[1]
                start = fields[2]
                end = fields[3]
                cn = fields[4]  # Adjust as per your CNV file format
                qual = fields[5]  # Quality score
                
                # Define VCF fields
                ID = '.'  # You can set a unique identifier here if available
                REF = '<CNV>'
                ALT = '<CNV>'  # Adjust according to your CNV representation
                QUAL = qual  # Quality score from the file
                FILTER = '.'  # Filter status, can be set based on CNV quality or other criteria
                INFO = f'SVTYPE=CNV;CN={cn};SAMPLE={sample}'  # Custom INFO field, adjust as needed
                
                # Create VCF-formatted line
                vcf_line = f"{chrom}\t{start}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\n"
                
                if sample not in sample_cnvs:
                    sample_cnvs[sample] = []
                sample_cnvs[sample].append(vcf_line)
        
        return sample_cnvs
    
    # Parse CLAMMS CNVs
    sample_cnvs = parse_clamms_cnv(cnv_file)
    
    # Write to individual VCF files
    for sample, cnvs in sample_cnvs.items():
        vcf_file = os.path.join(output_dir, f"{sample}_clamms.vcf")
        with open(vcf_file, 'w') as outfile:
            # Write VCF header
            outfile.write("##fileformat=VCFv4.2\n")
            outfile.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            outfile.write('##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number">\n')
            outfile.write('##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample ID">\n')
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Write each CNV line to VCF file
            for vcf_line in cnvs:
                outfile.write(vcf_line)

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python clamms_to_vcf.py <cnv_file> <output_dir>")
        sys.exit(1)

    cnv_file = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    convert_clamms_to_vcf(cnv_file, output_dir)