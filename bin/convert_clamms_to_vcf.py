import os

def convert_clamms_to_vcf(cnv_file, output_dir):
    def parse_clamms_cnv(cnv_file):
        cnvs = []
        with open(cnv_file, 'r') as infile:
            for line in infile:
                fields = line.strip().split('\t')
                chrom = fields[0]
                start = fields[1]
                end = fields[2]
                cn = fields[3]  # Adjust as per your CNV file format
                
                # Define VCF fields
                ID = '.'  # You can set a unique identifier here if available
                REF = '<CNV>'
                ALT = '<CNV>'  # Adjust according to your CNV representation
                QUAL = '.'  # Quality score, typically not applicable directly
                FILTER = '.'  # Filter status, can be set based on CNV quality or other criteria
                INFO = f'SVTYPE=CNV;CN={cn}'  # Custom INFO field, adjust as needed
                
                # Append VCF-formatted line to list
                vcf_line = f"{chrom}\t{start}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\n"
                cnvs.append(vcf_line)
        
        return cnvs
    
    # Parse CLAMMS CNVs
    sample_name = os.path.basename(cnv_file).split('.')[0]
    cnvs = parse_clamms_cnv(cnv_file)
    
    # Write to VCF file
    vcf_file = os.path.join(output_dir, f"{sample_name}_clamms.vcf")
    with open(vcf_file, 'w') as outfile:
        # Write VCF header
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Write each CNV line to VCF file
        for vcf_line in cnvs:
            outfile.write(vcf_line)

# Example usage for standalone testing
if __name__ == "__main__":
    import sys
    cnv_file = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    convert_clamms_to_vcf(cnv_file, output_dir)