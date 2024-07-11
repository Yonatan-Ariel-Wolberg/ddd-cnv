<<<<<<<<<<<<<<  ✨ Codeium Command ⭐ >>>>>>>>>>>>>>>>
import csv

def convert_canoes_csv_to_vcf(input_file, output_file):
    # Open input file in read mode
    with open(input_file, 'r') as csvfile:
        # Create a CSV reader object
        csvreader = csv.DictReader(csvfile)
        
        # Open output file in write mode
        with open(output_file, 'w') as vcf_file:
            # Write VCF header
            vcf_file.write("##fileformat=VCFv4.1\n")
            vcf_file.write("##source=CustomBedToVcfConverter\n")
            vcf_file.write("##ALT=<ID=CNV,Description=\"Copy number variation\">\n")
            vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the the structural variant\">\n")
            vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Float,Description=\"Length of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method used to identify the SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
            
            # Get the sample IDs from the CSV file
            sample_ids = csvreader.fieldnames[1:]
            
            # Write the sample-specific fields after the FORMAT field
            for sample_id in sample_ids:
                vcf_file.write(f"\t{sample_id}")
            
            vcf_file.write("\n")
            
            # Read each row from the CSV file
            for row in csvreader:
                # Extract fields from the row
                sample_id = row['SAMPLE']
                cnv = row['CNV']
                interval = row['INTERVAL']
                chr = row['CHR']
                mid_bp = row['MID_BP']
                pos = int(interval.split(':')[1].split('-')[0])
                end = int(interval.split('-')[1])
                cn = int(row['MLCN'])
                
                # Calculate SVLEN
                svlen = pos - end
                
                # Write the VCF record
                vcf_file.write(f"{chr}\t{pos}\t.\tN\t<CNV>{cn}</CNV>\t.\tEND={end};SVLEN={svlen};SVTYPE=CNV;SVMETHOD=CANOES\tGT")
                
                # Write the sample-specific genotype calls
                for sample_id in sample_ids:
                    vcf_file.write(f"\t{sample_id}:{cn}")
                
                vcf_file.write("\n")

# Example usage
input_file = 'input.csv'
output_file = 'output.vcf'
convert_canoes_csv_to_vcf(input_file, output_file)
<<<<<<<  c8175c51-08e3-4847-a7e8-99d37fc36deb  >>>>>>>