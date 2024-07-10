import argparse
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
            vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
            vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the the structural variant\">\n")
            vcf_file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
            vcf_file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n")
            vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Float,Description=\"Length of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Vector of samples supporting the SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n")
            vcf_file.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n")
            vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            
            # Read each row from the CSV file
            for row in csvreader:
                # Extract fields from the row
                sample = row['SAMPLE']
                cnv = row['CNV']
                interval = row['INTERVAL']
                kb = row['KB']
                chr = row['CHR']
                mid_bp = row['MID_BP']
                targets = row['TARGETS']
                num_targ = row['NUM_TARG']
                cn = row['MLCN']
                q_some = row['Q_SOME']
                
                # Check if the chr is equal to the number before ':' in the interval field
                if chr != interval.split(':')[0]:
                    print("Error: chr is not equal to the number before ':' in the INTERVAL field")
                    continue
                
                # Calculate start and end positions
                start, end = map(int, interval.split(':')[1].split('-'))
                
                # Calculate length of the SV
                svlen = end - start
                
                # Calculate confidence interval around POS and END
                ci_pos = (0, 0)
                ci_end = (0, 0)
                
                # Calculate STRANDS
                strands = '.'
                
                # Calculate FILTER status
                if int(kb) >= 100 and int(q_some) >= 80:
                    filter_status = 'PASS'
                elif int(kb) < 100 and int(q_some) < 80:
                    filter_status = 'LowQuality'
                else:
                    filter_status = '.'
                    
                # Calculate FORMAT values
                if cn == '3':
                    format_values = '0/1'
                elif cn == '1':
                    format_values = '0/1'
                elif cn == '2':
                    format_values = '0/0'
                else:
                    format_values = '.'
                    
                # Calculate QUAL
                qual = '.'
                
                # Write the VCF entry
                vcf_file.write(f"{chr}\t{start}\t.\tN\t{cnv}\t.\t{filter_status};{q_some_filter}\tEND={end};SVLEN={svlen};SVMETHOD=CANOES;SVTYPE={cnv};CIPOS={ci_pos};CIEND={ci_end};STRANDS={strands}\tGT:Q_SOME\t{format_values}:{q_some if q_some else 0}\n")

# Main function execution
if __name__ == "__main__":
    # Create an argument parser object with a description of the script
    parser = argparse.ArgumentParser(description="Convert CANOES CSV file to VCF file.")
    
    # Add a positional argument for the input CSV file
    parser.add_argument("input_file", help="Input CANOES CSV file")
    
    # Add a positional argument for the output VCF file
    parser.add_argument("output_file", help="Output VCF file")
    
    # Parse the command-line arguments and store them in args.
    args = parser.parse_args()
    
    # Call the function to convert the CSV file to VCF file
    convert_canoes_csv_to_vcf(args.input_file, args.output_file)
