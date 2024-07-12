#!/usr/bin/python

import argparse
import csv
import os

def convert_canoes_csv_to_vcf(input_file, output_file):
    
    # Create the directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # List to store all mutations
    all_mutations = []

    # Open input file in read mode
    with open(input_file, 'r') as csvfile:
        # Create a CSV reader object
        csvreader = csv.DictReader(csvfile, delimiter='\t')

        # Print CSV headers for debugging
        print("CSV headers:", csvreader.fieldnames)
        
        # Check if all required columns are present
        required_columns = ['SAMPLE', 'CNV', 'INTERVAL', 'KB', 'CHR', 'MID_BP', 'TARGETS', 'NUM_TARG', 'MLCN', 'Q_SOME']
        missing_columns = [col for col in required_columns if col not in csvreader.fieldnames]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")

        # Read each row from the CSV file
        for row in csvreader:
            try:
                sample_id = row['SAMPLE']
                mutation = {
                    'CNV': row['CNV'],
                    'INTERVAL': row['INTERVAL'],
                    'KB': float(row['KB']),
                    'CHR': row['CHR'],
                    'MID_BP': row['MID_BP'],
                    'TARGETS': row['TARGETS'],
                    'NUM_TARG': row['NUM_TARG'],
                    'MLCN': row['MLCN'],
                    'Q_SOME': row['Q_SOME']
                }

                all_mutations.append(mutation)

            except KeyError as e:
                print(f"Key error: {e}")
            except ValueError as e:
                print(f"Value error: {e}")

    # Write all mutations to a single VCF file
    write_vcf_file(output_file, all_mutations)

def write_vcf_file(output_file, mutations):
    # Open output file in write mode
    with open(output_file, 'w') as vcf_file:
        # Write VCF header
        vcf_file.write("##fileformat=VCFv4.1\n")
        vcf_file.write("##source=CustomBedToVcfConverter\n")
        vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
        vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
        vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
        vcf_file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
        vcf_file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n")
        vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Float,Description=\"Length of the SV\">\n")
        vcf_file.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Vector of samples supporting the SV\">\n")
        vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
        vcf_file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n")
        vcf_file.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n")
        vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")
        vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        
        # Write the sample-specific fields after the FORMAT field
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        
        # Write sample IDs as column headers
        sample_ids = set(m['SAMPLE'] for m in mutations)
        for sample_id in sample_ids:
            vcf_file.write(f"\t{sample_id}")
        vcf_file.write("\n")

        # Process each mutation
        for mutation in mutations:
            try:
                chr_num = mutation['CHR']
                cnv = mutation['CNV']
                interval = mutation['INTERVAL']
                kb = mutation['KB']
                mid_bp = mutation['MID_BP']
                targets = mutation['TARGETS']
                num_targ = mutation['NUM_TARG']
                cn = mutation['MLCN']
                q_some = mutation['Q_SOME']
            
                # Check if the chr_num is equal to the number before ':' in the interval field
                if chr_num != interval.split(':')[0]:
                    print("Error: chr_num is not equal to the number before ':' in the INTERVAL field")
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
                if kb >= 100 and int(q_some) >= 80:
                    filter_status = 'PASS'
                elif kb < 100 and int(q_some) < 80:
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
                
                # Write the VCF entry with 'chr' prefix
                vcf_file.write(f"chr{chr_num}\t{start}\t.\tN\t{cnv}\t.\t{filter_status}\tEND={end};SVLEN={svlen};SVMETHOD=CANOES;SVTYPE=CNV;CIPOS={ci_pos};CIEND={ci_end};STRANDS={strands}\tGT:Q_SOME")
                
                # Write FORMAT values for each sample
                for sample_id in sample_ids:
                    if sample_id == mutation['SAMPLE']:
                        vcf_file.write(f"\t{format_values}:{q_some if q_some else 0}")
                    else:
                        vcf_file.write("\t.")

                vcf_file.write("\n")

            except KeyError as e:
                print(f"Key error: {e}")
            except ValueError as e:
                print(f"Value error: {e}")

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
