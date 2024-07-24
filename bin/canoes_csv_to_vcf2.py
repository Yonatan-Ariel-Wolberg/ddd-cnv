#!/usr/bin/env python3

# Library used for parsing command-line arguments
import argparse

# Library used for reading and writing CSV files
import csv

# Library used for checking and creating directories if they don't exist
import os

# Library used for logging errors and other messages
import logging

# Sets up logging configuration
def setup_logging(log_file):
    logging.basicConfig(
        filename=log_file,  # Log messages will be written to this file
        level=logging.ERROR,  # Set the logging level to ERROR
        format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
    )

# Defines a function that converts CANOES output CSV file to VCF files
def convert_canoes_csv_to_vcf(input_file, output_dir, combined_vcf_name, log_file):
    setup_logging(log_file)  # Setup logging

    try:
        # Create the output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Dictionary to store mutations grouped by SAMPLE ID
        mutations_by_sample = {}

        # List to keep track of the order of sample IDs
        sample_order = []

        # Open input file in read mode
        with open(input_file, 'r') as csvfile:
            # Create a CSV reader object
            csvreader = csv.DictReader(csvfile, delimiter='\t')

            # Print CSV headers for debugging
            print("CSV headers:", csvreader.fieldnames)
            
            # List of required columns in the CSV file
            required_columns = ['SAMPLE', 'CNV', 'INTERVAL', 'KB', 'CHR', 'MID_BP', 'TARGETS', 'NUM_TARG', 'MLCN', 'Q_SOME']
            
            # Check if all required columns are present
            missing_columns = [col for col in required_columns if col not in csvreader.fieldnames]
            if missing_columns:
                # If some columns are missing, raise an error and log it
                error_msg = f"Missing required columns: {missing_columns}"
                logging.error(error_msg)
                raise ValueError(error_msg)

            # Read each row from the CSV file
            for row_number, row in enumerate(csvreader, start=1):
                try:
                    # Creates a dictionary called 'mutation' and stores the values from the current row under different keys
                    mutation = {
                        'SAMPLE': row['SAMPLE'],
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

                    # Extracts sample identifier from 'mutation' dictionary and assigns it to variable 'sample_id'
                    sample_id = mutation['SAMPLE']
                    
                    # Line checks whether 'sample_id' is already a key in 'mutations_by_sample'
                    # If not, it creates a new key and assigns an empty list to it
                    if sample_id not in mutations_by_sample:
                        mutations_by_sample[sample_id] = []
                        sample_order.append(sample_id)  # Track the order of sample IDs

                    # Append mutation to the corresponding sample list
                    mutations_by_sample[sample_id].append(mutation)

                # Catches 'KeyError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'KeyError'
                except KeyError as e:
                    error_msg = f"Key error in row {row_number}: {e} - Skipping row: {row}"
                    logging.error(error_msg)
                
                # Catches 'ValueError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'ValueError'
                except ValueError as e:
                    error_msg = f"Value error in row {row_number}: {e} - Skipping row: {row}"
                    logging.error(error_msg)

        # Write all mutations to a single unified VCF file
        if mutations_by_sample:
            unified_vcf_file = os.path.join(output_dir, combined_vcf_name)
            write_vcf_file(unified_vcf_file, mutations_by_sample, sample_order)
        else:
            print("No valid mutations found in the input file.")

        # Write separate VCF files for each sample
        write_vcf_file_per_sample(output_dir, mutations_by_sample, sample_order)

    # Catches 'IOError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'IOError'
    except IOError as e:
        error_msg = f"Error reading or writing file: {e}"
        logging.error(error_msg)
    
    # Catches 'ValueError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'ValueError'
    except ValueError as e:
        error_msg = f"Error: {e}"
        logging.error(error_msg)
   
    # Catches any other unexpected exceptions that occur within corresponding 'try' block, where 'e' is the unexpected issue
    except Exception as e:
        error_msg = f"Unexpected error: {e}"
        logging.error(error_msg)

# Defines a function that writes all mutations to a single VCF file
def write_vcf_file(output_file, mutations_by_sample, sample_order):
    try:
        # Open output file in write mode
        with open(output_file, 'w') as vcf_file:
            # Write VCF header
            vcf_file.write("##fileformat=VCFv4.1\n")
            vcf_file.write("##source=CANOESCSVToVcfConverter\n")
            vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
            vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
            vcf_file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
            vcf_file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n")
            vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Float,Description=\"Length of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Vector of samples supporting the SV\">\n")
            vcf_file.write("##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n")
            vcf_file.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n")
            vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Quality of the sample\">\n")
            
            # Write column headers
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

            # Write sample IDs as column headers in the order they appeared in the CSV file
            for sample_id in sample_order:
                vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # Process mutations for each sample in order
            for sample_id in sample_order:
                # Get all mutations for the current sample
                mutations = mutations_by_sample[sample_id]
                
                # Process each mutation for the current sample
                for mutation in mutations:
                    try:
                        # Extract various fields from 'mutation' dictionary
                        chr_num = mutation['CHR']
                        cnv = mutation['CNV']
                        interval = mutation['INTERVAL']
                        kb = mutation['KB']
                        mid_bp = mutation['MID_BP']
                        targets = mutation['TARGETS']
                        num_targ = mutation['NUM_TARG']
                        cn = mutation['MLCN']
                        q_some = mutation['Q_SOME']
                        
                        # Check if the chr_num is equal to the number before ':' in the interval field (also the chromosome)
                        if chr_num != interval.split(':')[0]:
                            error_msg = "Error: chr_num is not equal to the number before ':' in the INTERVAL field"
                            logging.error(error_msg)
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
                        # CNVs that are over 100 kb in size and have a Q_SOME score over 80 are retained
                        if kb >= 100 and int(q_some) >= 80:
                            filter_status = 'PASS'
            
                        # CNVs that are under 100 kb in size or have a Q_SOME score below 80 are filtered out
                        elif kb < 100 and int(q_some) < 80:
                            filter_status = 'LowQuality'
                        else:
                            filter_status = '.'
                        
                        # Calculate FORMAT values
                        if cn == '3':
                            # If CN = 3, this indicates a heterozygous duplication
                            format_values = '0/1'

                        # Write VCF line for the current mutation
                        vcf_file.write(f"{chr_num}\t{start}\t.\t.\t{cnv}\t.\t{filter_status}\tSVTYPE={cnv};SVLEN={svlen};CIPOS={ci_pos[0]},{ci_pos[1]};CIEND={ci_end[0]},{ci_end[1]};STRANDS={strands}\tGT:Q_SOME\t{format_values}:{q_some}\n")

                    # Catches exceptions that occur within corresponding 'try' block, where 'e' is the unexpected issue
                    except Exception as e:
                        error_msg = f"Error processing mutation {mutation}: {e}"
                        logging.error(error_msg)
    
    # Catches 'IOError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'IOError'
    except IOError as e:
        error_msg = f"Error reading or writing VCF file: {e}"
        logging.error(error_msg)

# Defines a function that writes separate VCF files for each sample
def write_vcf_file_per_sample(output_dir, mutations_by_sample, sample_order):
    for sample_id in sample_order:
        try:
            # Define output file name for the sample
            sample_vcf_file = os.path.join(output_dir, f"{sample_id}.vcf")
            # Write VCF file for the sample
            write_vcf_file(sample_vcf_file, {sample_id: mutations_by_sample[sample_id]}, [sample_id])
        except Exception as e:
            error_msg = f"Error writing VCF file for sample {sample_id}: {e}"
            logging.error(error_msg)

# Main function to handle command-line arguments and run the conversion
def main():
    parser = argparse.ArgumentParser(description="Convert CANOES CSV file to VCF format.")
    parser.add_argument("input_file", help="Path to the CANOES CSV file")
    parser.add_argument("output_dir", help="Directory to save the VCF files")
    parser.add_argument("combined_vcf_name", help="Name for the unified VCF file")
    parser.add_argument("log_file", help="Path to the log file for error reporting")
    args = parser.parse_args()

    convert_canoes_csv_to_vcf(args.input_file, args.output_dir, args.combined_vcf_name, args.log_file)

if __name__ == "__main__":
    main()
