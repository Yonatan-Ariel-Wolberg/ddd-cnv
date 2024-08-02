#!/usr/bin/python3

# Import argparse for parsing command-line arguments
import argparse

# Import os for file and directory operations
import os

# Import csv for reading and writing CSV files
import csv

# Import logging for logging errors and other messages
import logging

# Import datetime for getting the current date and time
from datetime import datetime

def setup_logging(log_file=None):
    """
    Set up the logging configuration.
    
    Parameters:
    - log_file: Optional log file to record errors.
    """
    if log_file:    # Check if a log file path is provided
        logging.basicConfig(
            filename=log_file,  # Log messages will be written to this file
            level=logging.ERROR,    # Set the logging level to ERROR
            format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
        )
    else:    # If no log file path is provided
        logging.basicConfig(
            level=logging.ERROR,    # Set the logging level to ERROR
            format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
        )

def read_sample_list(sample_file):
    """
    Function to read sample IDs from a .txt file.

    Parameters:
    - sample_file: The path to the .txt file containing sample IDs, one per line.

    Returns:
    - A list of sample IDs.
    """
    sample_set = set()
    try:
        with open(sample_file, 'r') as f:
            for line in f:
                # Read sample IDs, strip whitespace and convert to a set
                sample_id = line.strip()
                if sample_id:
                    sample_set.add(sample_id)

    # Catches any 'IOError' errors while reading sample IDs, where 'e' is the error object
    except IOError as e:
        logging.error(f"Error reading sample IDs from {sample_file}: {e}")
    return list(sample_set)


def convert_canoes_csv_to_vcf(input_file, output_dir, combined_vcf_name, log_file, sample_file=None):
    """
    Function to convert CANOES CSV output to VCF format.

    Parameters:
    - input_file: The path to the CANOES CSV output file.
    - output_dir: The path to the output directory.
    - combined_vcf_name: The name of the combined VCF file.
    - sample_file: Optional .txt file containing sample IDs, one per line.
    - log_file: Optional log file to record errors.
    """
    
    # Set up logging if a log file is specified
    if log_file:
        setup_logging(log_file)

    # Create the output directory if it does not exist
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        
        # Handle any errors that may occur while creating the directory
        except Exception as e:
            logging.error(f"Error creating directory {output_dir}: {e}")
            return

    # Read sample IDs from file if provided
    sample_list = read_sample_list(sample_file) if sample_file else []

    # Dictionary to store mutations grouped by SAMPLE ID
    mutations_by_sample = {}
    
    # List to keep track of the order of sample IDs
    sample_order = []

    try:
        # Open the input CSV file
        with open(input_file, 'r') as csvfile:
            # Create a CSV reader object
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            
            # Print CSV headers for debugging
            print("CSV headers:", csvreader.fieldnames)

            # List of required columns in CSV file
            required_columns = [
                'SAMPLE', 'CNV', 'INTERVAL', 'KB', 'CHR', 'MID_BP', 'TARGETS',
                'NUM_TARG', 'MLCN', 'Q_SOME'
            ]

            # Check if all required columns are present
            missing_columns = [col for col in required_columns if col not in csvreader.fieldnames]
            if missing_columns:
                # If any required columns are missing, raise an error and log it
                error_msg = f"Missing required columns: {missing_columns}"
                logging.error(error_msg)
                # Raise a ValueError to stop the script
                raise ValueError(error_msg)

            for row_number, row in enumerate(csvreader, start=1):
                try:
                    
                    # Extracts sample identifier from the row
                    sample_id = row['SAMPLE']

                    # Only process if sample_id is in sample_list
                    if sample_id not in sample_list:
                        continue

                    # Creates a dictionary called 'mutation' and stores values from the current row under different keys
                    mutation = {
                        'SAMPLE': sample_id,
                        'CNV': row['CNV'],
                        'INTERVAL': row['INTERVAL'],
                        'KB': row['KB'],
                        'CHR': row['CHR'],
                        'MID_BP': row['MID_BP'],
                        'TARGETS': row['TARGETS'],
                        'NUM_TARG': row['NUM_TARG'],
                        'MLCN': row['MLCN'],
                        'Q_SOME': row['Q_SOME']
                    }

                    # Appends 'sample_id' to 'sample_order' list
                    if sample_id not in sample_order:
                        sample_order.append(sample_id)

                    # Line checks if 'sample_id' is already in the 'mutations_by_sample' dictionary
                    if sample_id in mutations_by_sample:
                    # Line appends 'mutation' dictionary to the 'mutations_by_sample' dictionary at the key 'sample_id'    
                        mutations_by_sample[sample_id].append(mutation)
                    else:
                        # If 'sample_id' is not in the 'mutations_by_sample' dictionary, it creates a new key-value pair with the 'mutation' dictionary as the value
                        mutations_by_sample[sample_id] = [mutation]

                # Catches 'KeyError' exceptions within the corresponding 'try' block, where 'e' is the specific 'KeyError' object that was caught and 'row_number' is the line number of the exception that was caught
                except KeyError as e:
                    error_msg = f"Error processing row {row_number}: {e} - Skipping row: {row}"
                    logging.error(error_msg)
                
                # Catches 'ValueError' exceptions within the corresponding 'try' block, where 'e' is the specific 'ValueError' object that was caught and 'row_number' is the line number of the exception that was caught
                except ValueError as e:
                    error_msg = f"Error processing row {row_number}: {e} - Skipping row: {row}"
                    logging.error(error_msg)

    # Catches 'IOError' exceptions within the corresponding 'try' block, where 'e' is the specific 'IOError' object that was caught
    except IOError as e:
        error_msg = f"Error reading {input_file}: {e}"
        logging.error(error_msg)
        return
    
    # Catches 'ValueError' exceptions within the corresponding 'try' block, where 'e' is the specific 'ValueError' object that was caught
    except ValueError as e:
        error_msg = f"Value error in {input_file}: {e}"
        logging.error(error_msg)
        return
    
    # Catches unexpected errors within the corresponding 'try' block where 'e' is the specific 'Exception' object that was caught
    except Exception as e:
        error_msg = f"Unexpected error in {input_file}: {e}"
        logging.error(error_msg)
        return

    sample_order = [sample_id for sample_id in sample_order]

    # Check for duplicates in sample_order
    if len(sample_order) != len(set(sample_order)):
        logging.error("Duplicates found in sample_order!")
        print("Duplicates found in sample_order!")
        return

    # Write combined VCF file
    if combined_vcf_name:
        # Create combined VCF file name
        combined_vcf_file = os.path.join(output_dir, combined_vcf_name)
        write_vcf_file(combined_vcf_file, mutations_by_sample, sample_list, log_file, include_placeholders=False)

    # Write individual VCF files for each sample
    for sample_id in sample_list:
        if sample_id not in mutations_by_sample and len(mutations_by_sample) > 0:
            # Names each individual VCF file after the sample ID
            sample_vcf_file = os.path.join(output_dir, f"{sample_id}.vcf")
            write_vcf_file(sample_vcf_file, {sample_id: mutations_by_sample[sample_id]}, [sample_id], log_file, include_placeholders=True)

    # Print number of samples in combined VCF file or in individual VCF files
    if combined_vcf_name:
        print(f"Number of samples in combined VCF file: {len(sample_order)}")
    else:
        print(f"Number of samples in individual VCF files: {len(sample_order)}")


def write_vcf_file(output_file, mutations_by_sample, sample_order, log_file=None, include_placeholders=False):
    """
    Write mutations to a VCF file.

    Args:
        output_file (str): Output VCF file path.
        mutations_by_sample (dict): Dictionary mapping sample IDs to mutations.
        sample_order (list): List of sample IDs.
        log_file (str, optional): Log file path. Defaults to None.
    """
    
    # Setup logging if log_file is provided
    if log_file:
        setup_logging(log_file)

    try:
        # Open output VCF file for writing
        with open(output_file, 'w') as vcf_file:
            # Write VCF header
            vcf_file.write("##fileformat=VCFv4.1\n")
            vcf_file.write("##fileDate=" + datetime.now().strftime("%d%m%Y") + "\n")
            vcf_file.write("##source=CANOES-CSVToVcfConverter\n")
            vcf_file.write("##reference=GRCh37.p13\n")
            vcf_file.write("##phasing=partial\n")
            for i in range(1, 23):
                vcf_file.write(f"##contig=<ID=chr{i}>\n")
            vcf_file.write("##contig=<ID=chrX>\n")
            vcf_file.write("##contig=<ID=chrY>\n")
            vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
            vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
            vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Float,Description=\"Length of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Vector of samples supporting the SV\">\n")
            vcf_file.write("##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")
            vcf_file.write("##INFO=<ID=Q_START,Number=1,Type=Float,Description=\"Phred-scaled quality of left breakpoint of CNV\">\n") 
            vcf_file.write("##INFO=<ID=Q_STOP,Number=1,Type=Float,Description=\"Phred-scaled quality of right breakpoint of CNV\">\n")              
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Phred-scaled quality of some CNV event in the interval\">\n") 

            # Write VCF header row with sample IDs
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            
            for sample_id in sample_order:
                
                vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # For each sample in sample_order, write mutations
            for sample_id in sample_order:
                # Get list of mutations for current sample
                mutations = mutations_by_sample[sample_id]
                for mutation in mutations:
                    try:
                        # Extract relevant information from each field
                        chr_num = mutation['CHR']
                        cnv = mutation['CNV']
                        interval = mutation['INTERVAL']
                        kb = mutation['KB']
                        mid_bp = mutation['MID_BP']
                        targets = mutation['TARGETS']
                        num_targ = mutation['NUM_TARG']
                        cn = mutation['MLCN']
                        q_some = float(mutation['Q_SOME'])
                        
                        # Set reference allele to 'N'
                        ref = 'N'
                        
                        # Set alternate allele based on CNV type
                        alt = '<' + cnv + '>'
                        
                        # Checks if the chromosome number in the interval matches the specified chromosome number
                        if chr_num != interval.split(':')[0]:
                            # Raises an error if the chromosome number in the interval does not match the specified chromosome number
                            raise ValueError(f"Chromosome number in interval {interval} does not match the specified chromosome number {chr_num}.")
                        
                        # Start is the start of the interval and end is the end of the interval 
                        start = int(interval.split(':')[1].split('-')[0])
                        end = int(interval.split(':')[1].split('-')[1])
                        
                        # SVLEN is the length of the structural variant
                        svlen = end - start + 1                      
                        
                        # Set filter status based on quality scores
                        if q_some >= 80 or q_some != 'NA':
                            filter_status = 'PASS'
                        else:
                            filter_status = 'LowQuality'
                        
                        # Placeholder for strand information as CANOES does not provide strand information
                        strands = "."
                        
                        # Placeholder for genotype as CANOES does not provide genotype
                        gt = "./."

                        # Write the mutation to the VCF file
                        vcf_file.write(f"chr{chr_num}\t{start}\t.\t{ref}\t{alt}\t.\t{filter_status}\tEND={end};SVLEN={svlen};SVTYPE=CNV;CN={cn};SVMETHOD=CANOES;STRANDS={strands}\tGT:Q_SOME")

                        # Write 'FORMAT' values for each sample
                        for sample in sample_order:
                            if sample == sample_id:
                                vcf_file.write(f"\t{gt}:{q_some if q_some else '.'}")
                            else:
                                vcf_file.write(f"\t./.:.")
                        vcf_file.write("\n")
                    
                    # Catches 'KeyError' exceptions within the corresponding 'try' block, where 'e' is specifically the 'KeyError' object
                    except KeyError as e:
                        logging.error(f"Key error when writing mutation {mutation}: {e}")
                    
                    # Catches 'ValueError' exceptions within the corresponding 'try' block, where 'e' is specifically the 'ValueError' object
                    except ValueError as e:
                        logging.error(f"Value error when writing mutation {mutation}: {e}")

    # Catches 'IOError' exceptions within the corresponding 'try' block, where 'e' is specifically the 'IOError' object
    except IOError as e:
        logging.error(f"Error writing VCF to {output_file}: {e}")
    
    # Catches unexpected errors within the corresponding 'try' block, where 'e' is specifically the 'Exception' object
    except Exception as e:
        logging.error(f"Unexpected error writing VCF to {output_file}: {e}")



def main():
    """
    Main function of the script. It parses the command-line arguments and calls the convert_canoes_csv_to_vcf function.
    """
    # Create an argument parser with a description of the script
    parser = argparse.ArgumentParser(description='Convert CANOES .csv files to VCF format')
    
    # Add arguments to the parser for input file, output directory, combined VCF file name and log file
    parser.add_argument('-i', '--input', type=str, required=True, help='Input .csv file')
    parser.add_argument('-d', '--output_dir', type=str, help='Output directory for VCF files')
    parser.add_argument('-c', '--combined_vcf_name', type=str, help='Name for the combined VCF file')
    parser.add_argument('-l', '--log', type=str, help='Log file to record errors and issues')
    parser.add_argument('-s', '--sample_file', type=str, help='File containing sample IDs')

    # Parse command-line arguments and store them in args
    args = parser.parse_args()

    # Call the convert_xhmm_xcnv_to_vcf function with the parsed arguments
    convert_canoes_csv_to_vcf(args.input, args.output_dir, args.combined_vcf_name, args.log, args.sample_file)

# Main function execution
# Script will only run if executed and not run as a module
if __name__ == "__main__":
    main()