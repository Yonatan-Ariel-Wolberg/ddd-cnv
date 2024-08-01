#!/usr/bin/python3

# Import argparse for parsing command-line arguments
import argparse

# Import os for file and directory operations
import os

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

def safe_float(value):
    """
    Function to safely convert a string to a float.
    """
    try:
        return float(value)
    # If the value cannot be converted to a float, return None
    except ValueError:
        return float('nan')
    
def safe_int(value):
    """
    Function to safely convert a string to an integer.
    """
    try:
        return int(value)
    # If the value cannot be converted to an integer, return None
    except ValueError:
        return None

def convert_clamms_bed_to_vcf(input_file, output_dir, combined_vcf_name, log_file, sample_file=None):
    """
    Function to convert CLAMMS BED output to VCF format.

    Parameters:
    - input_file: The path to the CLAMMS BED output file.
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
        # Open the input BED file
        with open(input_file, 'r') as bedfile:
            # Iterate through each line of the BED file
            for line_number, line in enumerate(bedfile):
                try:
                    # Split the line into fields based on the tab delimiter
                    fields = line.strip().split('\t')
                    
                    # Extracts the sample ID from the fourth field
                    sample_id = fields[4]
                    
                    # Only process if sample_id is in sample_list
                    if sample_id not in sample_list:
                        continue
                
                    if len(fields) == 18:
                        # Create a dictionary for the mutation
                        mutation = {
                            'CHROM': fields[0],
                            'START': safe_int(fields[1]),
                            'END': safe_int(fields[2]),
                            'INTERVAL': fields[3],
                            'SAMPLE_ID': sample_id,
                            'CNV_TYPE': fields[5],
                            'MLCN': safe_int(fields[6]),
                            'NUM_WINDOWS': (fields[7]),
                            'Q_SOME': safe_float(fields[8]),
                            'Q_EXACT': safe_float(fields[9]),
                            'Q_LEFT_EXTEND': safe_float(fields[10]),
                            'LEFT_EXTEND_COORD': safe_int(fields[11]),
                            'Q_RIGHT_EXTEND': safe_float(fields[12]),
                            'RIGHT_EXTEND_COORD': safe_int(fields[13]),
                            'Q_LEFT_CONTRACT': safe_float(fields[14]),
                            'LEFT_CONTRACT_COORD': safe_int(fields[15]),
                            'Q_RIGHT_CONTRACT': safe_float(fields[16]),
                            'RIGHT_CONTRACT_COORD': safe_int(fields[17])
                        }
                
                    else:
                        logging.error(f"Invalid number of fields in line {line_number}: {line}")
                        continue

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

                # Catches 'KeyError' exceptions within the corresponding 'try' block, where 'e' is the specific 'KeyError' object that was caught and 'line_number' is the line number of the exception that was caught
                except KeyError as e:
                    error_msg = f"Error processing line {line_number}: {e} - Skipping line: {line}"
                    logging.error(error_msg)
                
                # Catches 'ValueError' exceptions within the corresponding 'try' block, where 'e' is the specific 'ValueError' object that was caught and 'line_number' is the line number of the exception that was caught
                except ValueError as e:
                    error_msg = f"Error processing line {line_number}: {e} - Skipping line: {line}"
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

    # Add placeholder mutations for samples that do not have any mutations
    for sample_id in sample_list:
        if sample_id not in mutations_by_sample:
            mutations_by_sample[sample_id] = []
            sample_order.append(sample_id)

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

    # Set combined VCF file name to 'XHMM_CNVs_filtered.vcf' if not provided
    if not combined_vcf_name:
        combined_vcf_file = 'CLAMMS_CNVs_filtered.vcf'
    
    # Write individual VCF files for each sample
    for sample_id in sample_list:
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
            vcf_file.write("##source=CLAMMS-BedToVcfConverter\n")
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
            vcf_file.write("##INFO=<ID=MLCN,Number=1,Type=Integer,Description=\"Most Likely Copy number of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")
            vcf_file.write("##INFO=<ID=NUM_WINDOWS,Number=1,Type=String,Description=\"Number of windows used to estimate copy number\">\n")
            vcf_file.write("##INFO=<ID=Q_LEFT_EXTEND,Number=1,Type=Float,Description=\"Phred-scaled quality of the left breakpoint\">\n")
            vcf_file.write("##INFO=<ID=Q_RIGHT_EXTEND,Number=1,Type=Float,Description=\"Phred-scaled quality of the right breakpoint\">\n")
            vcf_file.write("##INFO=<ID=LEFT_EXTEND_COORD,Number=1,Type=Integer,Description=\"Add this to the CNV start coordinate to get the start coordinate of the first window to the left of the called CNV\">\n")   
            vcf_file.write("##INFO=<ID=RIGHT_EXTEND_COORD,Number=1,Type=Integer,Description=\"Add this to the CNV end coordinate to get the end coordinate of the first window to the right of the called CNV\">\n")
            vcf_file.write("##INFO=<ID=Q_LEFT_CONTRACT,Number=1,Type=Float,Description=\"Phred-scaled quality of the left breakpoint\">\n")
            vcf_file.write("##INFO=<ID=Q_RIGHT_CONTRACT,Number=1,Type=Float,Description=\"Phred-scaled quality of the right breakpoint\">\n")
            vcf_file.write("##INFO=<ID=LEFT_CONTRACT_COORD,Number=1,Type=Integer,Description=\"Add this to the CNV start coordinate to get the start coordinate of the second window of the called CNV\">\n")
            vcf_file.write("##INFO=<ID=RIGHT_CONTRACT_COORD,Number=1,Type=Integer,Description=\"Add this to the CNV end coordinate to get the end coordinate of the second-to-last window of the called CNV\">\n")         
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Phred-scaled quality of some CNV event in the interval\">\n")
            vcf_file.write("##FORMAT=<ID=Q_EXACT,Number=1,Type=Float,Description=\"Phred-scaled quality of the exact CNV event along the entire interval\">\n")  
            
            # Write VCF header row with sample IDs
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            
            for sample_id in sample_order:
                # Write sample IDs as header columns
                vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # For each sample in sample_order, write mutations
            for sample_id in sample_order:
                # Get list of mutations for current sample
                mutations = mutations_by_sample[sample_id]
                if not mutations and include_placeholders:
                    placeholder_line = f".\t.\t.\t.\t.\t.\t.\t.\t.\t./.:.:."
                    vcf_file.write(placeholder_line + "\n")
                else:
                    for mutation in mutations:
                        try:
                            # Extract relevant information from each field
                            chr_num = mutation['CHROM']
                            cnv = mutation['CNV_TYPE']
                            interval = mutation['INTERVAL']
                            start = int(mutation['START'])
                            end = int(mutation['END'])
                            mlcn = int(mutation['MLCN'])
                            num_windows = int(mutation['NUM_WINDOWS'])
                            q_some = float(mutation['Q_SOME'])
                            q_exact = float(mutation['Q_EXACT'])
                            q_left_extend = float(mutation['Q_LEFT_EXTEND'])
                            left_extend_coord = int(mutation['LEFT_EXTEND_COORD'])
                            q_right_extend = float(mutation['Q_RIGHT_EXTEND'])
                            right_extend_coord = int(mutation['RIGHT_EXTEND_COORD'])
                            q_left_contract = float(mutation['Q_LEFT_CONTRACT'])
                            left_contract_coord = int(mutation['LEFT_CONTRACT_COORD'])
                            q_right_contract = float(mutation['Q_RIGHT_CONTRACT'])
                            right_contract_coord = int(mutation['RIGHT_CONTRACT_COORD'])
                            
                        
                            # Set reference allele to 'N'
                            ref = 'N'
                        
                            # Set alternate allele based on CNV type
                            alt = '<' + cnv + '>'
                        
                            # Checks if the chromosome number in the interval matches the specified chromosome number
                            if chr_num != interval.split(':')[0]:
                                # Raises an error if the chromosome number in the interval does not match the specified chromosome number
                                raise ValueError(f"Chromosome number in interval {interval} does not match the specified chromosome number {chr_num}.")
                        
                            # SVLEN is the length of the structural variant
                            svlen = end - start + 1                      
                        
                            # Set filter status based on quality scores
                            if q_some >= 500 and q_exact > 0:
                                filter_status = 'PASS'
                            else:
                                filter_status = 'LowQuality'
                        
                            # Placeholder for strand information as XHMM does not provide strand information
                            strands = "."
                        
                            # Placeholder for genotype as XHMM does not provide genotype
                            gt = "./."

                            # Write the mutation to the VCF file
                            vcf_file.write(f"{chr_num}\t{start}\t.\t{ref}\t{alt}\t.\t{filter_status}\tEND={end};SVLEN={svlen};SVTYPE={cnv};SVMETHOD=CLAMMS;STRANDS={strands};NUM_WINDOWS={num_windows};MLCN={mlcn};Q_LEFT_EXTEND={q_left_extend};Q_RIGHT_EXTEND={q_right_extend};LEFT_EXTEND_COORD={left_extend_coord};RIGHT_EXTEND_COORD={right_extend_coord};Q_LEFT_CONTRACT={q_left_contract};LEFT_CONTRACT_COORD={left_contract_coord};Q_RIGHT_CONTRACT={q_right_contract};RIGHT_CONTRACT_COORD={right_contract_coord}\tGT:Q_EXACT:Q_SOME:Q_EXACT")

                            # Write 'FORMAT' values for each sample
                            for sample in sample_order:
                                if sample == sample_id:
                                    vcf_file.write(f"\t{gt}:{q_some if q_some else '.'}:{q_exact if q_exact else '.'}")
                                else:
                                    vcf_file.write(f"\t./.:.:.")
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
    Main function of the script. It parses the command-line arguments and calls the convert_clamms_bed_to_vcf function.
    """
    # Create an argument parser with a description of the script
    parser = argparse.ArgumentParser(description='Convert CLAMMS .bed files to VCF format')
    
    # Add arguments to the parser for input file, output directory, combined VCF file name and log file
    parser.add_argument('-i', '--input', type=str, required=True, help='Input .bed file')
    parser.add_argument('-d', '--output_dir', type=str, help='Output directory for VCF files')
    parser.add_argument('-c', '--combined_vcf_name', type=str, help='Name for the combined VCF file')
    parser.add_argument('-l', '--log', type=str, help='Log file to record errors and issues')
    parser.add_argument('-s', '--sample_file', type=str, help='File containing sample IDs')

    # Parse command-line arguments and store them in args
    args = parser.parse_args()

    # Call the convert_xhmm_xcnv_to_vcf function with the parsed arguments
    convert_clamms_bed_to_vcf(args.input, args.output_dir, args.combined_vcf_name, args.log, args.sample_file)

# Main function execution
# Script will only run if executed and not run as a module
if __name__ == "__main__":
    main()