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

# Set up logging configuration
def setup_logging(log_file=None):
    if log_file:  # Check if a log file path is provided
        logging.basicConfig(
            filename=log_file,  # Log messages will be written to this file
            level=logging.ERROR,  # Set the logging level to ERROR
            format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
        )
    else:  # If no log file path is provided
        logging.basicConfig(
            level=logging.ERROR,  # Set the logging level to ERROR
            format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
        )

def convert_xhmm_xcnv_to_vcf(input_file, output_dir, combined_vcf_name, log_file):

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
            
            # List of required columns in the CSV file
            required_columns = ['SAMPLE', 'CNV', 'INTERVAL', 'KB', 'CHR', 'MID_BP', 'TARGETS', 'NUM_TARG', 'Q_EXACT', 'Q_SOME', 'Q_NON_DIPLOID', 'Q_START', 'Q_STOP', 'MEAN_RD', 'MEAN_ORIG_RD']
        
            # Check if all required columns are present
            missing_columns = [col for col in required_columns if col not in csvreader.fieldnames]
            if missing_columns:
                # If some columns are missing, raise an error and log it
                error_msg = f"Missing required columns: {missing_columns}"
                logging.error(error_msg)
                # Raise a ValueError to stop the script
                raise ValueError(error_msg)
        
            for row_number, row in enumerate(csvreader, start=1):
                try:
                   # Creates a dictionary called 'mutation' and stores values from the current row under different keys
                    mutation = {
                        'SAMPLE': row['SAMPLE'],
                        'CNV': row['CNV'],
                        'INTERVAL': row['INTERVAL'],
                        'KB': row['KB'],
                        'CHR': row['CHR'],
                        'MID_BP': row['MID_BP'],
                        'TARGETS': row['TARGETS'],
                        'NUM_TARG': row['NUM_TARG'],
                        'Q_EXACT': row['Q_EXACT'],
                        'Q_SOME': row['Q_SOME'],
                        'Q_NON_DIPLOID': row['Q_NON_DIPLOID'],
                        'Q_START': row['Q_START'],
                        'Q_STOP': row['Q_STOP'],
                        'MEAN_RD': row['MEAN_RD'],
                        'MEAN_ORIG_RD': row['MEAN_ORIG_RD']
                    }
                   
                    # Extracts sample identifier from 'mutation' dictionary and assigns it to the variable 'sample_id'
                    sample_id = mutation['SAMPLE']
                   
                    # Line appends 'sample_id' to the 'sample_order' list
                    sample_order.append(sample_id)
                   
                    # Line checks whether 'sample_id' is already a key in the 'mutations_by_sample' dictionary
                    if sample_id in mutations_by_sample:
                       # Line appends 'mutation' to the list associated with 'sample_id' in the 'mutations_by_sample' dictionary
                       mutations_by_sample[sample_id].append(mutation)
                    else:
                       # Line creates a new key-value pair in the 'mutations_by_sample' dictionary
                       mutations_by_sample[sample_id] = [mutation]
                
                # Catches 'KeyError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'KeyError' object        
                except KeyError as e:
                   error_msg = f"Error processing row {row_number}: {e} - Skipping row: {row}"
                   logging.error(error_msg)
                   
                # Catches 'ValueError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'ValueError' object
                except ValueError as e:
                   error_msg = f"Error processing row {row_number}: {e} - Skipping row: {row}"
                   logging.error(error_msg)
        
    # Catches 'IOError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'IOError' object
    except IOError as e:
        error_msg = f"Error reading {input_file}: {e}"
        # Log the error and exit the program
        logging.error(error_msg)
        return
    
    # Catches 'ValueError' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'ValueError' object
    except ValueError as e:
        error_msg = f"Value error in {input_file}: {e}" 
        # Log the error and exit the program   
        logging.error(error_msg) 
        return
    
    # Catches 'Exception' exceptions that occur within corresponding 'try' block, where 'e' is the specific 'Exception' object
    except Exception as e:
        error_msg = f"Unexpected error in {input_file}: {e}"
        # Log the error and exit the program
        logging.error(error_msg)
        return
    
    # Write combined VCF file
    if combined_vcf_name:
        # Create the combined VCF file name
        combined_vcf_file = os.path.join(output_dir, combined_vcf_name)
        write_vcf_file(combined_vcf_file, mutations_by_sample, sample_order, log_file)
        
    
    # Write individual VCF files for each sample
    for sample_id in sample_order:
        sample_vcf_file = os.path.join(output_dir, f"{sample_id}.vcf")
        write_vcf_file(sample_vcf_file, {sample_id: mutations_by_sample[sample_id]}, [sample_id], log_file)
        
    # Print the number of samples in the combined VCF file
    if combined_vcf_name:
        print(f"Number of samples in combined VCF file: {len(sample_order)}")
    else:
        print(f"Number of samples in individual VCF files: {len(sample_order)}")
    return

# Defines a function that writes all mutations to a single VCF file:
def write_vcf_file(output_file, mutations_by_sample, sample_order, log_file=None):
    """
    Write VCF file for a single sample.
    
    Args:
    - output_file (str): Path to the output VCF file.
    - mutations_by_sample (dict): Dictionary with sample IDs as keys and lists of mutations as values.
    - sample_order (list): Ordered list of sample IDs.
    - log_file (str, optional): Path to the log file. Defaults to None.
    """
    # Set up logging if a log file is specified
    if log_file:
        setup_logging(log_file)
    
    try:
        # Open the output VCF file for writing
        with open(output_file, 'w') as vcf_file:
            # Write VCF header
            vcf_file.write("##fileformat=VCFv4.1\n")
            vcf_file.write("##fileDate=" + datetime.now().strftime("%d%m%Y") + "\n")
            vcf_file.write("##source=XHMM-XCNVToVcfConverter\n")
            vcf_file.write("##reference==GRCh37.p13\n")
            vcf_file.write("##phasing=partial\n")
            vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
            vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
            vcf_file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
            vcf_file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n")
            vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Float,Description=\"Length of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Vector of samples supporting the SV\">\n")
            vcf_file.write("##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")
            vcf_file.write("##INFO=<ID=Q_START,Number=1,Type=Float,Description=\"Phred-scaled quality of left breakpoint of CNV\">\n") 
            vcf_file.write("##INFO=<ID=Q_STOP,Number=1,Type=Float,Description=\"Phred-scaled quality of right breakpoint of CNV\">\n")              
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Phred-scaled quality of some CNV event in the interval\">\n")
            vcf_file.write("##FORMAT=<ID=Q_EXACT,Number=1,Type=Float,Description=\"Phred-scaled quality of the exact CNV event along the entire interval\">\n")
            vcf_file.write("##FORMAT=<ID=Q_NON_DIPLOID,Number=1,Type=Float,Description=\"Phred-scaled quality of not being diploid, i.e., DEL or DUP event in the interval\">\n")   
            
            # Write the VCF header row with sample_IDs
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            
            for sample_id in sample_order:
                vcf_file.write("\t" + sample_id)
            vcf_file.write("\n")
            
            # Write mutations for each sample
            for sample_id, mutations in mutations_by_sample.items():
                for mutation in mutations:
                    try:
                        # Extract the relevant information from each field
                        chr_num = mutation['CHR']
                        cnv = mutation['CNV']
                        interval = mutation['INTERVAL']
                        kb = mutation['KB']
                        mid_bp = mutation['MID_BP']
                        targets = mutation['TARGETS']
                        num_targ = mutation['NUM_TARG']
                        q_start = mutation['Q_START']
                        q_stop = mutation['Q_STOP']
                        q_exact = mutation['Q_EXACT']
                        q_some = mutation['Q_SOME']
                        q_non_dip = mutation['Q_NON_DIPLOID']
                        mean_rd = mutation['MEAN_RD']
                        mean_orig_rd = mutation['MEAN_ORIG_RD']
                        
                        # Set reference allele to 'N'
                        ref = 'N'
                        # Set alternative allele based on CNV
                        alt = '<' + cnv + '>'
                        
                        if chr_num != interval.split(':')[0]:
                            # Raise an error if the chromosome number does not match the chromosome number in the interval
                            raise ValueError(f"Chromosome number in interval {interval} does not match the specified chromosome number {chr_num}")
                        
                        # Extract start and end positions from the interval field
                        start, end = map(int, interval.split(':')[1].split('-'))
                        
                        # Extract the length of the CNV
                        svlen = end - start + 1
                        
                        # Placeholder for strand information as neither CANOES, CLAMMS, nor XHMM provide it
                        strands = '.'
                        
                        # Placeholder for genotype informationas XHMM does not provide genotype information
                        gt = './.'
                        
                        # Set filter status depending on Q_SOME, Q_EXACT and Q_NON_DIPLOID
                        if q_some >= 60 and q_exact >= 60 and q_non_dip >= 60:
                            filter_status = 'PASS'
                        else:
                            filter_status = 'LowQuality'
                            
                        # Write the mutation to the VCF file
                        vcf_file.write(f"chr{chr_num}\t{start}\t.\t{ref}\t{alt}\t.\t{filter_status}\tEND={end};Q_START={q_start};Q_STOP={q_stop};SVLEN={svlen};SVMETHOD=XHMM;SVTYPE=CNV;CN={cnv};SVTYPE=CNV;STRANDS={strands}\tGT:Q_SOME:Q_EXACT:Q_NON_DIPLOID")
                        
                        # Write 'FORMAT' values for each sample
                        for sid in mutations_by_sample:
                            if sid == sample_id:
                                vcf_file.write(f"\t{gt}:{q_some}:{q_exact}:{q_non_dip}")
                            else:
                                # If the sample ID does not exist in the dictionary, write '.' for each sample
                                vcf_file.write(f"\t.")
                        
                        vcf_file.write("\n")
                    
                    # Catches 'KeyError' exceptions that occur within the corresponding 'try' block, where 'e' is the specific error
                    except KeyError as e:
                        error_msg = f"Key error in mutation: {e}"
                        logging.error(error_msg)
                    
                    # Catches 'ValueError' exceptions that occur within the corresponding 'try' block, where 'e' is the specific error
                    except ValueError as e:
                        error_msg = f"Value error in mutation: {e}"
                        logging.error(error_msg)
                        
    # Catches 'IOError' exceptions that occur within the corresponding 'try' block, where 'e' is the specific error
    except IOError as e:
        error_msg = f"IO error: {e}"
        logging.error(error_msg)

# Main function execution
# Script will only run if executed and not run as a module
if __name__ == "__main__":
    
    # Create an argument parser object with a description of the script
    parser = argparse.ArgumentParser(description="Convert XHMM XCNV file to VCF file.")
    
    # Add arguments for input file, output directory, combined VCF file name and log file
    parser.add_argument("-i", "--input_file", required=True, help="Input XHMM XCNV file")
    parser.add_argument("-d", "--output_dir", help="Output directory for VCF files")
    parser.add_argument("-o", "--combined_vcf_name", required=True, help="Combined VCF file name")
    parser.add_argument("-l", "--log_file", help="Log file name")
    
    # Parse the command-line arguments and store them in args.
    args = parser.parse_args()
    
    # Call the function to convert the XHMM XCNV file to VCF file
    convert_xhmm_xcnv_to_vcf(args.input_file, args.output_dir, args.combined_vcf_name, args.log_file)
        
                                
   
      
    
        
    
           
                   
