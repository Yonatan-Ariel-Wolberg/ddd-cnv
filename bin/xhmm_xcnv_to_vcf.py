#!/usr/bin/python3

# Import argparse for parsing command-line arguments
import argparse

# Import csv for reading and writing CSV files
import csv

# Import os for file and directory operations
import os

# Import logging for logging errors and other messages
import logging



def setup_logging(log_file=None):
    """
    Sets up logging configuration.
    """
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

# Define function that converts XHMM XCNV file to VCF file
def convert_xhmm_xcnv_to_vcf(input_file, output_dir, combined_vcf_name, log_file):
    """
    Converts XHMM XCNV file to VCF file and writes it to a combined VCF file.

        Parameters:
    - input_file: Path to the input XHMM XCNV file.
    - output_dir: Directory where output files will be saved.
    - combined_vcf_name: Name of the combined VCF file.
    - log_file: Optional log file to record errors.
    """
    
    setup_logging(log_file)  # Setup logging
    
    try:
        # If no output directory is specified, use the current working directory
        if not output_dir:
            output_dir = os.getcwd()

        # Create the output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Dictionary to store mutations grouped by SAMPLE ID
        mutations_by_sample = {}

        # List to keep track of the order of sample IDs
        sample_order = []

        try:
            # Open input file in read mode
            with open(input_file, 'r') as csvfile:
                # Create a CSV reader object
                csvreader = csv.DictReader(csvfile, delimiter='\t')

                # Print CSV headers for debugging
                print("CSV headers:", csvreader.fieldnames)

                # List of required columns in the CSV file
                required_columns = ['SAMPLE', 'CNV', 'INTERVAL', 'KB', 'CHR', 'MID_BP', 'TARGETS', 'NUM_TARG', 'Q_EXACT', 'Q_SOME', 'Q_NON_DIPLOID', 'Q_START', 'Q_STOP', 'MEAN_RD', 'MEAN_ORIG_RD']

                missing_columns = [col for col in required_columns if col not in csvreader.fieldnames]  # Check if all required columns are present
                if missing_columns:
                    # If columns are missing, raise an error
                    raise ValueError(f"Missing required columns: {missing_columns}")

                # Read each row in XCNV file
                for row_number, row in enumerate(csvreader, start=1):
                    try:
                        # Creates a dictionary called 'mutation' and stores the values from the current row under  different keys
                        mutation = {
                            'SAMPLE': row['SAMPLE'],
                            'CNV': row['CNV'],
                            'INTERVAL': row['INTERVAL'],
                            'KB': float(row['KB']),
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

                        # Extracts sample identifier from 'mutation' dictionary and assigns it to variable 'sample_id'
                        sample_id = mutation['SAMPLE']
                        # Line appends 'sample_id' to 'sample_order'
                        sample_order.append(sample_id)

                        # Line checks whether 'sample_id' is already a key in 'mutations_by_sample'
                        # If not, it creates a new key and assigns an empty list to it
                        if sample_id not in mutations_by_sample:
                            mutations_by_sample[sample_id] = []

                        # Line appends 'mutation' dictionary to the list of mutations for the current sample
                        mutations_by_sample[sample_id].append(mutation)

                    # Catches 'KeyError' execptions that occur within corresponding 'try' block, where 'e' is the specific 'KeyError'
                    except KeyError as e:
                        print(f"Key error in row {row_number}: {e} - Skipping row: {row}")

                    # Catches 'ValueError' execptions that occur within corresponding 'try' block, where 'e' is the specific 'ValueError'    
                    except ValueError as e:
                        print(f"Value error in row {row_number}: {e} - Skipping row: {row}")

    # Catches 'IOError' execptions that occur within corresponding 'try' block, where 'e' is the specific 'IOError'    
    except IOError as e:
        # Log the error and exit the program
        logging.error(f"Error reading {input_file}: {e}")
        return

    except ValueError as e:
        # Log the error and exit the program
        logging.error(f"Value error in {input_file}: {e}")
        return

    except Exception as e:
        # Log the error and exit the program
        logging.error(f"Unexpected error in {input_file}: {e}")
        return

    # Write the combined VCF file if specified
    if combined_vcf_name:
        combined_vcf_file = os.path.join(output_dir, combined_vcf_name)
        write_vcf_file(combined_vcf_file, mutations_by_sample, log_file)


    # Write separate VCF files for each sample
    for sample_id in mutations_by_sample:
        sample_vcf_file = os.path.join(output_dir, sample_id + ".vcf")
        write_vcf_file(sample_vcf_file, {sample_id: mutations_by_sample[sample_id]}, log_file)

    # Print the number of samples in the combined VCF file
    if combined_vcf_name:
        print(f"Number of samples in combined VCF file: {len(sample_order)}")
    else:
        print(f"Number of samples in separate VCF files: {len(sample_order)}")

    return


# Defines a function that writes all mutations to a single VCF file
def write_vcf_file(output_file, mutations_by_sample, log_file=None):
    try:    
        # Open output file for writing
        with open(output_file, 'w') as vcf_file:
            # Write header
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

            # Write the VCF header row with sample IDs
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

            sample_ids = list(mutations_by_sample.keys())

            for sample_id in sample_ids:
                vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # Write mutations for each sample
            for sample_id in sample_ids:
                mutations = mutations_by_sample[sample_id]
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
                        q_non_diploid = mutation['Q_NON_DIPLOID']
                        mean_rd = mutation['MEAN_RD']
                        mean_orig_rd = mutation['MEAN_ORIG_RD']

                        # Check if the chromosome number is the same as the number before the ':' in the interval field
                        if chr_num != interval.split(':')[0]:
                            # Raise an error if the chromosome number does not match the chromosome number in interval
                            raise ValueError("Chromosome number does not match chromosome number in interval field")

                        # Extract the start and end positions from the interval field
                        start, end = map(int, interval.split(':')[1].split('-'))
                        svlen = end - start + 1
                        ci_pos = (0, 0)
                        ci_end = (0, 0)
                        strands = '.'

                        # Filters variants that have Q_NON_DIPLOID, Q_EXACT and Q_SOME greater than or equal to 60
                        if float(q_non_diploid) >= 60 and float(q_exact) >= 60 and float(q_some) >= 60 :
                            filter_status = 'PASS'
        
                        else:
                            filter_status = 'LowQuality'

                        # Write the mutation information to the VCF file
                        vcf_file.write(f"chr{chr_num}\t{start}\t.\tN\t{cnv}\t.\t{filter_status}\tEND={end};SVLEN={svlen};SVMETHOD=XHMM;SVTYPE=CNV;CIPOS={ci_pos};CIEND={ci_end};STRANDS={strands}\tGT:Q_SOME:Q_EXACT:Q_NON_DIPLOID")