#!/usr/bin/python

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


def convert_xhmm_xcnv_to_vcf(input_file, output_dir, combined_vcf_name, log_file, sample_file=None, fai_file=None):
    """
    Function to convert XHMM XCNV files to VCF format.

    Parameters:
    - input_file: The path to the input XCNV file.
    - output_dir: The directory to save the output VCF file.
    - combined_vcf_name: The name of the combined VCF file.
    - log_file: Optional log file to record errors.
    - sample_file: Optional file containing sample IDs, one per line.
    - fai_file: FASTA index file for the reference genome.
    
    Returns:
    - The path to the combined VCF file.
    """
    # Set up logging if a log file is specified
    if log_file:
        setup_logging(log_file)

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        
        # Handles any 'Exception' errors while creating the directory
        except Exception as e:
            logging.error(f"Error creating directory {output_dir}: {e}")
            return
        
    # Read sample IDs from file if provided
    sample_list = read_sample_list(sample_file) if sample_file else []

    # Dictionary to store mutations grouped by sample ID
    mutations_by_sample = {}

    # List to keep track of order of sample IDs
    sample_order = []

    try:
        # Open the input CSV file
        with open(input_file, 'r') as bed_file:
            
            for line_number, line in enumerate(bed_file):
                try:
                    fields = line.strip().split('\t')
                    
                    # Extracts sample identifier from the row
                    sample = fields[4]

                    # Checks if the sample ID is in the sample list
                    if sample not in sample_list:
                        # If the sample ID is not in the sample list, skip the row
                        continue

                    # Creates a dictionary called 'mutation' and stores values from the current row under different keys
                    mutation = {
                         'CHROM': fields[0],
                            'START': safe_int(fields[1]),
                            'END': safe_int(fields[2]),
                            'INTERVAL': fields[3],
                            'SAMPLE_ID': sample,
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

                    # Appends 'sample' to 'sample_order' list
                    if sample not in sample_order:
                         sample_order.append(sample)

                    # Line checks if 'sample' is already in the 'mutations_by_sample' dictionary
                    if sample in mutations_by_sample:
                        # If 'sample' is already in the 'mutations_by_sample' dictionary, append 'mutation' to the corresponding list
                        mutations_by_sample[sample].append(mutation)
                    else:
                        # If 'sample' is not in the 'mutations_by_sample' dictionary, create a new list and append 'mutation' to it
                        mutations_by_sample[sample] = [mutation]

                
                # Catches any 'KeyError' errors within corresponding 'try' block, where 'e' is the error object
                except KeyError as e:
                    error_msg = f"Error processing row {line_number}: {e} - Skipping row: {line}"
                    logging.error(error_msg)

                # Catches any 'ValueError' errors within corresponding 'try' block, where 'e' is the error object 
                except ValueError as e:
                     error_msg = f"Error processing row {line_number}: {e} - Skipping row: {line}"
                     logging.error(error_msg)  


    # Catches any 'IOError' errors while reading the XCNV file, where 'e' is the error object
    except IOError as e:
            error_msg = f"Error readding {input_file}: {e}"
            logging.error(error_msg)
            return

    # Catches any 'ValueError' errors while reading the XCNV file, where 'e' is the error object
    except ValueError as e:
            error_msg = f"ValueError in {input_file}: {e}"
            logging.error(error_msg)
            return
    
    # Catches any 'Exception' errors while reading the XCNV file, where 'e' is the error object
    except Exception as e:
            error_msg = f"Unexpected error in {input_file}: {e}"
            logging.error(error_msg)
            return
    
    # Defines 'sample_order' as a list of sample IDs in 'sample_list' that are also in 'mutations_by_sample'
    sample_order = [sample for sample in sample_list if sample in mutations_by_sample]
            
    # Check for duplicates in sample_order
    if len(sample_order) != len(set(sample_order)):
        logging.error("Duplicates found in sample_order!")
        print("Duplicates found in sample_order!")
        return


    # Write combined VCF file
    if combined_vcf_name:
        # Create combined VCF file name
        combined_vcf_file = os.path.join(output_dir, combined_vcf_name)
        
        # Extract reference name from FASTA index file
        ref_name = extract_ref_name(fai_file)

        # Create VCF contig lines
        sorted_contig_lines = create_vcf_contig_lines(fai_file)
        
        # Write combined VCF file
        write_vcf_file(combined_vcf_file, ref_name, sorted_contig_lines, mutations_by_sample, sample_order, log_file)

    # Write individual VCF files for each sample
    for sample in sample_order:
        if sample not in mutations_by_sample:
            continue
        # Create individual VCF file name after sample ID
        individual_vcf_file = os.path.join(output_dir, f"{sample}.vcf")
        
        # Write individual VCF file
        write_vcf_file(individual_vcf_file, ref_name, sorted_contig_lines, {sample: mutations_by_sample[sample]}, [sample], log_file)

    # Print number of samples in combined VCF file or in indivivdual VCF files
    if combined_vcf_name:
        print(f"Number of samples in combined VCF file: {len(sample_order)}")
    else:
        print(f"Number of samples in individual VCF files: {len(sample_order)}")

def create_vcf_contig_lines(fai_file):
    """
    Create VCF contig lines from FASTA index file

    Parameters:
        fai_file (str): Path to FASTA index file

    Returns:
        list: List of VCF contig lines
    """
    # List of VCF contig lines
    contig_lines = []

    # List of valid contig names
    valid_contigs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 
                     'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY' 
                    ]
    
    # Set of valid contig names
    valid_contig_set = set(valid_contigs)

    # Read FASTA index file
    with open(fai_file) as fai:
        for line in fai:
            # Split line by tab and remove newline and whitespace
            parts = line.strip().split("\t")
            contig_name = parts[0]
            contig_length = int(parts[1])

            if contig_name in valid_contig_set:
                contig_line = f"##contig=<ID={contig_name},length={contig_length}>"
                contig_lines.append((valid_contigs.index(contig_name), contig_line))

    # Sort contig lines based on the index in 'valid_contigs'
    contig_lines.sort()

    # Create a list of contig lines, sorted by contig name in 'valid_contigs'
    sorted_contig_lines = [line for _, line in contig_lines]
    
    # Return list of VCF contig lines
    return sorted_contig_lines

def extract_ref_name(fai_file):
    """
    Extract reference name from FASTA index file

    Parameters:
        fai_file (str): Path to FASTA index file

    Returns:
        str: Reference name
    """
    # Get base name of FASTA index file
    base_name = os.path.basename(fai_file)

    # Remove last file extension
    base_name_no_ext = base_name.rsplit('.', 1)[0]
    # Remove second last file extension
    base_name_no_ext2 = base_name_no_ext.rsplit('.', 1)[0]
    # Remove third last file extension
    ref_name = base_name_no_ext2.rsplit('.', 1)[0]

    return ref_name

def write_vcf_file(output_file, ref_name, sorted_contig_lines, mutations_by_sample, sample_order, log_file=None):
    """
    Write mutations to a VCF file

    Parameters:
        - output_file (str): Output VCF file path
        - fai_file (str): Path to FASTA index file
        - sorted_contig_lines (list): List of VCF contig lines
        - mutations_by_sample (dict): Dictionary of mutations by sample
        - sample_order (list): List of sample IDs in the same order as mutations_by_sample
        - log_file (str, optional): Path to log file. Defaults to None.
    """
    
    # Setup logging if log_file is provided
    if log_file:
        setup_logging(log_file)

    try:
        # Open output VCF file for writing
        with open(output_file, 'w') as vcf_file:

            # Write header to VCF file
            vcf_file.write("##fileformat=VCFv4.1\n")
            vcf_file.write("##fileDate=" + datetime.now().strftime("%d%m%Y") + "\n")
            vcf_file.write("##source=XHMM-XCNVToVcfConverter\n")
            # Write reference to VCF file
            vcf_file.write(f"##reference={ref_name}\n")
            vcf_file.write("##phasing=partial\n")

            # Write contig lines to VCF file
            for line in sorted_contig_lines:
                vcf_file.write(line + "\n")

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
            
            # Write sample-specific fields after the FORMAT field
            for sample_id in sample_order:
                vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # For each sample in sample_order, write mutations
            for sample_id in sample_order:
                mutations = mutations_by_sample[sample_id]
                for mutation in mutations:
                    try:
                        # Extract relevant information from each field
                        chr_num = mutation['CHROM']
                        cnv = mutation['CNV_TYPE']
                        interval = mutation['INTERVAL']
                        start = mutation['START']
                        end = mutation['END']
                        mlcn = mutation['MLCN']
                        num_windows = mutation['NUM_WINDOWS']
                        q_some = mutation['Q_SOME']
                        q_exact = mutation['Q_EXACT']
                        q_left_extend = mutation['Q_LEFT_EXTEND']
                        left_extend_coord = mutation['LEFT_EXTEND_COORD']
                        q_right_extend = mutation['Q_RIGHT_EXTEND']
                        right_extend_coord = mutation['RIGHT_EXTEND_COORD']
                        q_left_contract = mutation['Q_LEFT_CONTRACT']
                        left_contract_coord = mutation['LEFT_CONTRACT_COORD']
                        q_right_contract = mutation['Q_RIGHT_CONTRACT']
                        right_contract_coord = mutation['RIGHT_CONTRACT_COORD']

                        # Set the reference allele to 'N'
                        ref = 'N'

                        # Set the alternate allele to the CNV type
                        alt = '<' + cnv + '>'

                        # Checks if the chromosome number in the interval matches the specified chromosome number
                        if chr_num != interval.split(':')[0]:
                            # Raise an error if the chromosome number in the interval does not match the specified chromosome number
                            raise ValueError(f"Chromosome number in interval {interval} does not match specified chromosome number {chr_num}")

                        # SVLEN is the length of the SV
                        svlen = end - start + 1

                        # Set filter status based on quality scores
                        if q_some >= 500 and q_exact > 0:
                            filter_status = 'PASS'
                        else:
                            filter_status = 'LowQuality'
                        
                        # Placeholder for strand information as XHMM does not provide strand information
                        strand = '.'

                        # Placeholder for genotype information as XHMM does not provide genotype information
                        gt = './.'

                        # Write the mutation to the VCF file
                        vcf_file.write(f"chr{chr_num}\t{start}\t.\t{ref}\t{alt}\t.\t{filter_status}\tEND={end};SVLEN={svlen};SVTYPE={cnv};SVMETHOD=CLAMMS;STRANDS={strand};NUM_WINDOWS={num_windows};MLCN={mlcn};Q_LEFT_EXTEND={q_left_extend};Q_RIGHT_EXTEND={q_right_extend};LEFT_EXTEND_COORD={left_extend_coord};RIGHT_EXTEND_COORD={right_extend_coord};Q_LEFT_CONTRACT={q_left_contract};LEFT_CONTRACT_COORD={left_contract_coord};Q_RIGHT_CONTRACT={q_right_contract};RIGHT_CONTRACT_COORD={right_contract_coord}\tGT:Q_SOME:Q_EXACT")

                        # Write 'FORMAT' values for each sample
                        for sample_id in sample_order:
                            mutation_found = False
                            if sample_id in mutations_by_sample:
                                # Get list of mutations for current sample
                                sample_mutations = mutations_by_sample[sample_id]
                                # Check if the current mutation is in the list of mutations for the current sample
                                for sample_mutation in sample_mutations:
                                    # Check if the current mutation is the same as the corresponding mutation in the list of mutations
                                    if sample_mutation == mutation:
                                        vcf_file.write(f"\t{gt}:{q_some if q_some else '.'}:{q_exact if q_exact else '.'}")
                                        mutation_found = True
                                        break
                            if not mutation_found:
                                vcf_file.write("\t./.:.:.")
                        vcf_file.write("\n")
                    
                    # Catches 'KeyError' exceptions in the corresponding 'try' block, where 'e' is specifically the 'KeyError' object
                    except KeyError as e:
                        logging.error(f"KeyError when writing mutation {mutation}: {e}")

                    # Catches 'ValueError' exceptions in the corresponding 'try' block, where 'e' is specifically the 'ValueError' object
                    except ValueError as e:
                        logging.error(f"ValueError when writing mutation {mutation}: {e}")

    # Catches 'IOError' exceptions in the corresponding 'try' block, where 'e' is specifically the 'IOError' object               
    except IOError as e:
        logging.error(f"Error writing VCF to {output_file}: {e}")

    # Catches any other exceptions in the corresponding 'try' block, where 'e' is the exception object
    except Exception as e:
        logging.error(f"Unexpected error writing VCF to {output_file}: {e}")


def main():
    """
    Main function of the script.
    """

    # Create an argument parser with a description of the script
    parser = argparse.ArgumentParser(description="Convert XHMM XCNV output to VCF format")

    # Add arguments to the parser for the input file, output directory, combined VCF file name, log file and sample_file
    parser.add_argument('-i', '--input', type=str, help='Input XCNV file', required=True)
    parser.add_argument('-d', '--output_dir', type=str, help='Output directory for VCF files')
    parser.add_argument('-c', '--combined_vcf_name', type=str, help='Name for the combined VCF file')
    parser.add_argument('-l', '--log_file', type=str, help='Name for the log file to record errors and issues')
    parser.add_argument('-s', '--sample_file', type=str, help='File containing sample IDs')
    parser.add_argument('-f', '--fai_file', type=str, help='FASTA index file for the reference genome')

    # Parse the arguments and store them in args
    args = parser.parse_args()

    # Call the convert_xhmm_xcnv_to_vcf function with the parsed arguments
    convert_xhmm_xcnv_to_vcf(args.input, args.output_dir, args.combined_vcf_name, args.log_file, args.sample_file, args.fai_file)

# Main function execution
# Script will only run if exectuted and not run as a module
if __name__ == "__main__":
    main()