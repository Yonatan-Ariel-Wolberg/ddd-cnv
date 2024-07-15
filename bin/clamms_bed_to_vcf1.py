#!/usr/bin/env python3

# Import 'argparse' module for parsing command-line arguments specifying input and output files
import argparse

# Import 'os' module for file and directory operations
import os

# Import 'logging' module to provide informative messages about the scripts exectution, such as errors, warnings and general information
import logging



def setup_logging():
    """
    Set up logging configuration
    """
    # Set logging level to 'INFO', ignoring 'DEBUG' messages
    # 'format' specifies the time when the log message is created, the severity level of the log and an actual log message provided by the code
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')



def validate_input_file(input_file):
    """
    Validate that the input file exists and is readable.
    """
    # Function checks if path specified by input file is an existing regular file
    if not os.path.isfile(input_file):
        
        # If the path doesn't exist or is not a regular file, raise an error
        logging.error(f"Input file '{input_file}' not found.")
        
        # Raise FileNotFoundError exception
        raise FileNotFoundError(f"Input file '{input_file}' not found.")



def create_output_directory(output_file):
    """
    Create the output directory if it doesn't exist.
    """
    # Returns the directory name of the path
    output_dir = os.path.dirname(output_file)

    # If the output directory doesn't exist, create it
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)



def convert_clamms_bed_to_vcf(input_file, output_file):
    """
    Convert CLAMMS BED format to VCF format.

    Args:
    - input_file (str): Path to the input CLAMMS BED file.
    - output_file (str): Path to the output VCF file.

    """
    try:
        # Create the output directory if it doesn't exist
        create_output_directory(output_file)

        # Dictionary to store mutations grouped by SAMPLE ID
        mutations_by_sample = {}

        # Open input file in read mode
        with open(input_file, 'r') as bedfile:
            # Read each line from the BED file
            for line_number, line in enumerate(bedfile, start=1):
                try:
                    # Split the line into fields
                    fields = line.strip().split('\t')
                    
                    # Extract fields from the BED format
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    interval = fields[3]  # Assuming the interval format is chr:start-end
                    sample_id = fields[4]
                    cnv_type = fields[5]
                    mlcn = int(fields[6])
                    num_windows = int(fields[7])
                    q_some = float(fields[8])
                    q_exact = float(fields[9])  # Add Q_EXACT field
                    
                    # Filter based on Q_EXACT and Q_SOME
                    if q_exact > 0 and q_some >= 500:
                        filter_status = 'PASS'
                    else:
                        filter_status = 'LowQuality'

                    # Initialize the mutation dictionary
                    mutation = {
                        'CHROM': chrom,
                        'START': start,
                        'END': end,
                        'INTERVAL': interval,
                        'SAMPLE_ID': sample_id,
                        'CNV': cnv_type,
                        'MLCN': mlcn,
                        'NUM_WINDOWS': num_windows,
                        'Q_SOME': q_some,
                        'Q_EXACT': q_exact,
                        'FILTER_STATUS': filter_status
                    }

                    # Check if sample ID is already in the dictionary
                    if sample_id not in mutations_by_sample:

                        # Initializes a new list for that sample ID
                        mutations_by_sample[sample_id] = []

                    # Append mutation to the corresponding sample list
                    mutations_by_sample[sample_id].append(mutation)


                # Catch any IndexError exceptions and log them
                except IndexError:
                    logging.error(f"Error: Insufficient fields in line {line_number}: {line}")
                
                # Catch any ValueError or TypeError exceptions and log them
                except (ValueError, TypeError) as e:
                    logging.error(f"Error: Invalid field value in line {line_number}: {e}")


        # Write all mutations to a single VCF file maintaining sample order
        if mutations_by_sample:
            write_vcf_file(output_file, mutations_by_sample)
        else:
            # If no mutations were found, log a warning
            logging.warning("No mutations found in the input file.")

    # Catch any IOError exceptions and log them
    except IOError as e:
        logging.error(f"Error reading or writing file: {e}")

    # Catch any Exception exceptions and log them
    except Exception as e:
        logging.error(f"Unexpected error: {e}")



def write_vcf_file(output_file, mutations_by_sample):
    """
    Write mutations grouped by SAMPLE ID to a VCF file.

    Args:
    - output_file (str): Path to the output VCF file.
    - mutations_by_sample (dict): Dictionary containing mutations grouped by SAMPLE ID.

    """
    try:
        # Open output file in write mode
        with open(output_file, 'w') as vcf_file:
            # Write VCF header
            vcf_file.write("##fileformat=VCFv4.1\n")
            vcf_file.write("##source=CLAMMSBedToVcfConverter\n")
            vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
            vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
            vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
            vcf_file.write("##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Tool used for calling SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n")
            vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strand orientation of the SV\">\n")
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Phred-scaled quality of the CNV call\">\n")
            vcf_file.write("##FORMAT=<ID=Q_EXACT,Number=1,Type=Float,Description=\"Non-Phred-scaled quality of the CNV call\">\n")
            
            # Write sample-specific fields after the FORMAT field
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

            # Write sample IDs as column headers in the order they appear in the BED file
            for sample_id in mutations_by_sample:
                vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # Process mutations for each sample in order
            for sample_id in mutations_by_sample:
                # For each sample ID, this line retrieves a list of mutations associated with that sample ID and assigns it to the variable 'mutations'
                mutations = mutations_by_sample[sample_id]
                
                # Process each mutation for the current sample
                for mutation in mutations:
                    try:
                        chrom = mutation['CHROM']
                        start = mutation['START']
                        end = mutation['END']
                        interval = mutation['INTERVAL']
                        sample_id = mutation['SAMPLE_ID']
                        cnv_type = mutation['CNV']
                        mlcn = mutation['MLCN']
                        num_windows = mutation['NUM_WINDOWS']
                        q_some = mutation['Q_SOME']
                        q_exact = mutation['Q_EXACT']
                        
                        # Calculate SV length
                        svlen = end - start
                        
                        # Calculate confidence intervals
                        ci_pos = (0, 0)  # Placeholder for now
                        strands = '.'    # Placeholder for now
                        
                        # Determine SV type based on CNV type
                        svtype = 'CNV' if cnv_type == 'DEL' or cnv_type == 'DUP' else cnv_type
                        
                        # Calculate FILTER status based on Q_EXACT and Q_SOME
                        if q_exact > 0 and q_some >= 500:
                            filter_status = 'PASS'
                        else:
                            filter_status = 'LowQuality'
                        
                        # Determine FORMAT value based on mlcn
                        if mlcn == 3:
                            format_value = '0/1'  # Example, adjust as per your logic
                        elif mlcn == 1:
                            format_value = '0/1'  # Example, adjust as per your logic
                        elif mlcn == 2:
                            format_value = '0/0'  # Example, adjust as per your logic
                        elif mlcn == 0:
                            format_value = '1/1'  # Example, adjust as per your logic
                        else:
                            format_value = '.'    # Example, adjust as per your logic
                        
                        # Write the VCF entry
                        vcf_file.write(f"{chrom}\t{start}\t.\tN\t<{svtype}>\t.\t{filter_status}\tEND={end};SVLEN={svlen};CN={mlcn};SVMETHOD=CLAMMS;SVTYPE={svtype};CIPOS={ci_pos};STRANDS={strands}\tGT:Q_SOME:Q_EXACT")
                        
                        # Write FORMAT values for each sample
                        for sid in mutations_by_sample:
                            if sid == sample_id:
                                vcf_file.write(f"\t{format_value}:{q_some}:{q_exact}")
                            else:
                                vcf_file.write("\t0/0:0:0")

                        vcf_file.write("\n")

                    # Catch any KeyError exceptions and log them
                    except KeyError as e:
                        logging.error(f"Key error in mutation: {e}")

                    # Catch any ValueError exceptions and log them
                    except ValueError as e:
                        logging.error(f"Value error in mutation: {e}")
   
    # Catch any IOError exceptions and log them
    except IOError as e:
        logging.error(f"Error reading or writing file: {e}")
   
    # Catch any other unexpected exceptions and log them
    except Exception as e:
        logging.error(f"Unexpected error: {e}")


# Main function
# Script must be executed directly
if __name__ == "__main__":
    # Calls 'setup_logging()' to set up logging configuration
    setup_logging()


    # Create an argument parser object with a description of the script
    parser = argparse.ArgumentParser(description="Convert CLAMMS BED file to VCF file.")
    
    # Add a positional argument for the input BED file
    parser.add_argument("input_file", help="Input CLAMMS BED file")
    
    # Add a positional argument for the output VCF file
    parser.add_argument("output_file", help="Output VCF file")
    
    # Parse the command-line arguments and store them in args.
    args = parser.parse_args()
    
    # Call the function to convert the BED file to VCF file
    convert_clamms_bed_to_vcf(args.input_file, args.output_file)

    # Validate the input file
    validate_input_file(args.input_file)

    # Call the function to convert the BED file to VCF file
    convert_clamms_bed_to_vcf(args.input_file, args.output_file)