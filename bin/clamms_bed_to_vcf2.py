#!/usr/bin/env python3

# Import necessary modules
import argparse  # For parsing command-line arguments
import os  # For file and directory operations
import logging  # For logging information and errors

def setup_logging():
    """
    Set up logging configuration.
    """
    # Set logging level to 'INFO' and define the log message format
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def validate_input_file(input_file):
    """
    Validate that the input file exists and is readable.
    """
    if not os.path.isfile(input_file):
        # Log an error and raise FileNotFoundError if the file does not exist
        logging.error(f"Input file '{input_file}' not found.")
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

def create_output_directory(output_dir):
    """
    Create the output directory if it doesn't exist.
    """
    if not os.path.exists(output_dir):
        # Create the directory if it doesn't exist
        os.makedirs(output_dir)

def convert_clamms_bed_to_vcf(input_file, output_dir):
    """
    Convert CLAMMS BED format to VCF format and create both unified and separate VCF files.
    """
    try:
        # Ensure the output directory exists
        create_output_directory(output_dir)

        # Dictionaries to store mutations by sample and a list for all mutations
        mutations_by_sample = {}
        all_mutations = []

        # Open and read the input BED file
        with open(input_file, 'r') as bedfile:
            for line_number, line in enumerate(bedfile, start=1):
                try:
                    # Split the line into fields based on tab delimiter
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    interval = fields[3]
                    sample_id = fields[4]
                    cnv_type = fields[5]
                    mlcn = int(fields[6])
                    num_windows = int(fields[7])
                    q_some = float(fields[8])
                    q_exact = float(fields[9])

                    # Determine the filter status based on quality metrics
                    if q_exact > 0 and q_some >= 500:
                        filter_status = 'PASS'
                    else:
                        filter_status = 'LowQuality'

                    # Create a mutation dictionary with the extracted fields
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

                    # Add mutation to the dictionary for the corresponding sample ID
                    if sample_id not in mutations_by_sample:
                        mutations_by_sample[sample_id] = []
                    mutations_by_sample[sample_id].append(mutation)
                    all_mutations.append(mutation)

                except IndexError:
                    # Log an error if there are not enough fields in the line
                    logging.error(f"Error: Insufficient fields in line {line_number}: {line}")
                except (ValueError, TypeError) as e:
                    # Log an error if there is an invalid field value
                    logging.error(f"Error: Invalid field value in line {line_number}: {e}")

        # Write the unified VCF file with all mutations
        if all_mutations:
            write_vcf_file(output_dir, "combined.vcf", all_mutations, include_samples=True)
        else:
            # Log a warning if no mutations were found
            logging.warning("No mutations found in the input file.")

        # Write separate VCF files for each sample
        if mutations_by_sample:
            write_vcf_file_per_sample(output_dir, mutations_by_sample)

    except IOError as e:
        # Log an error if there is an I/O issue with file operations
        logging.error(f"Error reading or writing file: {e}")
    except Exception as e:
        # Log any unexpected errors
        logging.error(f"Unexpected error: {e}")

def write_vcf_file(output_dir, filename, mutations, include_samples=False):
    """
    Write mutations to a VCF file.

    Args:
    - output_dir (str): Directory to save the VCF file.
    - filename (str): Name of the VCF file.
    - mutations (list): List of mutations to write.
    - include_samples (bool): Whether to include all samples in the unified VCF file.
    """
    try:
        # Open the VCF file for writing
        output_file = os.path.join(output_dir, filename)
        with open(output_file, 'w') as vcf_file:
            # Write VCF headers
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
            
            # Write column headers
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            
            # Write sample IDs as column headers if creating a unified VCF
            if include_samples:
                sample_ids = set(m['SAMPLE_ID'] for m in mutations)
                for sample_id in sample_ids:
                    vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # Write each mutation to the VCF file
            for mutation in mutations:
                try:
                    chrom = mutation['CHROM']
                    start = mutation['START']
                    end = mutation['END']
                    cnv_type = mutation['CNV']
                    mlcn = mutation['MLCN']
                    q_some = mutation['Q_SOME']
                    q_exact = mutation['Q_EXACT']
                    filter_status = mutation['FILTER_STATUS']
                    svlen = end - start
                    ci_pos = (0, 0)
                    strands = '.'
                    svtype = 'CNV' if cnv_type in ['DEL', 'DUP'] else cnv_type

                    if mlcn == 3:
                        format_value = '0/1'
                    elif mlcn == 1:
                        format_value = '0/1'
                    elif mlcn == 2:
                        format_value = '0/0'
                    elif mlcn == 0:
                        format_value = '1/1'
                    else:
                        format_value = '.'

                    # Create VCF entry
                    line = f"chr{chrom}\t{start}\t.\tN\t<{cnv_type}>\t.\t{filter_status}\tEND={end};SVLEN={svlen};CN={mlcn};SVMETHOD=CLAMMS;SVTYPE={svtype};CIPOS={ci_pos};STRANDS={strands}\tGT:Q_SOME:Q_EXACT"
                    
                    # Add FORMAT values for each sample if creating a unified VCF
                    if include_samples:
                        sample_ids = set(m['SAMPLE_ID'] for m in mutations)
                        for sample_id in sample_ids:
                            if sample_id == mutation['SAMPLE_ID']:
                                line += f"\t{format_value}:{q_some}:{q_exact}"
                            else:
                                line += "\t0/0:0:0"
                    else:
                        line += f"\t{format_value}:{q_some}:{q_exact}"

                    # Write the VCF entry
                    vcf_file.write(line + "\n")

                except KeyError as e:
                    logging.error(f"Key error in mutation: {e}")
                except ValueError as e:
                    logging.error(f"Value error in mutation: {e}")

    except IOError as e:
        logging.error(f"Error reading or writing file: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")

def write_vcf_file_per_sample(output_dir, mutations_by_sample):
    """
    Write separate VCF files for each sample.
    
    Args:
    - output_dir (str): Directory to save the VCF files.
    - mutations_by_sample (dict): Dictionary with sample IDs as keys and lists of mutations as values.
    """
    for sample_id, mutations in mutations_by_sample.items():
        # Write a VCF file for each sample
        output_file = os.path.join(output_dir, f"{sample_id}.vcf")
        write_vcf_file(output_dir, f"{sample_id}.vcf", mutations)

if __name__ == "__main__":
    # Set up logging
    setup_logging()

    # Create argument parser
    parser = argparse.ArgumentParser(description="Convert CLAMMS BED file to unified VCF and separate VCF files per sample.")
    
    # Define command-line arguments
    parser.add_argument("input_file", help="Input CLAMMS BED file")
    parser.add_argument("output_dir", help="Output directory for VCF files")

    # Parse command-line arguments
    args = parser.parse_args()

    # Validate input file existence
    validate_input_file(args.input_file)

    # Convert BED file to VCF and create output files
    convert_clamms_bed_to_vcf(args.input_file, args.output_dir)