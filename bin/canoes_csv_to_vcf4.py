#!/usr/bin/env python3

import argparse
import csv
import os
import logging

# Configure logging
def setup_logging(log_file=None):
    """
    Set up the logging configuration.
    
    Parameters:
    - log_file: Optional log file to record errors.
    """
    logging.basicConfig(
        filename=log_file,
        level=logging.ERROR,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def convert_canoes_csv_to_vcf(input_file, output_dir, combined_vcf_name, log_file=None):
    """
    Convert CANOES CSV file to VCF format.

    Parameters:
    - input_file: Path to the input CANOES CSV file.
    - output_dir: Directory where output files will be saved.
    - combined_vcf_name: Name of the combined VCF file.
    - log_file: Optional log file to record errors.
    """
    
    # Set up logging if a log file is specified
    if log_file:
        setup_logging(log_file)
    
    # Create the output directory if it does not exist
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except Exception as e:
            logging.error(f"Error creating directory {output_dir}: {e}")
            return

    # Dictionary to store mutations grouped by SAMPLE ID
    mutations_by_sample = {}

    try:
        # Open the input CSV file
        with open(input_file, 'r') as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')

            # Check if all required columns are present
            required_columns = ['SAMPLE', 'CNV', 'INTERVAL', 'KB', 'CHR', 'MID_BP', 'TARGETS', 'NUM_TARG', 'MLCN', 'Q_SOME']
            missing_columns = [col for col in required_columns if col not in csvreader.fieldnames]
            if missing_columns:
                raise ValueError(f"Missing required columns: {missing_columns}")

            # Process each row in the CSV file
            for row_number, row in enumerate(csvreader, start=1):
                try:
                    # Extract mutation information from the row
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

                    sample_id = mutation['SAMPLE']

                    # Initialize the list for this sample if not already present
                    if sample_id not in mutations_by_sample:
                        mutations_by_sample[sample_id] = []

                    # Append mutation to the list for this sample
                    mutations_by_sample[sample_id].append(mutation)

                except KeyError as e:
                    logging.error(f"Key error in row {row_number}: {e} - Skipping row: {row}")
                except ValueError as e:
                    logging.error(f"Value error in row {row_number}: {e} - Skipping row: {row}")

    except IOError as e:
        logging.error(f"Error reading input file {input_file}: {e}")
        return
    except ValueError as e:
        logging.error(f"Value error: {e}")
        return
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        return

    # Write the combined VCF file if specified
    if combined_vcf_name:
        combined_vcf_file = os.path.join(output_dir, combined_vcf_name) if output_dir else combined_vcf_name
        write_vcf_file(combined_vcf_file, mutations_by_sample, log_file)

    # Write separate VCF files for each sample
    for sample_id in mutations_by_sample:
        sample_vcf_file = os.path.join(output_dir, f"{sample_id}.vcf") if output_dir else f"{sample_id}.vcf"
        write_vcf_file(sample_vcf_file, {sample_id: mutations_by_sample[sample_id]}, log_file)

def write_vcf_file(output_file, mutations_by_sample, log_file=None):
    """
    Write mutations to a VCF file.

    Parameters:
    - output_file: Path to the output VCF file.
    - mutations_by_sample: Dictionary of mutations grouped by sample.
    - log_file: Optional log file to record errors.
    """

    try:
        # Open the VCF file for writing
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
                        chr_num = mutation['CHR']
                        cnv = mutation['CNV']
                        interval = mutation['INTERVAL']
                        kb = mutation['KB']
                        mid_bp = mutation['MID_BP']
                        targets = mutation['TARGETS']
                        num_targ = mutation['NUM_TARG']
                        cn = mutation['MLCN']
                        q_some = mutation['Q_SOME']

                        if chr_num != interval.split(':')[0]:
                            raise ValueError("chr_num does not match the interval field")

                        start, end = map(int, interval.split(':')[1].split('-'))
                        svlen = end - start
                        ci_pos = (0, 0)
                        ci_end = (0, 0)
                        strands = '.'

                        if kb >= 100 and int(q_some) >= 80:
                            filter_status = 'PASS'
                        elif kb < 100 and int(q_some) < 80:
                            filter_status = 'LowQuality'
                        else:
                            filter_status = '.'

                        if cn == '3' or cn == '1':
                            format_values = '0/1'
                        elif cn == '2':
                            format_values = '0/0'
                        elif cn == '0':
                            format_values = '1/1'
                        else:
                            format_values = '.'

                        vcf_file.write(f"chr{chr_num}\t{start}\t.\tN\t{cnv}\t.\t{filter_status}\tEND={end};SVLEN={svlen};CN={cn};SVMETHOD=CANOES;SVTYPE=CNV;CIPOS={ci_pos};CIEND={ci_end};STRANDS={strands}\tGT:Q_SOME")
                        
                        for sid in sample_ids:
                            if sid == sample_id:
                                vcf_file.write(f"\t{format_values}:{q_some if q_some else 0}")
                            else:
                                vcf_file.write("\t0/0:0")

                        vcf_file.write("\n")

                    except KeyError as e:
                        logging.error(f"Key error in mutation: {e}")
                    except ValueError as e:
                        logging.error(f"Value error in mutation: {e}")

    except IOError as e:
        logging.error(f"Error writing VCF file {output_file}: {e}")

# Main function execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CANOES CSV file to VCF file.")
    
    # Add arguments for input file, output directory, combined VCF file name, and log file
    parser.add_argument("-i", "--input_file", required=True, help="Input CANOES CSV file")
    parser.add_argument("-d", "--output_dir", help="Optional output directory")
    parser.add_argument("-o", "--combined_vcf_name", required=True, help="Name for the unified VCF file")
    parser.add_argument("-l", "--log_file", help="Optional log file")

    # Parse arguments
    args = parser.parse_args()
    
    # Convert CSV to VCF
    convert_canoes_csv_to_vcf(args.input_file, args.output_dir, args.combined_vcf_name, args.log_file)
