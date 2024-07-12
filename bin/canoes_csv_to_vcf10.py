#!/usr/bin/env python3

# Library used for parsing command-line arguments
import argparse

# Library used for reading and writing CSV files
import csv

# Library used for checking and creating directories if they don't exist
import os

# Defines a function that converts CANOES output CSV file to VCF file
def convert_canoes_csv_to_vcf(input_file, output_file):
    try:
        # Create the directory if it doesn't exist
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Dictionary to store mutations grouped by SAMPLE ID
        mutations_by_sample = {}

        # Open input file in read mode
        with open(input_file, 'r') as csvfile:
            # Create a CSV reader object
            csvreader = csv.DictReader(csvfile, delimiter='\t')

            # Print CSV headers for debugging
            print("CSV headers:", csvreader.fieldnames)
            
            # Check if all required columns are present
            required_columns = ['SAMPLE', 'CNV', 'INTERVAL', 'KB', 'CHR', 'MID_BP', 'TARGETS', 'NUM_TARG', 'MLCN', 'Q_SOME']
            
            # List of columns that are actually preseht in CSV file and stores column names in 'fieldnames'
            missing_columns = [col for col in required_columns if col not in csvreader.fieldnames]
            if missing_columns:
                # If some columns are missing, raise an error
                raise ValueError(f"Missing required columns: {missing_columns}")

            # Read each row from the CSV file
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
                        'MLCN': row['MLCN'],
                        'Q_SOME': row['Q_SOME']
                    }

                    # Extracts sample identifier from 'mutation' dictinary and assigns it to variable 'sample_id'
                    sample_id = mutation['SAMPLE']
                    
                    # Line checks whether 'sample_id' is already a key in 'mutations_by_sample'
                    # If not, it creates a new key and assigns an empty list to it
                    if sample_id not in mutations_by_sample:
                        mutations_by_sample[sample_id] = []

                    # Append mutation to the corresponding sample list
                    mutations_by_sample[sample_id].append(mutation)

                # Catches 'KeyError' execptions that occur within corresponding 'try' block, where 'e' is the specific 'KeyError'
                except KeyError as e:
                    print(f"Key error in row {row_number}: {e} - Skipping row: {row}")
                
                # Catches 'ValueError' execptions that occur within corresponding 'try' block, where 'e' is the specific 'ValueError'
                except ValueError as e:
                    print(f"Value error in row {row_number}: {e} - Skipping row: {row}")

        # Write all mutations to a single VCF file maintaining sample order
        # Checks 'mutations_by_sample' dictionary is not empty
        if mutations_by_sample:
            # If 'mutations_by_sample' is not empty, write all mutations to a single VCF file
            write_vcf_file(output_file, mutations_by_sample)
        else:
            print("No valid mutations found in the input file.")

    # Catches 'IOError' execptions that occur within corresponding 'try' block, where 'e' is the specific 'IOError'
    except IOError as e:
        print(f"Error reading or writing file: {e}")
    
    # Catches 'ValueError' execptions that occur within corresponding 'try' block, where 'e' is the specific 'ValueError'
    except ValueError as e:
        print(f"Error: {e}")
   
    # Catches any other unexpected execptions that occur within corresponding 'try' block, where 'e' is the unexpected issue
    except Exception as e:
        print(f"Unexpected error: {e}")


# Defines a function that writes all mutations to a single VCF file
def write_vcf_file(output_file, mutations_by_sample):
    try:
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
            vcf_file.write("##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of the SV\">\n")
            vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
            vcf_file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n")
            vcf_file.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n")
            vcf_file.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint\">\n")
            vcf_file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf_file.write("##FORMAT=<ID=Q_SOME,Number=1,Type=Float,Description=\"Quality of the sample\">\n")
            
            # Write the sample-specific fields after the FORMAT field
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

            # Write sample IDs as column headers in the order they appeared in the CSV file
            for sample_id in mutations_by_sample:
                vcf_file.write(f"\t{sample_id}")
            vcf_file.write("\n")

            # Process mutations for each sample in order
            for sample_id in mutations_by_sample:
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

                        elif cn == '1':
                            # If CN = 1, this indicates a heterozygous deletion
                            format_values = '0/1'
                        
                        elif cn == '2':
                            # If CN = 2, this likely indicates a diploid condition, but it could represent a heterozygous deletion and duplication
                            format_values = '0/0'

                        elif cn == '0':
                            # If CN = 0, this indicates a homozygous deletion
                            format_values = '1/1'
                        
                        else:
                            # If CN > or = to 4, this could indicate a heterozygous triplication or a homozygous duplication
                            format_values = '.'
                        
                        # Write the VCF entry with 'chr' prefix
                        vcf_file.write(f"chr{chr_num}\t{start}\t.\tN\t{cnv}\t.\t{filter_status}\tEND={end};SVLEN={svlen};CN={cn};SVMETHOD=CANOES;SVTYPE=CNV;CIPOS={ci_pos};CIEND={ci_end};STRANDS={strands}\tGT:Q_SOME")
                        
                        # Write FORMAT values for each sample
                        for sid in mutations_by_sample:
                            # Check if the sample ID matches the current sample ID
                            if sid == sample_id:
                                # Writes genotype (format_values) followed by the quality (q_some)
                                # If q_some is 'None', it writes '0'
                                vcf_file.write(f"\t{format_values}:{q_some if q_some else 0}")
                           
                            else:
                                # If the sample ID does not match the current sample ID, it writes '0/0' for genotype as sample is assumed to be homozygous for no mutation and quality score is '0'
                                vcf_file.write("\t0/0:0")

                        vcf_file.write("\n")

                    # Catches 'KeyError' execptions that occur within corresponding 'try' block, where 'e' is the specific error
                    except KeyError as e:
                        print(f"Key error in mutation: {e}")
                    
                    # Catches 'ValueError' execptions that occur within corresponding 'try' block, where 'e' is the specific error
                    except ValueError as e:
                        print(f"Value error in mutation: {e}")

    # Catches 'IOError' execptions that occur within corresponding 'try' block, where 'e' is the specific error
    except IOError as e:
        print(f"Error writing VCF file: {e}")

# Main function execution
# Script will only run if executed and not run as a module
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
