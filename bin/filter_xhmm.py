#!/usr/bin/python3

# Import argparse for parsing command-line arguments
import argparse

# Import os for file and directory operations
import os


def parse_args():
    """
    Parse command-line arguments
    """

    # Create an argument parser
    parser = argparse.ArgumentParser(description="Filter XHMM output VCF file based on SQ, EQ and NDQ values in sample fields.")

    # Add command-line arguments for input and output files and output directory
    parser.add_argument("-i", "--input", required=True, help="Unfiltered Input VCF file")
    parser.add_argument("-d", "--directory", help="Directory for Output VCF files")
    parser.add_argument("-o", "--output", required=True, help="Filtered Output VCF file")
    
    # Parse the command-line arguments
    args =parser.parse_args()

    if args.directory:
        # If an output directory is provided, create it if it doesn't exist
        if not os.path.exists(args.directory):
            os.makedirs(args.directory)

    else:
        # If no output directory is provided, use the current working directory
        args.directory = os.getcwd()

    return args

def filter_vcf(input_vcf, output_vcf):
    """
    Filters VCF file based on SQ, EQ and NDQ values in sample fields 
    """

    # Open the input and output VCF files
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            # Write header lines to output file
            if line.startswith('#'):
                outfile.write(line)
            
            else:
                # Split the line into fields
                fields = line.strip().split('\t')
                format_field = fields[8] # FORMAT field
                sample_fields = fields[9:] # Sample fields

            # Split the FORMAT field into a list of format labels
            format_labels = format_field.split(':')

            # Get the indices for NDQ, SQ and EQ in the Sample fields
            try:
                # Try to get the indices for NDQ, SQ and EQ
                ndq_index = format_labels.index('NDQ')
                eq_index = format_labels.index('EQ')
                sq_index = format_labels.index('SQ')

            except ValueError as e:
                # If an error occurs, print an error message and continue
                raise ValueError(f"Required format labels not found: {e}")
                continue

            # Initialize keep_line to False
            keep_line = False

            for sample in sample_fields:
                sample_data = sample.split(':') # Split the sample field into a list of sample values
                
                # Check if the sample has NDQ, EQ and SQ
                if len(sample_data) > max(ndq_index, eq_index, sq_index):
                   
                    # Extract NDQ, EQ and SQ values
                    ndq = float(sample_data[ndq_index]) if sample_data[ndq_index] != '.' else 0
                    eq = float(sample_data[eq_index]) if sample_data[eq_index] != '.' else 0
                    sq = float(sample_data[sq_index]) if sample_data[sq_index] != '.' else 0
                    
                    # Filter based on NDQ, EQ and SQ
                    if ndq >= 60 and eq >= 60 and sq >= 60:
                        # If all three values are >= 60, keep the line
                        keep_line = True
                        break
            
            # Write the line to the output file if it passes the filter
            if keep_line:
                outfile.write(line)

def main():
    """
    Main function to parse command-line arguments and call the filter_vcf function
    """

    # Parse command-line arguments
    args = parse_args()

    # Call the filter_vcf function
    filter_vcf(args.input, args.output)


# Call the main function
if __name__ == "__main__":
    main()
        