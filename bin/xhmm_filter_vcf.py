#!/usr/bin/env python3

import argparse

def parse_sample_field(format_field, sample_field):
    """
    Parse the sample field based on the format field.

    Args:
        format_field (str): The format field (e.g., 'GT:NDQ:DQ:EQ:SQ:NQ:LQ:RQ:PL:RD:ORD:DSCVR').
        sample_field (str): The sample field data (e.g., '0:0:81:0,0:0,0:81,93:0,0:0,0:0,255,255:-0.25:29.32:N').

    Returns:
        dict: A dictionary mapping each format key to its corresponding sample value.
    """
    # Split the format field into individual format keys
    format_keys = format_field.split(':')
    # Split the sample field into individual sample values
    sample_values = sample_field.split(':')
    
    # Create a dictionary mapping each format key to its corresponding sample value
    format_dict = {format_keys[i]: sample_values[i] for i in range(len(format_keys))}
    
    return format_dict

def filter_vcf(input_vcf, output_vcf, min_sq=60, min_eq=60, min_ndq=60):
    """
    Filter the VCF file based on specified criteria for SQ, EQ, and NDQ.

    Args:
        input_vcf (str): Path to the input VCF file.
        output_vcf (str): Path to the output VCF file.
        min_sq (int): Minimum SQ value to keep a sample (default is 60).
        min_eq (int): Minimum EQ value to keep a sample (default is 60).
        min_ndq (int): Minimum NDQ value to keep a sample (default is 60).
    """
    # Open the input VCF file for reading and the output VCF file for writing
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        # Iterate over each line in the input file
        for line in infile:
            # If the line starts with '#', it's a header line; write it as-is to the output file
            if line.startswith('#'):
                outfile.write(line)
            else:
                # Split the line into fields based on tab character
                fields = line.strip().split('\t')
                # Extract the format field (8th field in zero-indexed list)
                format_field = fields[8]
                # Extract all sample fields (starting from 9th field)
                sample_fields = fields[9:]
                
                # List to hold the processed sample fields
                new_samples = []
                
                # Process each sample field
                for sample_field in sample_fields:
                    # Parse the sample field based on the format field
                    sample_dict = parse_sample_field(format_field, sample_field)
                    
                    # Extract SQ, EQ, and NDQ values from the parsed sample dictionary
                    # If the value is not present, set it to 0
                    sq = int(sample_dict.get('SQ', 0))
                    eq = int(sample_dict.get('EQ', 0))
                    ndq = int(sample_dict.get('NDQ', 0))
                    
                    # Apply the filtering criteria
                    if sq >= min_sq and eq >= min_eq and ndq >= min_ndq:
                        new_samples.append(sample_field)
                    else:
                        # If the criteria are not met, set the sample to '0/0' with default values
                        new_samples.append('0/0:' + ':'.join(['0'] * (len(sample_dict) - 1)))
                
                # Write the line with potentially modified sample fields to the output file
                outfile.write('\t'.join(fields[:9] + new_samples) + '\n')

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Filter VCF file based on SQ, EQ, and NDQ values.')
    
    # Add arguments for input file, output file, and optional parameters
    parser.add_argument('-i', '--input', required=True, help='Path to the input VCF file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output VCF file.')
    parser.add_argument('-d', '--directory', help='Optional directory to create for the output file.')
    parser.add_argument('--min-sq', type=int, default=60, help='Minimum SQ value to keep a sample (default is 60).')
    parser.add_argument('--min-eq', type=int, default=60, help='Minimum EQ value to keep a sample (default is 60).')
    parser.add_argument('--min-ndq', type=int, default=60, help='Minimum NDQ value to keep a sample (default is 60).')
    
    # Parse arguments
    args = parser.parse_args()
    
    # If a directory is specified, create it if it does not exist
    if args.directory:
        # Import the os module for file and directory operations
        import os
        if not os.path.exists(args.directory):
            os.makedirs(args.directory)
            print(f"Created directory: {args.directory}")
    
    # Call the filter_vcf function with parsed arguments
    filter_vcf(args.input, args.output, args.min_sq, args.min_eq, args.min_ndq)

# Call the main function
if __name__ == "__main__":
    main()
