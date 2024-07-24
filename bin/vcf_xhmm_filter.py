#!/usr/bin/env python3

# Import argparse for parsing command-line arguments
import argparse

# Import os for file and directory operations
import os

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



def filter_vcf(input_vcf, output_vcf, min_sq, min_eq, min_ndq):
    """
    Filter the VCF file based on SQ, EQ, and NDQ criteria.

    Args:
        input_vcf (str): Path to the input VCF file.
        output_vcf (str): Path to the output filtered VCF file.
        min_sq (int): Minimum SQ value to keep.
        min_eq (int): Minimum EQ value to keep.
        min_ndq (int): Minimum NDQ value to keep.
    """
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                # Write header lines directly to the output file
                outfile.write(line)
            else:
                # Split line into fields
                parts = line.strip().split('\t')
                # Check if the line has at least 10 fields
                if len(parts) < 10:
                    continue

                # Extract the format field and sample fields
                format_field = parts[8]
                sample_fields = parts[9:]

                # Initialize a flag to check if any sample meets the criteria
                keep_line = False

                for sample_field in sample_fields:
                    # Parse the format and sample fields
                    format_dict = parse_sample_field(format_field, sample_field)

                    # Extract SQ, EQ, and NDQ values from the sample field
                    sq = format_dict.get('SQ', '0').split(',')[0]
                    eq = format_dict.get('EQ', '0').split(',')[0]
                    ndq = format_dict.get('NDQ', '0').split(',')[0]

                    try:
                        sq = float(sq)
                        eq = float(eq)
                        ndq = float(ndq)
                    except ValueError:
                        continue

                    # Check if the sample meets the criteria
                    if sq >= min_sq and eq >= min_eq and ndq >= min_ndq:
                        keep_line = True
                        break

                # Write the line to the output file if any sample meets the criteria
                if keep_line:
                    outfile.write(line)


def main():
    """
    Main function to parse command-line arguments and call the filter_vcf function
    """
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Filter VCF file based on SQ, EQ, and NDQ values.')

    # Add arguments for input file, output file, and optional parameters
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    parser.add_argument('-d', '--directory', help='Directory for output file (optional)')
    parser.add_argument('--min_sq', type=int, default=60, help='Minimum SQ value to keep (default: 60)')
    parser.add_argument('--min_eq', type=int, default=60, help='Minimum EQ value to keep (default: 60)')
    parser.add_argument('--min_ndq', type=int, default=60, help='Minimum NDQ value to keep (default: 60)')

    # Parse the command-line arguments and assign them to variable arg
    args = parser.parse_args()

    # Create directory if specified and does not exist
    if args.directory:
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created directory: {output_dir}")

    # Call the filter function
    filter_vcf(args.input, args.output, args.min_sq, args.min_eq, args.min_ndq)

# Call the main function
if __name__ == '__main__':
    main()
