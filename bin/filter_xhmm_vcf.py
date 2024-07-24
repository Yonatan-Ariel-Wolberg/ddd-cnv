#!/usr/bin/env python3

import argparse
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
    format_keys = format_field.split(':')
    sample_values = sample_field.split(':')
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
    try:
        with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    # Write header lines directly to the output file
                    outfile.write(line)
                else:
                    parts = line.strip().split('\t')
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
                            sq = int(sq)
                            eq = int(eq)
                            ndq = int(ndq)
                        except ValueError:
                            continue

                        # Check if the sample meets the criteria
                        if sq >= min_sq and eq >= min_eq and ndq >= min_ndq:
                            keep_line = True
                            break

                    # Write the line to the output file if any sample meets the criteria
                    if keep_line:
                        outfile.write(line)
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except IOError as e:
        print(f"IO Error: {e}")

def main():
    parser = argparse.ArgumentParser(description='Filter VCF file based on SQ, EQ, and NDQ values.')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    parser.add_argument('-d', '--directory', help='Directory for output file (optional)')
    parser.add_argument('--min_sq', type=int, default=60, help='Minimum SQ value to keep (default: 60)')
    parser.add_argument('--min_eq', type=int, default=60, help='Minimum EQ value to keep (default: 60)')
    parser.add_argument('--min_ndq', type=int, default=60, help='Minimum NDQ value to keep (default: 60)')

    args = parser.parse_args()

    # Construct the full output file path
    if args.directory:
        output_dir = args.directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created directory: {output_dir}")
        output_vcf = os.path.join(output_dir, args.output)
    else:
        output_vcf = args.output

    print(f"Filtering VCF file: {args.input}")
    print(f"Output file: {output_vcf}")

    # Call the filter function
    filter_vcf(args.input, output_vcf, args.min_sq, args.min_eq, args.min_ndq)

if __name__ == '__main__':
    main()

