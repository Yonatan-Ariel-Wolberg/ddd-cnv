#!/usr/bin/python3

# Import argparse for parsing command-line arguments
import argparse

# Import os for file and directory operations
import os

def filter_vcf(input_vcf, output_vcf):
    """
    Filters VCF file based on SQ, EQ and NDQ values in sample fields
    """
    
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            # Write header lines to output file
            if line.startswith('#'):
                outfile.write(line)
            else:
                # Split the line into fields
                parts = line.strip().split('\t')
                
                # Skip lines with less than 10 fields
                if len(parts) < 10:
                    continue
                
                # Extract fields
                info_field = parts[7]
                format_field = parts[8]
                sample_fields = parts[9:]
                
                # Extract SQ, EQ, and NDQ
                # Split info field into a dictionary where the item before the '=' is the key
                info_dict = {item.split('=')[0]: item.split('=')[1] for item in info_field.split(';') if '=' in item}
                
                # Split format field into a dictionary
                
                format_dict = {item.split(':')[0]: item.split(':')[1:] for item in format_field.split(':')}
                
                filter_passed = True
                for sample in sample_fields:
                    sample_data = sample.split(':')
                    if len(sample_data) != len(format_dict):
                        filter_passed = False
                        break
                    
                    format_values = {fmt: val for fmt, val in zip(format_dict.keys(), sample_data)}
                    
                    SQ = int(format_values.get('SQ', 0))
                    EQ = int(format_values.get('EQ', 0))
                    NDQ = int(format_values.get('NDQ', 0))
                    
                    if SQ < 60 or EQ < 60 or NDQ < 60:
                        filter_passed = False
                        break
                
                if filter_passed:
                    outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Filter VCF file based on SQ, EQ, and NDQ values.")
    parser.add_argument('-i', '--input', required=True, help='Input VCF file.')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file.')
    parser.add_argument('-d', '--directory', help='Optional directory to create for the output file.')

    args = parser.parse_args()

    if args.directory:
        # If directory is specified, create it if it does not exist
        if not os.path.exists(args.directory):
            os.makedirs(args.directory)
        output_vcf = os.path.join(args.directory, os.path.basename(args.output))
    else:
        # If no directory is specified, use the output file path as is
        output_vcf = args.output

    # Ensure the directory of the output file exists
    output_dir = os.path.dirname(output_vcf)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filter_vcf(args.input, output_vcf)

if __name__ == '__main__':
    main()