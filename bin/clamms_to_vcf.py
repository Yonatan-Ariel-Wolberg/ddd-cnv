import os

def convert_clamms_to_vcf(bed_file, vcf_file):
    with open(bed_file, 'r') as infile, open(vcf_file, 'w') as outfile:
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in infile:
            fields = line.strip().split('\t')
            chrom, start, end, interval, sample, type, cn, num_window, q_some, q_exact, q_left_extend, left_extend_coord, q_right_extend, right_extend_coord, q_left_contract, left_contract_coord, q_right_contract, right_contract_coord = fields
            start = int(start)
            end = int(end)
            id_ = fields[4]
            ref = 'N'
            alt = '<' + type + '>'
            qual = q_some
            filter = 'PASS' if q_exact >= 0 else 'LowQuality'
            info = f'SVTYPE=CNV;CN={cn};SAMPLE={sample};NUM_WINDOW={num_window};Q_SOME={q_some};Q_EXACT={q_exact};Q_LEFT_EXTEND={q_left_extend};LEFT_EXTEND_COORD={left_extend_coord};Q_RIGHT_EXTEND={q_right_extend};RIGHT_EXTEND_COORD={right_extend_coord};Q_LEFT_CONTRACT={q_left_contract};LEFT_CONTRACT_COORD={left_contract_coord};Q_RIGHT_CONTRACT={q_right_contract};RIGHT_CONTRACT_COORD={right_contract_coord}'
            outfile.write(f'{chrom}\t{start}\t{id_}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n')

convert_clamms_to_vcf(bed_file, vcf_file)
