#script for converting reformatted SNP data into multilocus genotype (MLG) file
#Author:Alex Samano, 2024
import sys

def process_file(input_file, output_file):
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                # strip newline and split into columns
                cols = line.strip().split('\t')
                
                # skip lines where the count of missing genotypes is greater than 2
                if sum(col == '.' for col in cols[4:]) > 2:
                    continue

                # process the genotypes
                processed_line = []
                for genotype in cols[4:]:
                    if genotype == '0/0':
                        processed_line.append(cols[2])  # reference allele
                    elif genotype == '1/1':
                        processed_line.append(cols[3])  # alternate allele
                    elif genotype in ('0/1', '1/0'):
                        processed_line.append('.')
                    else:
                        break  # skip line if any genotype doesn't match the criteria
                else:
                    # write the processed line to the output file
                    outfile.write('\t'.join(cols[:4] + processed_line) + '\n')

input_file = 'kochi_X23_onlysnps_filt.reformat'
output_file = 'kochi_X23_onlysnps.mlg'
process_file(input_file, output_file)
