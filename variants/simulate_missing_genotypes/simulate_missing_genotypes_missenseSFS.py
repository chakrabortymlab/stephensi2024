import random
import numpy as np
import argparse

# script to simulate missing data, calculate site frequency spectra of missense variants, and evaluate the effect of error rate on SFS 

# randomly sprinkles missing genotypes at a given error rate
# data converted from vcf to simple text file: position, annotation, Sample1_genotype,Sample2_genotype, ...

def add_missing_genotypes(data, error_rate):
    for row in data:
        for i in range(2, len(row)):  
            if random.random() < error_rate:
                row[i] = '.'
    return data

# function to filter data based on 'missense_variant' annotation and the minimum number of nonmissing haplotypes
def filter_data(data, min_non_missing=90):
    filtered_data = []
    for row in data:
        # filter for missense variants
        if row[1] == 'missense_variant':
            non_missing = sum(1 for x in row[2:] if x != '.')
            # keep the row only if the number of non-missing haplotypes is at least min_non_missing
            if non_missing >= min_non_missing:
                filtered_data.append(row)
    return filtered_data

# function to compute the folded site frequency spectrum 
def compute_sfs(data):
    sfs = []
    for row in data:
        genotypes = row[2:]  # assuming genotypes start at column 10
        allele_counts = {'0': 0, '1': 0}  # count of each genotypes

        for genotype in genotypes:
            if genotype != '.':
                alleles = genotype.split('/')
                for allele in alleles:
                    if allele in allele_counts:
                        allele_counts[allele] += 1
        
        # compute the site frequency 
        non_ref_count = allele_counts['1']
        sfs.append(non_ref_count)
    
    # fold the SFS 
    max_freq = max(sfs)
    folded_sfs = [0] * (max_freq + 1)
    for freq in sfs:
        folded_sfs[min(freq, max_freq - freq)] += 1

    return folded_sfs


#  handle argument parsing
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str)
    parser.add_argument('output_file', type=str)
    return parser.parse_args()

# function to read data from a file
def read_data(input_file):
    data = []
    with open(input_file, 'r') as file:
        for line in file:
            data.append(line.strip().split('\t'))
    return data

# main function to run the process
def main():
    # parse arguments
    args = parse_args()

    # read data from the input file
    data = read_data(args.input_file)

    # Parameters
    error_rate = 0.1  # adjust to 10,20,30 ...
    min_non_missing = 90  # min number of nonmissing haplotypes

    #apply random missing genotypes
    data_with_missing = add_missing_genotypes(data, error_rate)

    # rilter data to only keep missense_variant and rows with sufficient nonmissing genotypes
    filtered_data = filter_data(data_with_missing, min_non_missing)

    # compute folded site frequency spectrum
    sfs = compute_sfs(filtered_data)

    #output the folded site frequency spectrum to the specified file
    with open(args.output_file, 'w') as output_file:
        output_file.write("Folded SFS:\n")
        output_file.write(" ".join(map(str, sfs)) + "\n")

# run 
if __name__ == "__main__":
    main()
