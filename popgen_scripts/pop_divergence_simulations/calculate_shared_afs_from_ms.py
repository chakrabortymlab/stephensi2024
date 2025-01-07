import numpy as np
import argparse
from scipy.stats import pearsonr

def calculate_mutation_frequencies_and_r2(file_path):
    """
    calculate allele frequencies in two populations from an msms output file
    and compute the correlation coefficient between their frequencies
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # parse  file to find the //
    samples_start = False
    samples = []
    for line in lines:
        line = line.strip()
        if line.startswith("//"):  # skip lines until the //
            samples_start = True
            continue
        if not samples_start or line.startswith("segsites:") or line.startswith("positions:") or not line:
            continue
        #collect sampled chromosomes
        samples.append(line)

    
    # convert to numpy
        matrix = np.array([list(sample) for sample in samples], dtype=int)

    # divide into two populations 
    pop1 = matrix[:36, :]  # first 36 chromosomes are coastal (kochi or tvm)
    pop2 = matrix[36:156, :]  # next 120 chromosomes are island

    # calculate allele frequencies in each pop
    pop1_frequencies = pop1.mean(axis=0)
    pop2_frequencies = pop2.mean(axis=0)

    # calculate R2 
    r, _ = pearsonr(pop1_frequencies, pop2_frequencies)
    r_squared = r ** 2

    # put the result in a dictionry
    result = {
        "Pop1 ": pop1_frequencies.tolist(),
        "Pop2 ": pop2_frequencies.tolist(),
        "R2": r_squared
    }

    return result

def main():
    # set up argument parsing
    parser = argparse.ArgumentParser(description="calculate allele frequencies in two populations from an msms output file")
    parser.add_argument("file", type=str, help="name of ms output file.")
    
    args = parser.parse_args()
    file_path = args.file

    # calculate frequencies and R2, throw error if input has weird header
    try:
        result = calculate_mutation_frequencies_and_r2(file_path)
        print("R2:", result["R2"])
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
