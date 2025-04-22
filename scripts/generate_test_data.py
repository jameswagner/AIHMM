import random
import numpy as np

def generate_test_data(expression_levels_file):

    num_individuals = None
    expression_means = None
    imbalance_means = None
    num_snps = None

    # Read expression levels and imbalance ratios from the provided file
    with open(expression_levels_file, 'r') as file:
        num_individuals = int(file.readline().strip())
        expression_means = list(map(int, file.readline().strip().split(',')))
        imbalance_means = list(map(float, file.readline().strip().split(',')))
        num_snps = len(expression_means)

    # Initialize data structure to store test data
    test_data = []

    # Generate test data for each SNP
    for snp_index in range(num_snps):
        snp_data = [ f'rs{snp_index + 1}', 1,  (snp_index + 1) * 100] #snp id, chromosome, location

        for _ in range(num_individuals):
            individual_data = []

            # Generate heterozygosity status (1 for heterozygous, 0 for homozygous)
            heterozygous = random.random() < 0.25
            individual_data.append(1 if heterozygous else 0)

            # Randomly sample expression level
            expression_level = np.random.normal(expression_means[snp_index], 0.25)
            individual_data.append(round(expression_level, 4))

            # If individual is heterozygous, sample imbalance ratio, else set to 0
            if heterozygous:
                imbalance_ratio = np.random.normal(imbalance_means[snp_index], 0.25)
                individual_data.append(round(imbalance_ratio, 4))
            else:
                individual_data.append(0)

            snp_data.extend(individual_data)

        test_data.append(snp_data)
    print(num_individuals)
    print(test_data)
    return num_individuals, test_data

# Example usage
expression_levels_file = 'expression_levels.txt'
num_individuals, test_data = generate_test_data(expression_levels_file)

for snp_data in test_data:
    print('\t'.join(map(str, snp_data)))
