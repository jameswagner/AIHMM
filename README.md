# AIHMM: Allelic Imbalance Hidden Markov Model

A tool for detecting and analyzing allelic imbalance using Hidden Markov Models.

## Background

This tool was originally developed to analyze allelic imbalance in gene expression as discussed in PLoS 2010 Wagner et al. paper. The implementation uses Hidden Markov Models to detect and smooth allelic imbalance signals across chromosomes.

In this context:
- Allelic imbalance refers to unequal expression of alleles from the maternal and paternal chromosomes
- The HMM helps identify regions of consistent allelic imbalance by considering neighboring SNPs and data from multiple individuals
- The model includes states for varying degrees of allelic imbalance, as well as a state for SNPs with low/no expression

## Project Structure

```
AIHMM/
├── src/                     # Source code
│   ├── core/                # Core data structures
│   ├── hmm/                 # HMM implementation
│   ├── io/                  # Input/output utilities
│   └── main.cpp             # Main application
├── models/                  # HMM parameter files
├── data/                    # Input data files
├── output/                  # Output directory
├── scripts/                 # Utility scripts
├── build/                   # Build directory
└── docs/                    # Documentation
```

## Requirements

- C++ compiler with C++14 support (GCC 5+ or equivalent)
- Boost C++ Libraries (specifically program_options)
- CMake 3.10 or higher

## Building the Project

### Using CMake (Recommended)

```bash
# Create a build directory
mkdir -p build
cd build

# Configure and build
cmake ..
make

# Or use the provided script
./build.sh
```

## Running AIHMM

After building, you can run AIHMM with the following options:

```bash
# From build directory
./bin/AIHMM --numChromosomes 22 --numIndividuals 55 --startFile ceu8StateHMMParams.txt --fileNamePrefix ../data/ceuRevisedRatiosChrom --outputFile ../output/outputFile.txt
```

### Command Line Options

- `--numChromosomes`: Number of chromosomes to process (typically 22)
- `--numIndividuals`: Number of samples with allelic imbalance data
- `--startFile`: HMM parameter file
- `--fileNamePrefix`: Prefix for input allelic imbalance data files
- `--outputFile`: Path for output files

## Input File Format

### Allelic Imbalance Data Files

One file per chromosome with the following format:
- First line: `#SAMPLE` followed by sample IDs
- Subsequent lines: SNP ID, chromosome, location, followed by three numbers for each individual:
  1. Heterozygosity (0=homozygous, 1=heterozygous)
  2. Total intensity level (log2 of cDNA/gDNA ratios)
  3. Allelic imbalance ratio (log2 of normalized ratio between alleles)

For each individual, in each line, the first number is 0 if the SNP is homozygous in that individual and 1 if heterozygous. The second number indicates the total intensity level of the two alleles for that SNP (sum of the two cDNA intensities divided by the sum of the two gDNA intensities, then take log base 2). The third number indicates the allelic imbalance (or allelic expression) "ratio of ratios". That is the ratio of the first allele intensity to the second allele intensity, normalized with respect to gDNA. Then the log base 2 of this ratio is taken, so that 0 would indicate perfectly equal intensities, +1 would indicate that the first allele has double the intensity of the second allele, etc.

Note that the definition of "first" and "second" allele is somewhat arbitrary but does depend on haplotype phasing. We expect if consecutive SNPs have the same directionality of allelic imbalance, and there are no phasing errors, then they will all have the same sign for the allelic imbalance ratio in this region in a given individual.

### HMM Parameter File

- First line: Number of states in the model (an even number ≥ 4)
- Subsequent lines: One per state with start probability, mean, variance, and transition probabilities

The first number indicates the start probability (probability the first SNP of a chromosome in any particular individual will start in that state), then the mean, variance, normal distribution, then transition probabilities to all states in the model.

In a typical 8-state model, the 8th state represents the state of low/no expression at a given SNP (very low cDNA intensities). The first 7 states indicate varying degrees of allelic imbalance, ranging from very negative to very positive, with the middle state indicating that at a given SNP there is expression but no allelic imbalance.

## Output Files

- Main output file: Trained model parameters, in the same format as the input model file
- Files with suffix `.viterbi`: Most likely state for each SNP (bedgraph format)
- Files with suffix `.evals`: Smoothed allelic expression values (bedgraph format)

The viterbi file represents the output after running the left-to-right HMM and then running the Viterbi algorithm to determine the most likely state in each SNP, for each sample. It is in bedgraph format and for each individual, there will be a track line, as well as a line for each SNP in that individual. For each particular line there is a value chrXX to indicate the chromosome it is on, the location of the SNP, the location of the SNP plus 1, then the viterbi state that SNP is in for that individual. The states are numbered 0-7, with 7 being non-expression, 3 indicates expression but no allelic expression, 0 and 6 are the highest magnitude of allelic expression, 1 and 5 are moderate AE, and 2 and 4 are weak AI.

The files with the `.evals` suffix represent a sort of smoothed value for allelic expression, after having run the HMM which will take into account AI present at neighboring SNPs and in other individuals. They are the same format except the final column of lines for each SNP represents the smoothed e-value, which is obtained by summing for each state from 0-7 the probability that SNP is in that state multiplied by the mean allelic expression level at that state. No distinction is made in this type of file between SNPs that have expression but no allelic expression, and those SNPs which are not expressed at all. These e-values represent the LOG of (cDNA allelic expression / gDNA allelic expression), so as with the input files a value of 0 represents no allelic expression, and highly positive or negative values indicate that one or the other allele is expressed more than the other.

For more detailed information, see the documentation in the `docs` directory. 