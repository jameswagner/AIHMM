#Input files
- One file for each chromosome (presently code only considers 22 human autosomal chromosomes). First line of each file contains the tag #SAMPLE followed
sample ids, then one line for each SNP, containing the SNP id, chromosome and location, followed by three numbers for each individual
(can be tab or space separated). For each individual, in each line, the first number is 0 if the SNP is homozygous in that
individual and 1 if heterozygous, the second number indicates the total intensity level of the two alleles for that SNP (sum of the
two cDNA intensities divided by the sum of the two gDNA intensities, then take log base 2). The third number indicates the
allelic imbalance (or allelic expression) "ratio of ratios". That is the ratio of the first allele intensity to the second allele intensity,
normalized with respect to gDNA. Then the log base 2 of this ratio is taken, so that 0 would indicate perfectly equal intensities, +1 would 
indicate that the first allele has double the intensity of the second allele, etc. Note that the definition of
"first" and "second" allele is somewhat arbitrary but does depend on haplotype phasing. We expect if consecutive SNPs
have the same directionality of allelic imbalance, and there are no phasing errors, then they will all have the same sign for the allelic
imbalance ratio in this region in a given individual.

Input files are given for CEU LCL's as discussed in PLoS 2010 Wagner et. al. paper, in the directory data.


Also required is a file describing the parameters for an HMM. A suitable file that can be used as initial parameters (still require Baum Welch optimization)
would be ceu8StateHMMParams.txt. 
The first line indicates the number of states in the model (an even number >= 4).
Then one line for each state, the first number indicates the start probability (probability the first SNP of a chromosome in any particular individual will
start in that state), then the mean, variance, normal distribution, then transition probabilities to all states in the model. 
Note that in this example, there are 8 states, the 8th state represents the state of low/no expression at a given SNP (very low cDNA intensities). The first
7 states indicate varying degrees of allelic imbalance, ranging from very negative to very positive, with the middle state indicating that at a given 
SNP there is expression but no allelic imbalance. 


# Code structure

Code is written in C++, a number of .cpp and .h files exist for data structures, to represent an Individual, a Chromosome, a SNP and (allelic) ExpressionInfo 
for a particular SNP. FileReader reads in the files of the format specified above and stores the data in C++ data structures such that HMM's can be applied to
 them. HMMSNP contains the code for Baum Welch learning of parameters, and the "left to right" HMM described in the paper, including output of HMM-smoothed allelic imbalance values. allelicImbalance.cpp contains the main function and will call the necessary functions in other classes to read in the input files, train an HMM and 
output results.

# Building and running code (note that code was only implemented and tested on Ubuntu systems. It is expected to work on other Linux-based systems or Mac 
but has in no way been tested or designed for Windows).

cd to allelic/Debug
make clean  
make all 
now from the same directory can run the command ./AIHMM  -numChromosomes 22 -numIndividuals 55 -startFile ceu8StateHMMParams.txt -fileNamePrefix  ../../data/ceuRevisedRatiosChrom  -outputFile outputFile.txt 

numChromosomes is typically 22 but can be a smaller number if for testing purposes you want to just run a smaller subset of the genome and have results faster
numIndividuals is, as the name suggests, the number of samples for which there is AI data
the startFile is the HMM parameter file
fileNamePrefix is the location of the input AI data, not that it is only a prefix and a number (between 1 and 22) will be appended
for the chromosome number followed by the suffix ".txt", eg. ceuRevisedRatiosChrom15.txt for the allelic imbalance data for SNPS located on chromosome 15.
the outputFile is a path indicating the names of files where the final Baum-Welch optimized model file will be stored as well as information about the HMM output.

For a full genome run with the LCL data at its present size the HMM training and then the left to right training and testing are expected to take several hours each on a typical processor. Unfortunately
it requires quite a bit of RAM in its present implementation.


# outputFile format

-the output file specified by the outputFile parameter will contain the final trained model parameters, in the same format as the input model file.
(Note that this outputFile can be used for future runs as the startFile, and the line in allelicImbalance.cpp in which baumWelch is called
can be commented out)

- Suffixes will be added to this output file name for bed files containing information about the allelic expression after having been smoothed by the HMM. 
With the current implementation a total of 4 such files will be formed  for each format (it is a good idea to specify another destination directory as these files can really pile up!) The files containing "iter0" represent information at the start of left to right HMM 
(i.e. just using a typical "Ergodic" HMM and not yet with any left-to-right specific training of parameters). What is probably more of interest would be
the "iter99" files that represent the information at the end of 100 iterations of the left to right HMM. The two files for "iter99" have suffixes
"viterbi" and "evals". The viterbi file represents the output after running the left to right HMM and then running the Viterbi algorithm to determine
the most likely state in each SNP, for each sample. It is in bedgraph format and for each individual, there will be a track line, as well as a line for each SNP in 
that individual. For each particular line there is a value chrXX to indicate the chromosome it is on, the location of the SNP, the location of the SNP plus 1, then
the viterbi state that SNP is in for that individual. The states are numbered 0-7, with 7 being non-expression, 3 indicates expression but no allelic expression
0 and 6 are the highest magnitude of allelic expression, 1 and 5 are moderate AE and 2 and 4 are weak AI. 

The files with the ".evals" suffix represent a sort of smoothed value for allelic expression, after having run the HMM which will take into account AI present
at neighboring SNPs and in other individuals. They are the same format except the final column of lines for each SNP represents the smoothed
e-value, which is obtained by summing for each state from 0-7 the probability that SNP is in that SNP multiplied by the mean allelic expression level at that state.
Note that no distinction is made in this type file of between SNPs that have expression but no allelic expression, and those SNPs which are not expressed at all.
Note that these e-values represent the LOG of (cDNA allelic expression / gDNA allelic expression), so as with the input files a value of
0 represents no allelic expression, and highly positive or negative values indicate that one or the other allele is expressed more than the other. 

#outputFiles

- precomputed output files are available in the allelic/Debug/outputFiles directory.
