#ifndef HMMSNP_H_
#define HMMSNP_H_
#include <vector>
#include "Chromosome.h"

#include "Individual.h"
#include <fstream>
#include <map>
#include <iostream>


class HMMSNP
{
public:

	virtual ~HMMSNP();

	HMMSNP(char* fileName);
	
	void baumWelch(Individual** individuals, char* outputFile);
	void initializeTransitionValues( Individual** individuals, double **** tempTransitionValues);
	void perChromosomeSet(Individual** individuals, int** maxStates, double**** tempTransitionValues,  float** iarray, int chromosomeIndex);
	void setEmissionProbabilities(  Individual** individuals, double **emissionProbabilities, vector <ExpressionInfo*> expressionSuperVector, int individualIndex, int chromosomeIndex,  float summationRatio);
	void forward(Individual** individuals, double** alphas, double **emissionProbabilities, double* scales, double*** transitionsLocal, int **viterbiStates, int size, int individualIndex, vector<ExpressionInfo*> expressionSuperVector);
	void backward(Individual **individuals, double** betas, double* scales, double*** transitionsLocal, double** emissionProbabilities, int size, int individualIndex, vector<ExpressionInfo*> expressionSuperVector);
	void viterbiTraceback(double** alphas, int** maxStates, int** viterbiStates, int size, int individualIndex);
	void fileOutput(Individual** individuals, ofstream &myViterbiFile,  ofstream &myfile, int chromosomeIndex, vector<ExpressionInfo*> expressionSuperVector, float** iarray, int** maxStates, int individualIndex, char* fileName);
	void addToTempTransitions (double**** tempTransitionValues,double** alphas, double**emissionProbabilities, double** betas, double***transitionsLocal, vector<ExpressionInfo*> expressionSuperVector, int chromosomeIndex);
	void updateTransitions(double**** tempTransitionValues, float weight, int size, int chromosomeIndex);
	void leftToRight (Individual ** individuals, char* outputFileName, double weight);
	void deleteLTOR(Individual **individuals);
	
	
 private:
	bool useLow;
	
	bool useSimplified;
	
	int numStates;
	
	static int * ilogsum_lookup;
	
	double runifValue;
	
	double ** transitionProbabilities;
	
	
		
	double * ratioMeans;
	
	double *ratioVariances;
	
	double **** ltorTransitions;

	double **oneStepTransitions;
	
	double *startProbs;
	
	double computeProbability(double mean, double variance, double value);

	
};


#endif /*HMMSNP_H_*/
