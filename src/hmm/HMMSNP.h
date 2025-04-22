#ifndef HMMSNP_H_
#define HMMSNP_H_

#include <vector>
#include "Chromosome.h"
#include "Individual.h"
#include <fstream>
#include <map>
#include <iostream>
#include <memory>
#include "constants.h"

class HMMSNP {
public:
    HMMSNP(const std::string& fileName);
    virtual ~HMMSNP();
    void baumWelch(std::vector<std::unique_ptr<Individual>>& individuals, const std::string& outputFile);
    void leftToRight(std::vector<std::unique_ptr<Individual>>& individuals, const std::string outputFileName, double weight);

private:
    bool useLow;
    bool useSimplified;
    int numStates;
    double runifValue;
    std::vector<std::vector<double>> transitionProbabilities;
    std::vector<double> ratioMeans;
    std::vector<double> ratioVariances;
    std::vector<std::vector<std::vector<std::vector<double>>>> ltorTransitions;
    std::vector<std::vector<double>> oneStepTransitions;
    std::vector<double> startProbs;


    static std::vector<int> ilogsum_lookup;
    
    void initializeTransitionValues(std::vector<std::unique_ptr<Individual>>& individuals, std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues);
    void perChromosomeSet(std::vector<std::unique_ptr<Individual>>& individuals, std::vector<std::vector<int>>& maxStates, std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues, std::vector<std::vector<float>>& iarray, int chromosomeIndex);
    void setEmissionProbabilities(std::vector<std::unique_ptr<Individual>>& individuals, std::vector<std::vector<double>>& emissionProbabilities, 
    std::vector<ExpressionInfo*>& expressionSuperVector, int individualIndex, int chromosomeIndex, float summationRatio);
    void forward(std::vector<std::vector<double>>& alphas, std::vector<std::vector<double>>& emissionProbabilities, std::vector<double>& scales, 
    std::vector<std::vector<std::vector<double>>>& transitionsLocal, std::vector<std::vector<int>>& viterbiStates);
    void backward(std::vector<std::unique_ptr<Individual>>& individuals, std::vector<std::vector<double>>& betas, std::vector<double>& scales, std::vector<std::vector<std::vector<double>>>& transitionsLocal, std::vector<std::vector<double>>& emissionProbabilities, int size, int individualIndex, std::vector<ExpressionInfo*>& expressionSuperVector);
    void viterbiTraceback(std::vector<std::vector<double>>& alphas, std::vector<std::vector<int>>& maxStates, std::vector<std::vector<int>>& viterbiStates, int individualIndex);
    void fileOutput(std::vector<std::unique_ptr<Individual>>& individuals, std::ofstream& myViterbiFile, std::ofstream& myfile, int chromosomeIndex, std::vector<ExpressionInfo*>& expressionSuperVector, std::vector<std::vector<float>>& iarray, std::vector<std::vector<int>>& maxStates, int individualIndex, const std::string& fileName);
    void addToTempTransitions(std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues, std::vector<std::vector<double>>& alphas, std::vector<std::vector<double>>& emissionProbabilities, std::vector<std::vector<double>>& betas, std::vector<std::vector<std::vector<double>>>& transitionsLocal, std::vector<ExpressionInfo*>& expressionSuperVector, int chromosomeIndex);
    void updateTransitions(std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues, float weight, int size, int chromosomeIndex);
    void initializeArrays(int numStates, std::vector<std::vector<double>>& alphas, std::vector<std::vector<double>>& betas, std::vector<std::vector<double>>& emissionProbabilities);
    void computeEmissionProbabilities(std::vector<ExpressionInfo*>& expressionSuperVector, int numStates, float* expressedTableLocal, float* unexpressedTableLocal, std::vector<double>& ratioMeans, std::vector<double>& ratioVariances, double runifValue, std::vector<std::vector<double>>& emissionProbabilities);


    double computeProbability(double mean, double variance, double value);
};

#endif /* HMMSNP_H_ */
