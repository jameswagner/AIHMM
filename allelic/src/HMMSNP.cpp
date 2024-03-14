#include "HMMSNP.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ExpressionInfo.h"
#include "constants.h"
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <map>
#include <sstream>


const double PI =  3.14159265;


HMMSNP::~HMMSNP() {

}

HMMSNP::HMMSNP(const std::string& fileName) {

    useLow = true;
    runifValue = 0.01;

    std::ifstream infile(fileName);
    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    std::string line;
    std::getline(infile, line);
    numStates = std::stoi(line);

    ratioMeans.resize(numStates);
    ratioVariances.resize(numStates);
    transitionProbabilities.resize(numStates);
    for (int i = 0; i < numStates; i++) {
        transitionProbabilities[i].resize(numStates);
    }
    startProbs.resize(numStates);

    for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        std::getline(infile, line);
        std::istringstream iss(line);

        iss >> startProbs[stateIndex];
        iss >> ratioMeans[stateIndex];
        iss >> ratioVariances[stateIndex];

        for (int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
            iss >> transitionProbabilities[stateIndex][innerStateIndex];
        }
    }
    infile.close();
}

double HMMSNP::computeProbability(double mean, double variance, double value) {
    double std = sqrt(variance);
    return (1 / (std * sqrt(2. * PI)) * exp(-1 * (value - mean) * (value - mean) / (2 * variance)));
}

 

void HMMSNP::initializeArrays(int numStates, std::vector<std::vector<double>>& alphas, std::vector<std::vector<double>>& betas, std::vector<std::vector<double>>& emissionProbabilities) {
    // Resize vectors
    alphas.resize(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    betas.resize(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    emissionProbabilities.resize(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
}

void HMMSNP::computeEmissionProbabilities(std::vector<ExpressionInfo*>& expressionSuperVector, int numStates, float* expressedTableLocal, float* unexpressedTableLocal, std::vector<double>& ratioMeans, std::vector<double>& ratioVariances, double runifValue, std::vector<std::vector<double>>& emissionProbabilities) {
    for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
        int logRatioIndex = (int)((expressionSuperVector[expressionIndex]->getExpressionLogRatio() + 10) *10);
        if(logRatioIndex < 0) {
            logRatioIndex = 0;
        }
        if(logRatioIndex > 139) {
            logRatioIndex = 139;
        }
        for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
            emissionProbabilities[expressionIndex][emitIndex] = 1.0;

            if(emitIndex < numStates-1 && (numStates %2 == 0)) {
                emissionProbabilities[expressionIndex][emitIndex]  *= expressedTableLocal[logRatioIndex];
            }
            else if(numStates %2 == 0) {
                emissionProbabilities[expressionIndex][emitIndex]  *= unexpressedTableLocal[logRatioIndex];
            }

            if(expressionSuperVector[expressionIndex]->isHeterozygote()) {
                emissionProbabilities[expressionIndex][emitIndex] *= computeProbability(ratioMeans[emitIndex], ratioVariances[emitIndex], expressionSuperVector[expressionIndex]->getRatioOfRatios());
            }

            emissionProbabilities[expressionIndex][emitIndex] *= 10;
            emissionProbabilities[expressionIndex][emitIndex] += runifValue;
        }
    }
}



void forwardPart( int numStates, const std::vector<double>& startProbs, const std::vector<std::vector<double>>& transitionProbabilities, 
const std::vector<std::vector<double>>& emissionProbabilities, std::vector<std::vector<double>>& alphas, std::vector<double>& scales) {
    double sumAlphas = 0.0;
    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        alphas[0][stateIndex] = startProbs[stateIndex] * emissionProbabilities[0][stateIndex];
        sumAlphas += alphas[0][stateIndex];
    }
    for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
        alphas[0][emitIndex] /= sumAlphas;
    }
    scales[0] = sumAlphas;
    for(unsigned int alphaIndex = 1; alphaIndex < emissionProbabilities.size(); alphaIndex++) {
        sumAlphas = 0.0;
        for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
            double alphaVal = 0.0;
            for(int innerEmitIndex = 0; innerEmitIndex < numStates; innerEmitIndex++) {
                alphaVal += transitionProbabilities[innerEmitIndex][emitIndex] * alphas[alphaIndex-1][innerEmitIndex];
            }
            alphaVal *= emissionProbabilities[alphaIndex][emitIndex];
            alphas[alphaIndex][emitIndex] = alphaVal;
            sumAlphas += alphaVal;
        }
        scales[alphaIndex] = sumAlphas;
        for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
            alphas[alphaIndex][emitIndex] /= sumAlphas;
        }
    }
}



void HMMSNP::forward(std::vector<std::vector<double>>& alphas, std::vector<std::vector<double>>& emissionProbabilities, std::vector<double>& scales, 
std::vector<std::vector<std::vector<double>>>& transitionsLocal, std::vector<std::vector<int>>& viterbiStates) {
    double sumAlphas = 0.0;
    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        alphas[0][stateIndex] = startProbs[stateIndex] * emissionProbabilities[0][stateIndex];
        sumAlphas += alphas[0][stateIndex];
    }
    scales[0] = sumAlphas;
    for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
        alphas[0][emitIndex] /= sumAlphas;
    }
    for(unsigned int alphaIndex = 1; alphaIndex < emissionProbabilities.size(); alphaIndex++) {
        int previousAlphaIndex = alphaIndex - 1;

        sumAlphas = 0.0;
        for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
            double alphaVal = 0.0;
            float maxAlpha = 0;
            for(int innerEmitIndex = 0; innerEmitIndex < numStates; innerEmitIndex++) {
                alphaVal += transitionsLocal[previousAlphaIndex][innerEmitIndex][emitIndex] * alphas[previousAlphaIndex][innerEmitIndex];
            
                if(transitionsLocal[previousAlphaIndex][innerEmitIndex][emitIndex] * alphas[previousAlphaIndex][innerEmitIndex]  > maxAlpha) {
                    maxAlpha = transitionsLocal[previousAlphaIndex][innerEmitIndex][emitIndex] * alphas[previousAlphaIndex][innerEmitIndex];
                    viterbiStates[alphaIndex][emitIndex] = innerEmitIndex;
                }
            }
            alphaVal *= emissionProbabilities[alphaIndex][emitIndex];
            alphas[alphaIndex][emitIndex] = alphaVal;
            sumAlphas += alphaVal;
        }
        scales[alphaIndex] = sumAlphas;
        for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
            alphas[alphaIndex][emitIndex] /= sumAlphas;
            if(isnan(alphas[alphaIndex][emitIndex])) {
                printf("This is a problem in alpha %d %d from emission %f previous alpha %f\n", alphaIndex, emitIndex,  emissionProbabilities[alphaIndex][emitIndex],  alphas[alphaIndex-1][emitIndex] );
            }
        }
    } 
}




// Backward Part
void backwardPart(const std::vector<ExpressionInfo*>& expressionSuperVector, int numStates, const std::vector<std::vector<double>>& transitionProbabilities, 
const std::vector<std::vector<double>>& emissionProbabilities, std::vector<std::vector<double>>& betas, const std::vector<double>& scales) {
    
    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        betas[expressionSuperVector.size()-1][stateIndex] = 1.0/numStates;
    }
    for (int betaIndex = expressionSuperVector.size() - 2; betaIndex >= 0; betaIndex--) {
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            betas[betaIndex][stateIndex] = 0.0;
            for (int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
                betas[betaIndex][stateIndex] += transitionProbabilities[stateIndex][innerStateIndex] * betas[betaIndex + 1][innerStateIndex] * emissionProbabilities[betaIndex + 1][innerStateIndex];
            }
        }
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            betas[betaIndex][stateIndex] /= scales[betaIndex];
        }
    }
}

// E Step
void eStep(const std::vector<ExpressionInfo*>& expressionSuperVector, int numStates, const std::vector<std::vector<double>>& alphas, 
const std::vector<std::vector<double>>& betas, std::vector<std::vector<double>>& iarray, std::vector<double>& sumWeightHets, 
double& iarray0Sum, const std::vector<std::vector<double>>& transitionProbabilities, const std::vector<std::vector<double>>& emissionProbabilities, 
std::vector<std::vector<double>>& transitionSums
) {
    for (unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
        double sumIndex = 0.0;
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            iarray[expressionIndex][stateIndex] = alphas[expressionIndex][stateIndex] * betas[expressionIndex][stateIndex];
            sumIndex += iarray[expressionIndex][stateIndex];
        }
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            iarray[expressionIndex][stateIndex] /= sumIndex;
        }
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            if (expressionSuperVector[expressionIndex]->isHeterozygote() && iarray[expressionIndex][stateIndex] > 0.05) {
                sumWeightHets[stateIndex] += iarray[expressionIndex][stateIndex];
            }
        }
    }
    
    iarray0Sum = 0.0;
    for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        iarray0Sum += iarray[0][stateIndex];
    }

    for (int firstState = 0; firstState < numStates; firstState++) {
        for (int secondState = 0; secondState < numStates; secondState++) {
            transitionSums[firstState][secondState] = transitionProbabilities[firstState][secondState] * 1e-05;
        }
    }
    
    for (unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size() - 1; expressionIndex++) {
        for (int firstState = 0; firstState < numStates; firstState++) {
            if (iarray[expressionIndex][firstState] > 0.05 || firstState == numStates - 1) {
                for (int secondState = 0; secondState < numStates; secondState++) {
                    transitionSums[firstState][secondState] += alphas[expressionIndex][firstState] * transitionProbabilities[firstState][secondState] * emissionProbabilities[expressionIndex + 1][secondState] * betas[expressionIndex + 1][secondState];
                }
            }
        }
    }
}


void mStep(std::vector<ExpressionInfo*>& expressionSuperVector, int numStates, std::vector<std::vector<double>>& iarray,
std::vector<double>& sumWeightHets, std::vector<double>& sumStartProbs, std::vector<double>& sumRatioMeans,
std::vector<double>& sumRatioVariances, double iarray0Sum, double sumWeightAllSNPs, std::vector<std::vector<double>>& transitionSums,
std::vector<std::vector<double>>& sumTransitionProbabilities) {

    std::vector<double> transitionFullSums(numStates);

    for (int firstState = 0; firstState < numStates; firstState++) {
        for (int secondState = 0; secondState < numStates; secondState++) {
            transitionFullSums[firstState] += transitionSums[firstState][secondState];
        }
    }

    for (int firstState = 0; firstState < numStates; firstState++) {
        for (int secondState = 0; secondState < numStates; secondState++) {
            sumTransitionProbabilities[firstState][secondState] += transitionSums[firstState][secondState] / transitionFullSums[firstState] * expressionSuperVector.size();
        }
    }

    std::vector<double> sumOfMeans(numStates);
    std::vector<double> sumOfVariances(numStates);

    for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        sumOfMeans[stateIndex] = 0.0;
        sumOfVariances[stateIndex] = 0.001;

        for (unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
            if (expressionSuperVector[expressionIndex]->isHeterozygote() && iarray[expressionIndex][stateIndex] > 0.05) {
                sumOfMeans[stateIndex] += expressionSuperVector[expressionIndex]->getRatioOfRatios() * iarray[expressionIndex][stateIndex] / sumWeightHets[stateIndex];
            }
            else if (expressionSuperVector[expressionIndex]->isHeterozygote() && stateIndex == numStates - 1 && iarray[expressionIndex][stateIndex] > 0.0) {
                sumOfMeans[stateIndex] += expressionSuperVector[expressionIndex]->getRatioOfRatios() * iarray[expressionIndex][stateIndex] / sumWeightHets[stateIndex];
            }
        }
    }

    for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        for (unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
            if (expressionSuperVector[expressionIndex]->isHeterozygote() && iarray[expressionIndex][stateIndex] > 0.05) {
                sumOfVariances[stateIndex] += (expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex]) * (expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex]) * iarray[expressionIndex][stateIndex] / sumWeightHets[stateIndex];
            }

            if (expressionSuperVector[expressionIndex]->isHeterozygote() && stateIndex == numStates - 1 && iarray[expressionIndex][stateIndex] > 0.0) {
                sumOfVariances[stateIndex] += (expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex]) * (expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex]) * iarray[expressionIndex][stateIndex] / sumWeightHets[stateIndex];
            }
        }
    }

    for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        sumStartProbs[stateIndex] += iarray[0][stateIndex] / iarray0Sum * sumWeightAllSNPs;
        sumRatioMeans[stateIndex] += sumOfMeans[stateIndex] * sumWeightHets[stateIndex];
        sumRatioVariances[stateIndex] += sumOfVariances[stateIndex] * sumWeightHets[stateIndex];
    }
}

bool updateParameters(std::vector<double>& startProbs, std::vector<double>& ratioMeans, std::vector<double>& ratioVariances,
std::vector<std::vector<double>>& transitionProbabilities, const std::vector<double>& sumStartProbs, const std::vector<double>& sumRatioMeans,
const std::vector<double>& sumRatioVariances, const std::vector<std::vector<double>>& sumTransitionProbabilities, const std::vector<double>& sumWeightHetsAllChrom, int numStates, int totalSNPS, const double maxChange) {
    bool changed = false;
    
    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        startProbs[stateIndex] = sumStartProbs[stateIndex] / totalSNPS;
        if(fabs(ratioMeans[stateIndex] - sumRatioMeans[stateIndex] / sumWeightHetsAllChrom[stateIndex]) > maxChange || fabs(ratioVariances[stateIndex] - sumRatioVariances[stateIndex] / sumWeightHetsAllChrom[stateIndex]) > maxChange)
            changed = true;
        ratioMeans[stateIndex] = sumRatioMeans[stateIndex] / sumWeightHetsAllChrom[stateIndex];
        ratioVariances[stateIndex] = sumRatioVariances[stateIndex] / sumWeightHetsAllChrom[stateIndex];
        printf("%d The values are start %f mean %f var %f %f\n", stateIndex, startProbs[stateIndex], ratioMeans[stateIndex], ratioVariances[stateIndex], sumWeightHetsAllChrom[stateIndex]);

        for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
            if(fabs(transitionProbabilities[stateIndex][innerStateIndex] - sumTransitionProbabilities[stateIndex][innerStateIndex] / totalSNPS) > maxChange)
                changed = true;
            transitionProbabilities[stateIndex][innerStateIndex] = sumTransitionProbabilities[stateIndex][innerStateIndex] / totalSNPS;
            printf("%f\t", transitionProbabilities[stateIndex][innerStateIndex]);
        }
        printf("\n");
    }
    
    return changed;
}



void HMMSNP::baumWelch(std::vector<Individual*>& individuals, const std::string& outputFile) {
    float* unexpressedTableLocal=0;
    float* expressedTableLocal=0;
    double maxChange = 1e-05;  
    std::vector <ExpressionInfo*>  expressionSuperVector;
 
    std::vector<std::vector<double>> iarray(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    std::vector<std::vector<double>> alphas(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    std::vector<std::vector<double>> betas(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    std::vector<std::vector<double>> emissionProbabilities(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    std::vector<double> scales(MAX_EXPRESSION_SIZE);

    for(int outer = 0; outer < 100; outer++) {
        bool changed = 0;
        std::vector<double> sumStartProbs(numStates);
        std::vector<std::vector<double>> sumTransitionProbabilities(numStates, std::vector<double>(numStates));
        std::vector<double> sumRatioMeans(numStates);
        std::vector<double> sumRatioVariances(numStates);
        std::vector<double> sumWeightHetsAllChrom(numStates);
        
        int totalSNPS = 0;
 		std::vector<double> sumWeightHets(numStates);
   
        for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            sumRatioMeans[stateIndex] = 0.001;
            sumRatioVariances[stateIndex] = 0.001;
            sumStartProbs[stateIndex] = 0.001;
            for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
                sumTransitionProbabilities[stateIndex][innerStateIndex] = 0.001;
            }
        }

        for(int individualIndex = 0 ; individualIndex < numIndividuals; individualIndex++) {
            expressedTableLocal = expressedTable;
            unexpressedTableLocal = unexpressedTable;
            std::vector <Chromosome*> chromosomes = individuals[individualIndex]->getChromosomes();
            for(unsigned int chromosomeIter = 0; chromosomeIter < chromosomes.size(); chromosomeIter++) {
                std::vector <ExpressionInfo*> expressionInfos = chromosomes[chromosomeIter]->getFullExpressions(1);
                expressionSuperVector.clear();
                expressionSuperVector.insert(expressionSuperVector.end(), expressionInfos.begin(), expressionInfos.end());
    
                for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
                    for(int stateIndex  = 0; stateIndex < numStates; stateIndex++) {
                        iarray[expressionIndex][stateIndex] = 0.;
                    }
                }
                int sumWeightAllSNPs = expressionInfos.size();
                
                totalSNPS += sumWeightAllSNPs;
                if(expressionSuperVector.size() < 1) {
                    chromosomes.erase (chromosomes.begin()+chromosomeIter);
                    chromosomeIter--;
                    continue;
                }
                for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
                    scales[expressionIndex] = 0.0;
                    std::fill(betas[expressionIndex].begin(), betas[expressionIndex].end(), 0.0);
                    std::fill(alphas[expressionIndex].begin(), alphas[expressionIndex].end(), 0.0);
                    std::fill(emissionProbabilities[expressionIndex].begin(), emissionProbabilities[expressionIndex].end(), 0.0);
                }
                std::vector<double> sumsHetero(numStates);
                double iarray0Sum = 0.;
                std::vector<std::vector<double>> transitionSums(numStates, std::vector<double>(numStates));
                
                computeEmissionProbabilities(expressionSuperVector, numStates, expressedTableLocal, unexpressedTableLocal, ratioMeans, ratioVariances, runifValue, emissionProbabilities);
                forwardPart( numStates, startProbs, transitionProbabilities, emissionProbabilities, alphas, scales);
                backwardPart(expressionSuperVector, numStates, transitionProbabilities, emissionProbabilities, betas, scales);
                eStep(expressionSuperVector, numStates, alphas, betas, iarray, sumWeightHets, iarray0Sum, transitionProbabilities, emissionProbabilities, transitionSums);
                mStep(expressionSuperVector, numStates, iarray, sumWeightHets, sumStartProbs, sumRatioMeans, sumRatioVariances, iarray0Sum, sumWeightAllSNPs, transitionSums, sumTransitionProbabilities);
                for (int index = 0; index < numStates; index++) {
                    sumWeightHetsAllChrom[index] += sumWeightHets[index];
                }
            }	
        }

        changed = updateParameters(startProbs, ratioMeans, ratioVariances, transitionProbabilities, sumStartProbs, sumRatioMeans, sumRatioVariances, sumTransitionProbabilities, sumWeightHetsAllChrom, numStates, totalSNPS, maxChange);
        if(!changed) {
            printf("NO MORE CHANGE iteration %d\n", outer);
            break;
        }		
    }
 

    std::ofstream myfile(outputFile);
    
    if (myfile.is_open()) {
        myfile << numStates << "\n";
        
        for (int outerState = 0; outerState < numStates; outerState++) {
            myfile << startProbs[outerState] << "\t" << ratioMeans[outerState] << "\t" << ratioVariances[outerState];
            for (int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
                myfile << "\t" << transitionProbabilities[outerState][innerStateIndex];
            }
            myfile << std::endl;
        }
        
        myfile.close();
        std::cout << "File written successfully." << std::endl;
    } 
    else {
        std::cerr << "Unable to open file: " << outputFile << std::endl;
    }

    oneStepTransitions.resize(numStates, std::vector<double>(numStates, 0.0));
    int expressedStates = numStates-1;
    for(int firstIndex = 0; firstIndex <= expressedStates/2; firstIndex++) {
        int mirrorFirst = expressedStates-firstIndex-1;
        for(int secondIndex = 0; secondIndex < expressedStates; secondIndex++) {
            int mirrorSecond = expressedStates-secondIndex-1;
            double sum = oneStepTransitions[firstIndex][secondIndex] + oneStepTransitions[mirrorFirst][mirrorSecond];
            oneStepTransitions[firstIndex][secondIndex] = sum/2;
            oneStepTransitions[mirrorFirst][mirrorSecond] = sum/2;
	     }
    }
    for(int firstIndex = 0; firstIndex <= expressedStates/2; firstIndex++)  {
        int mirrorFirst = expressedStates-firstIndex-1;
        double sum = oneStepTransitions[firstIndex][expressedStates] + oneStepTransitions[mirrorFirst][expressedStates];
        oneStepTransitions[firstIndex][expressedStates] = sum/2;
        oneStepTransitions[mirrorFirst][expressedStates] = sum/2;
        sum = oneStepTransitions[expressedStates][firstIndex] + oneStepTransitions[expressedStates][mirrorFirst];
        oneStepTransitions[expressedStates][firstIndex] = sum/2;
        oneStepTransitions[expressedStates][mirrorFirst] = sum/2;
	}

}


void HMMSNP::initializeTransitionValues(std::vector<Individual*>& individuals, std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues) {

    std::vector<Chromosome*> tempChromosomes = individuals[0]->getChromosomes();
    ltorTransitions = std::vector<std::vector<std::vector<std::vector<double>>>>(tempChromosomes.size());
    for (unsigned int chromosomeIndex = 0; chromosomeIndex < tempChromosomes.size(); chromosomeIndex++) {
        ltorTransitions[chromosomeIndex] = std::vector<std::vector<std::vector<double>>>();
    }

    for (unsigned int chromosomeIter = 0; chromosomeIter < tempChromosomes.size(); chromosomeIter++) {

        int size = individuals[0]->getChromosomes()[chromosomeIter]->getFullExpressions(1).size();
        tempTransitionValues[chromosomeIter] = std::vector<std::vector<std::vector<double>>>(size);

        for (int SNPiter = 0; SNPiter < size; SNPiter++) {

            ltorTransitions[chromosomeIter][SNPiter] = std::vector<std::vector<double>>(numStates);
            tempTransitionValues[chromosomeIter][SNPiter] = std::vector<std::vector<double>>(numStates);

            for (int k = 0; k < numStates; k++) {
                ltorTransitions[chromosomeIter][SNPiter][k] = std::vector<double>(numStates);
                for (int l = 0; l < numStates; l++) {
                    ltorTransitions[chromosomeIter][SNPiter][k][l] = 0;
                }

                tempTransitionValues[chromosomeIter][SNPiter][k] = std::vector<double>(numStates);

                for (int l = 0; l < numStates; l++) {
                    if (transitionProbabilities[k][l] == 0) {
                        transitionProbabilities[k][l] = 1e-5;
                    }
                    if (ltorTransitions[chromosomeIter][SNPiter][k][l] == 0) {
                        ltorTransitions[chromosomeIter][SNPiter][k][l] = transitionProbabilities[k][l];
                    }
                    tempTransitionValues[chromosomeIter][SNPiter][k][l] = 0.0;
                }
            }
        }
    }
}

void HMMSNP::perChromosomeSet(std::vector<Individual*>& individuals, std::vector<std::vector<int>>& maxStates, std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues, std::vector<std::vector<float>>& iarray, int chromosomeIndex) {
    int size = individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1).size();

    maxStates.resize(numIndividuals);
    for (int i = 0; i < numIndividuals; i++) {
        maxStates[i].resize(size);
    }

    std::vector<ExpressionInfo*> expressions = individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1);

    for (int expressionIndex = 0; expressionIndex < size; expressionIndex++) {
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            iarray[expressionIndex][stateIndex] = 0.0;
        }

        for (int innerIndex = 0; innerIndex < numStates; innerIndex++) {
            for (int outerIndex = 0; outerIndex < numStates; outerIndex++) {
                tempTransitionValues[chromosomeIndex][expressionIndex][innerIndex][outerIndex] = 0.0;
            }
        }
    }
}

  

void HMMSNP::setEmissionProbabilities(std::vector<Individual*>& individuals, std::vector<std::vector<double>>& emissionProbabilities, 
std::vector<ExpressionInfo*>& expressionSuperVector, int individualIndex, int chromosomeIndex, float summationRatio) {
    runifValue = summationRatio;

    for (unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
        int logRatioIndex = static_cast<int>((expressionSuperVector[expressionIndex]->getExpressionLogRatio() + 10) * 10);
        if (logRatioIndex < 0) {
            logRatioIndex = 0;
        }
        if (logRatioIndex > 139) {
            logRatioIndex = 139;
        }
        double emissionSummation = 0.0;
        for (int emitIndex = 0; emitIndex < numStates; emitIndex++) {
            emissionProbabilities[expressionIndex][emitIndex] = 1.0;

            if (emitIndex < numStates - 1 && (numStates % 2 == 0)) {
                emissionProbabilities[expressionIndex][emitIndex] *= expressedTable[logRatioIndex];
            } else if (numStates % 2 == 0) {
                emissionProbabilities[expressionIndex][emitIndex] *= unexpressedTable[logRatioIndex];
            }

            if (expressionSuperVector[expressionIndex]->isHeterozygote()) {
                emissionProbabilities[expressionIndex][emitIndex] *= computeProbability(ratioMeans[emitIndex], ratioVariances[emitIndex], expressionSuperVector[expressionIndex]->getRatioOfRatios());
            }

            emissionProbabilities[expressionIndex][emitIndex] *= 10;
            emissionSummation += emissionProbabilities[expressionIndex][emitIndex];
            emissionProbabilities[expressionIndex][emitIndex] += runifValue;
        }
    }
}

      
 


void HMMSNP::backward(std::vector<Individual*>& individuals, std::vector<std::vector<double>>& betas, std::vector<double>& scales, 
std::vector<std::vector<std::vector<double>>>& transitionsLocal, std::vector<std::vector<double>>& emissionProbabilities, int size, 
int individualIndex, std::vector<ExpressionInfo*>& expressionSuperVector) {
    int individualToPrint = -6;
    int minToPrint = 0;
    int maxToPrint = 100000000;

    for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        betas[size - 1][stateIndex] = 1;
    }

    for (int orderingIndex = size - 2; orderingIndex >= 0; orderingIndex--) {
        int betaIndex = orderingIndex;
        int nextBetaIndex = orderingIndex + 1;
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            betas[betaIndex][stateIndex] = 0.0;

            for (int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
                betas[betaIndex][stateIndex] += transitionsLocal[betaIndex][stateIndex][innerStateIndex] * betas[nextBetaIndex][innerStateIndex] * emissionProbabilities[nextBetaIndex][innerStateIndex];
                if (individualIndex == individualToPrint && expressionSuperVector[orderingIndex]->getSNP()->getLocation() >= minToPrint && expressionSuperVector[orderingIndex]->getSNP()->getLocation() <= maxToPrint) {
                    printf("individual %d For  emit %d innerEmit %d betaIndex %d  and   from loc %d we now have ratio %f exp %f transition %f beta next %f \n", individualIndex, stateIndex, innerStateIndex, betaIndex, expressionSuperVector[betaIndex]->getSNP()->getLocation(), expressionSuperVector[nextBetaIndex]->getRatioOfRatios(), expressionSuperVector[nextBetaIndex]->getExpressionLogRatio(), transitionsLocal[betaIndex][stateIndex][innerStateIndex], betas[nextBetaIndex][innerStateIndex]);
                }
            }
        }

        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
            betas[betaIndex][stateIndex] /= scales[betaIndex];

            if (std::isnan(betas[betaIndex][stateIndex])) {
                printf("This is a problem in beta %d %d\n", betaIndex, stateIndex);
            } else if (betas[betaIndex][stateIndex] == 0.0) {
                printf("we are at EXACTLY 0  %d %d from %f %f now at %f \n", betaIndex, stateIndex, betas[betaIndex + 1][stateIndex], scales[betaIndex], betas[betaIndex][stateIndex]);
            }
        }
    }
}



void HMMSNP::viterbiTraceback(std::vector<std::vector<double>>& alphas, std::vector<std::vector<int>>& maxStates, std::vector<std::vector<int>>& viterbiStates, int individualIndex) {
    int maxState = 1;
    float maxVal = -1;

    for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {
        if (alphas.back()[stateIndex] > maxVal) {
            maxState = stateIndex;
            maxVal = alphas.back()[stateIndex];
        }
    }

    for (int orderingIndex = alphas.size() - 1; orderingIndex >= 1; orderingIndex--) {
        maxStates[individualIndex][orderingIndex] = maxState;
        maxState = viterbiStates[orderingIndex][maxState];
    }
    maxStates[individualIndex][0] = maxState;
}




void HMMSNP::fileOutput(std::vector<Individual*>& individuals, std::ofstream& myViterbiFile, std::ofstream& myfile, int chromosomeIndex, 
std::vector<ExpressionInfo*>& expressionSuperVector, std::vector<std::vector<float>>& iarray, std::vector<std::vector<int>>& maxStates, 
int individualIndex, const std::string& fileName) {
    myfile << "track type=bedGraph name=ind" << individuals[individualIndex]->getName() << "chr" << chromosomeIndex << std::endl;
    myViterbiFile << "track type=bedGraph name=ind" << individuals[individualIndex]->getName() << "viterbi" << "chr" << chromosomeIndex << std::endl;
    for (unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
        myViterbiFile << "chr";
        myViterbiFile << chromosomeIndex << " ";
        myViterbiFile << expressionSuperVector[expressionIndex]->getSNP()->getLocation() << " ";
        myViterbiFile << expressionSuperVector[expressionIndex]->getSNP()->getLocation() + 1 << " ";
        myViterbiFile << maxStates[individualIndex][expressionIndex];

        myfile << "chr" << chromosomeIndex << " " << expressionSuperVector[expressionIndex]->getSNP()->getLocation() << " " << 
            expressionSuperVector[expressionIndex]->getSNP()->getLocation() + 1 << " ";
        float eval = 0;
        for (int stateIndex = 0; stateIndex < numStates; stateIndex++) {

            eval += iarray[expressionIndex][stateIndex] * ratioMeans[stateIndex];
        }
        myfile << eval << "\n";
        myViterbiFile << "\n";
    }
}


void HMMSNP::addToTempTransitions(std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues, std::vector<std::vector<double>>& alphas, 
std::vector<std::vector<double>>& emissionProbabilities, std::vector<std::vector<double>>& betas, std::vector<std::vector<std::vector<double>>>& transitionsLocal, 
std::vector<ExpressionInfo*>& expressionSuperVector, int chromosomeIndex) {
    for (unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size() - 1; expressionIndex++) {
        if (expressionSuperVector[expressionIndex + 1]->isHeterozygote() || numStates % 2 == 0) {
            std::vector<std::vector<double>> transitionSums(numStates, std::vector<double>(numStates));
            std::vector<double> transitionFullSums(numStates, 0.0);
            for (int firstState = 0; firstState < numStates; firstState++) {
                for (int secondState = 0; secondState < numStates; secondState++) {
                    transitionSums[firstState][secondState] = alphas[expressionIndex][firstState] * transitionsLocal[expressionIndex][firstState][secondState] * emissionProbabilities[expressionIndex + 1][secondState] * betas[expressionIndex + 1][secondState];
                    transitionFullSums[firstState] += transitionSums[firstState][secondState];
                    tempTransitionValues[chromosomeIndex][expressionIndex][firstState][secondState] += transitionSums[firstState][secondState];
                }
            }
        }
    }
}


void HMMSNP::updateTransitions(std::vector<std::vector<std::vector<std::vector<double>>>>& tempTransitionValues, float weight, int size, int chromosomeIndex) {
    for(int expressionIndex = 0; expressionIndex < size; expressionIndex++) {
        std::vector<std::vector<double>>& distanceDependentTransitions = transitionProbabilities;
        std::vector<double> sums(numStates);
        int expressedStates = numStates - 1;
        for(int firstIndex = 0; firstIndex < expressedStates; firstIndex++) {
            sums[firstIndex] = 0.;
            for(int secondIndex = 0; secondIndex < numStates; secondIndex++) {
                sums[firstIndex] +=  tempTransitionValues[chromosomeIndex][expressionIndex][firstIndex][secondIndex];
            }
        }
        for(int firstIndex = 0; firstIndex <= expressedStates/2; firstIndex++) {
            int mirrorFirst = expressedStates-firstIndex-1;
            for(int secondIndex = 0; secondIndex < expressedStates; secondIndex++) {
                int mirrorSecond = expressedStates-secondIndex-1;
                ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] = tempTransitionValues[chromosomeIndex][expressionIndex][firstIndex][secondIndex] + tempTransitionValues[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond];
                ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] += weight * (distanceDependentTransitions[firstIndex][secondIndex] + distanceDependentTransitions[mirrorFirst][mirrorSecond]);
                ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] /= (sums[firstIndex] + sums[mirrorFirst] + weight*2);
                ltorTransitions[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond] = ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex];
            }
        }
    }
}


void HMMSNP::leftToRight(std::vector<Individual*>& individuals, std::string outputFileName, double weight) {
    std::ofstream myfile;
    std::ofstream myViterbiFile;

    int maxOuter = 100;

    std::vector<Chromosome*> tempChromosomes = individuals[0]->getChromosomes();
    std::vector<std::vector<std::vector<int>>> minMaxIndices(tempChromosomes.size(), std::vector<std::vector<int>>(MAX_EXPRESSION_SIZE, std::vector<int>(2)));
    std::vector<std::vector<std::vector<std::vector<double>>>> tempTransitionValues(tempChromosomes.size(), std::vector<std::vector<std::vector<double>>>(numStates, std::vector<std::vector<double>>(numStates, std::vector<double>(numStates))));
    initializeTransitionValues(individuals, tempTransitionValues);

    std::vector<std::vector<float>> iarray(MAX_EXPRESSION_SIZE, std::vector<float>(numStates));
    std::vector<std::vector<double>> alphas(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    std::vector<std::vector<double>> betas(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    std::vector<std::vector<std::vector<int>>> viterbiStates(numIndividuals, std::vector<std::vector<int>>(MAX_EXPRESSION_SIZE, std::vector<int>(numStates)));
    std::vector<std::vector<double>> emissionProbabilities(MAX_EXPRESSION_SIZE, std::vector<double>(numStates));
    std::vector<double> scales(MAX_EXPRESSION_SIZE);

    for (int individualIndex = 0; individualIndex < numIndividuals; individualIndex++) {
        viterbiStates[individualIndex] = std::vector<std::vector<int>>(MAX_EXPRESSION_SIZE, std::vector<int>(numStates));
    }

    for (int i = 0; i < MAX_EXPRESSION_SIZE; i++) {
        iarray[i] = std::vector<float>(numStates);
        alphas[i] = std::vector<double>(numStates);
        betas[i] = std::vector<double>(numStates);
        emissionProbabilities[i] = std::vector<double>(numStates);
    }

    for (int outer = 0; outer < maxOuter; outer++) {
        std::cout << "ltor loop:" << outer << std::endl;

        for (unsigned int chromosomeIndex = 0; chromosomeIndex < tempChromosomes.size(); chromosomeIndex++) {
            std::vector<std::vector<std::vector<double>>> transitionsLocal = ltorTransitions[chromosomeIndex];
            int sizeOfArray = individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1).size();
            if (sizeOfArray < 1) {
                continue;
            }
            std::vector<std::vector<int>> maxStates(numIndividuals, std::vector<int>());
            perChromosomeSet(individuals, maxStates, tempTransitionValues, iarray, chromosomeIndex);
            for (int individualIndex = 0; individualIndex < numIndividuals; individualIndex++) {
                std::vector<Chromosome*> chromosomes = individuals[individualIndex]->getChromosomes();
                std::vector<ExpressionInfo*> expressionSuperVector = chromosomes[chromosomeIndex]->getFullExpressions(1);
                for (unsigned int index = 0; index < expressionSuperVector.size(); index++) {
                    scales[index] = 0.0;
                    for (int innerIndex = 0; innerIndex < numStates; innerIndex++) {
                        betas[index][innerIndex] = 0.0;
                        alphas[index][innerIndex] = 0.0;
                        emissionProbabilities[index][innerIndex] = 0.0;
                        viterbiStates[individualIndex][index][innerIndex] = 0;
                    }
                }
                setEmissionProbabilities(individuals, emissionProbabilities, expressionSuperVector, individualIndex, chromosomeIndex, runifValue);
                forward(alphas, emissionProbabilities, scales, transitionsLocal, viterbiStates[individualIndex]);
                backward(individuals, betas, scales, transitionsLocal, emissionProbabilities, expressionSuperVector.size(), individualIndex, expressionSuperVector);
                viterbiTraceback(alphas, maxStates, viterbiStates[individualIndex], individualIndex);
                addToTempTransitions(tempTransitionValues, alphas, emissionProbabilities, betas, transitionsLocal, expressionSuperVector, chromosomeIndex);
                if (!outputFileName.empty() && (outer == maxOuter - 1 || outer % 99 == 0)) {
                    std::ofstream myfile;
                    std::ofstream myViterbiFile;

                    if (!outputFileName.empty() && (outer == maxOuter - 1 || outer % 99 == 0)) {
                        std::ofstream myfile;
                        std::ofstream myViterbiFile;
                        std::ostringstream evalNameStream;
                        evalNameStream << outputFileName << "chr" << chromosomeIndex << "iter" << outer << ".evals";
                        std::string evalName = evalNameStream.str();

                        std::ostringstream viterbiNameStream;
                        viterbiNameStream << outputFileName << "chr" << chromosomeIndex << "iter" << outer << ".viterbi";
                        std::string viterbiName = viterbiNameStream.str();

                        myfile.open(evalName);
                        myViterbiFile.open(viterbiName);

                        for (int individualIndex = 0; individualIndex < numIndividuals; ++individualIndex) {
                            fileOutput(individuals, myViterbiFile, myfile, chromosomeIndex, expressionSuperVector,
                                        iarray, maxStates, individualIndex, outputFileName);
                        }

                        myfile.close();
                        myViterbiFile.close();
                    }

                    updateTransitions(tempTransitionValues, weight, individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1).size(), chromosomeIndex);
                }
            }
        }
    }
}
