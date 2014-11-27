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

int numChromosomes2;
using namespace std;


const int PI =  3.14159265;


double HMMSNP::computeProbability(double mean, double variance, double value) {
  double std = sqrt(variance);
  
  return(1/(std*sqrt(2*PI))*exp(-1*(value-mean)*(value-mean)/(2*variance)));


}




HMMSNP::~HMMSNP()
{

	for(int i = 0 ; i < numStates; i++) {
 
	delete[]transitionProbabilities[i];

	}
	delete []ratioMeans;
	delete []ratioVariances;
	
	delete []transitionProbabilities;
	delete []startProbs;


}





HMMSNP::HMMSNP(char* fileName) {
  ltorTransitions = 0; 
  oneStepTransitions = 0;
  useLow = 1;  

 ifstream infile;
  runifValue = 0.01;
  //    printf("Called with %s\n", fileName);
  infile.open(fileName, ifstream::in);
  char line[100000];

  infile.getline(line, 100000);
  numStates = atoi(strtok(line, " \t"));
  
  cout << "numstates " << numStates << endl;
  ratioMeans = new double[numStates];
  ratioVariances = new double[numStates];
  transitionProbabilities = new double*[numStates];
  for(int i = 0; i < numStates; i++) {
    transitionProbabilities[i] = new double[numStates];
  } 
  startProbs = new double[numStates];
  
  for(int stateIndex = 0 ; stateIndex < numStates; stateIndex++) {
    infile.getline(line,100000);
    cout << "line " << line << endl;
    startProbs[stateIndex] = atof(strtok(line, " \t"));
    cout << "set start probs " << startProbs[stateIndex] << endl;
    ratioMeans[stateIndex] = atof(strtok(NULL, " \t"));
    cout << "mean " << ratioMeans[stateIndex] << endl;
    ratioVariances[stateIndex] = atof(strtok(NULL, " \t"));
    cout << "variance " << ratioVariances[stateIndex] << endl;
    for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
      double val = atof(strtok(NULL, " \t"));
      cout << "obtained transition " << val << endl;
      transitionProbabilities[stateIndex][innerStateIndex] = val;
      
    }
  }
}
 


void HMMSNP::baumWelch(Individual** individuals, char* outputFile) {


 float* unexpressedTableLocal=0;
 float* expressedTableLocal=0;
 double maxChange = 1e-05;  
 vector <ExpressionInfo*>  expressionSuperVector;
 
 int size = individuals[0]->getChromosomes().size();
 
 float  **iarray = new float*[90000];
 double **alphas = new double*[90000]   ;
 double **betas = new double*[90000];
 double **emissionProbabilities = new double*[90000];
 double *scales  = new double[90000];

 for(int i = 0; i < 90000; i++) {
   iarray[i] = new float[numStates];
   alphas[i] = new double[numStates];
   betas[i] = new double[numStates];
   emissionProbabilities[i] = new double[numStates];
 }


 
 
 for(int outer = 0; outer < 100; outer++) {
   bool changed = 0;
   double sumStartProbs[numStates];
   double sumTransitionProbabilities[numStates][numStates];
   double sumRatioMeans[numStates];
   double sumRatioVariances[numStates];
   int totalSNPS = 0;
 		
   float* totalHets = new float[numStates]();
   
   for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
     sumRatioMeans[stateIndex] = 0.001;
     totalHets[stateIndex] = 0.1;
     sumRatioVariances[stateIndex] = 0.001;
     sumStartProbs[stateIndex] = 0.001;
     for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {

       sumTransitionProbabilities[stateIndex][innerStateIndex] = 0.001;
       
       
     }
     
   }

	       
   int currentIndividualToUse = 0; 
   for(unsigned int individualIndex = 0 ; individualIndex < numIndividuals; individualIndex++) {
      
       expressedTableLocal = expressedTable;
       unexpressedTableLocal = unexpressedTable;
   
       vector <Chromosome*> chromosomes = individuals[individualIndex]->getChromosomes();
       for(unsigned  int chromosomeIter = 0; chromosomeIter < chromosomes.size(); chromosomeIter++) {
	 
	 vector <ExpressionInfo*> expressionInfos = chromosomes[chromosomeIter]->getFullExpressions(1);
	 expressionSuperVector.clear();
	 expressionSuperVector.insert(expressionSuperVector.end(), expressionInfos.begin(), expressionInfos.end());
	 
	 for(int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
	   
	   for(int stateIndex  = 0; stateIndex < numStates; stateIndex++) {
	     iarray[expressionIndex][stateIndex] = 0.;
	   }
	 }
	 int   sumWeightAllSNPs = expressionInfos.size();
	 float* sumWeightHets = new float[numStates]();
	 
	 totalSNPS += sumWeightAllSNPs;
	 if(expressionSuperVector.size() < 1) {
	   
	 chromosomes.erase (chromosomes.begin()+chromosomeIter);
	 chromosomeIter--;
	 continue;
       }
	 
	 
       
	 
	 for( int index = 0; index < expressionSuperVector.size(); index++) {
	   
	 scales[index] = 0.0;
	 for( int innerIndex = 0; innerIndex < numStates; innerIndex++) {
	   betas[index][innerIndex] = 0.0;
	   
	   alphas[index][innerIndex] = 0.0;
	   emissionProbabilities[index][innerIndex] = 0.0;
	   
	 }
	 
       }
	 for( int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
	   
	   
	   
	   
	   
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
	       
	       emissionProbabilities[expressionIndex][emitIndex] *= computeProbability(ratioMeans[emitIndex], ratioVariances[emitIndex], expressionSuperVector[expressionIndex]->getRatioOfRatios() ) ;
			
	       
	       
	     }
	     
	     emissionProbabilities[expressionIndex][emitIndex] *= 10;
	     emissionProbabilities[expressionIndex][emitIndex] += runifValue;
	     
	     
	     
	   }
	   
	 }

	

      
	 
       //FORWARD PART
	 double sumAlphas = 0.0;
	 for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	   alphas[0][stateIndex] = startProbs[stateIndex] * emissionProbabilities[0][stateIndex];
	   betas[expressionSuperVector.size()-1][stateIndex] = 1.0/numStates;
	   sumAlphas += alphas[0][stateIndex];
	 }
	 scales[0] = sumAlphas;
	 for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
	   alphas[0][emitIndex] /= sumAlphas;
		    
	 }
	 for(unsigned int alphaIndex = 1; alphaIndex < expressionSuperVector.size(); alphaIndex++) {
	   double sumAlphas = 0.0;
	   
	   
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
	 

	 // backward part
	 for( int betaIndex = expressionSuperVector.size()-2; betaIndex >= 0; betaIndex--) {
	   
	   for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	     betas[betaIndex][stateIndex] = 0.0;
	     for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
			
	       betas[betaIndex][stateIndex] += transitionProbabilities[stateIndex][innerStateIndex]*betas[betaIndex+1][innerStateIndex]*emissionProbabilities[betaIndex+1][innerStateIndex];
	       
	     }
		      

	   }
	   for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	     betas[betaIndex][stateIndex] /= scales[betaIndex];
 		
	   }
	 }
 		
		  
 		  // E STEP
	 double sumsHetero[numStates];
	 
	 for(unsigned int i = 0; i < expressionSuperVector.size(); i++) {
	   
	   
	   for(int j = 0; j < numStates; j++) {
	     iarray[i][j] = 0.00;
	     sumsHetero[j] = 0.0;
	     
	   }
	   

	 }

		  
	 for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
	   
	   double sumIndex = 0.0;
	   for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	     
	     iarray[expressionIndex][stateIndex] = alphas[expressionIndex][stateIndex] * betas[expressionIndex][stateIndex] ;
	     
	     sumIndex += iarray[expressionIndex][stateIndex];
	     
	   }
	   
	   for(int stateIndex  =0; stateIndex < numStates; stateIndex++) {
	     iarray[expressionIndex][stateIndex] /= sumIndex;
	     
	     
	   }	  
	
	   
	  
	   
	   for(int stateIndex  =0; stateIndex < numStates; stateIndex++) {
	     
	     if(expressionSuperVector[expressionIndex]->isHeterozygote() && iarray[expressionIndex][stateIndex] > 0.05) {
	       sumsHetero[stateIndex] += iarray[expressionIndex][stateIndex];
	       sumWeightHets[stateIndex]+= iarray[expressionIndex][stateIndex];
	       totalHets[stateIndex]+= iarray[expressionIndex][stateIndex];
		
	     }
		
	   }
			  
			  
	 }
			
	 double iarray0Sum = 0.0;
	 for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	   
	   iarray0Sum += iarray[0][stateIndex];
	 }
	 double transitionSums[numStates][numStates];
	 double transitionFullSums[numStates];
	 for(int firstState = 0 ; firstState < numStates; firstState++) {
	   transitionFullSums[firstState] = 0.0;
	   for(int secondState = 0; secondState < numStates; secondState++) {
	     
	     transitionSums[firstState][secondState] = transitionProbabilities[firstState][secondState] * 1e-05;
	     
	   }
	   
	 }
	 for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size()-1; expressionIndex++) {
	   for(int firstState = 0; firstState < numStates; firstState++) {
	     if( iarray[expressionIndex][firstState] > 0.05 || firstState == numStates-1) {
	       for(int secondState = 0; secondState < numStates; secondState++) {
		 
		 transitionSums[firstState][secondState] += alphas[expressionIndex][firstState]*transitionProbabilities[firstState][secondState] *emissionProbabilities[expressionIndex+1][secondState]*betas[expressionIndex+1][secondState];
				  
	       }
	       
				
	     }
	     
	   }
	   
	 }
       
	 // M STEP


			       
	 for(int firstState = 0; firstState < numStates; firstState++) {
	   
	   for(int secondState = 0; secondState < numStates; secondState++) {
			  
	     transitionFullSums[firstState] += transitionSums[firstState][secondState];
	     
	   }
	 }
	 for(int firstState = 0; firstState < numStates; firstState++) {
	   for(int secondState = 0; secondState < numStates; secondState++) {
	     
	     sumTransitionProbabilities[firstState][secondState] += transitionSums[firstState][secondState]/transitionFullSums[firstState] * expressionSuperVector.size();
					
	     
	   }
	   
	 }
				
	 double sumOfMeans[numStates];
	 double sumOfVariances[numStates];
	 double sumOfMeansExpressions[numStates];
	 double sumOfVariancesExpressions[numStates];
	 for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	   sumOfMeans[stateIndex] = 0.0;
	   sumOfVariances[stateIndex] = 0.001;
	   sumOfMeansExpressions[stateIndex] = 0.0;
	   sumOfVariancesExpressions[stateIndex] = 0.0;
	   for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
	     if(expressionSuperVector[expressionIndex]->isHeterozygote() && iarray[expressionIndex][stateIndex] > 0.05 ) {
	       sumOfMeans[stateIndex]+= expressionSuperVector[expressionIndex]->getRatioOfRatios() * iarray[expressionIndex][stateIndex] 	/ sumsHetero[stateIndex];
		
	     }
	     else if(expressionSuperVector[expressionIndex]->isHeterozygote() && stateIndex == numStates-1 && iarray[expressionIndex][stateIndex] > 0.0 ) {
				  sumOfMeans[stateIndex]+= expressionSuperVector[expressionIndex]->getRatioOfRatios() * iarray[expressionIndex][stateIndex]  / sumsHetero[stateIndex];
				  
	     }
	     
	     
	   }
	   
	 }

				
	 for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	   for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
	     if(expressionSuperVector[expressionIndex]->isHeterozygote() && iarray[expressionIndex][stateIndex] > 0.05) {
	       sumOfVariances[stateIndex] +=(expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex]) * (expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex])  * iarray[expressionIndex][stateIndex] / sumsHetero[stateIndex];
	     }

	     if(expressionSuperVector[expressionIndex]->isHeterozygote() && stateIndex == numStates-1 && iarray[expressionIndex][stateIndex] > 0.0) {
	       sumOfVariances[stateIndex] += (expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex]) * (expressionSuperVector[expressionIndex]->getRatioOfRatios() - sumOfMeans[stateIndex]) * iarray[expressionIndex][stateIndex] / sumsHetero[stateIndex];
	       

	       
	       
	     }

	   }
	   
	 }
       
	 for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
	   
	   sumStartProbs[stateIndex] += iarray[0][stateIndex]/iarray0Sum * sumWeightAllSNPs;
	   
	   
	   sumRatioMeans[stateIndex] += sumOfMeans[stateIndex] * sumWeightHets[stateIndex];
		
	   sumRatioVariances[stateIndex] += sumOfVariances[stateIndex] * sumWeightHets[stateIndex];
	   
		
	 }
						
 		
		
       }
	
   }
		
		
   for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
     startProbs[stateIndex] = sumStartProbs[stateIndex] / totalSNPS;
     if( fabs(ratioMeans[stateIndex] - sumRatioMeans[stateIndex]/totalHets[stateIndex]) > maxChange || fabs(ratioVariances[stateIndex] - sumRatioVariances[stateIndex]/totalHets[stateIndex]) > maxChange) 
       changed = 1;                    
     ratioMeans[stateIndex] = sumRatioMeans[stateIndex] / totalHets[stateIndex];
     ratioVariances[stateIndex] = sumRatioVariances[stateIndex] / totalHets[stateIndex];
     printf("%d The values are start %f mean %f var %f %f\n", stateIndex, startProbs[stateIndex], ratioMeans[stateIndex], ratioVariances[stateIndex], totalHets[stateIndex]);  
	
     
     for(int innerStateIndex =0; innerStateIndex < numStates; innerStateIndex++) {
       if(fabs(transitionProbabilities[stateIndex][innerStateIndex] - sumTransitionProbabilities[stateIndex][innerStateIndex]/totalSNPS) > maxChange) 
	 changed = 1;
       transitionProbabilities[stateIndex][innerStateIndex] = sumTransitionProbabilities[stateIndex][innerStateIndex] / totalSNPS;
       printf("%f\t", transitionProbabilities[stateIndex][innerStateIndex]);
     }
     printf("\n");
     
		  
   }
		

   if(!changed) {
 	
     printf("NO MORE CHANGE iteration %d\n", outer);
     break;
   }
		  
		
 }
 
 for(int i = 0; i < 90000; i++) {
   delete[] iarray[i];
   delete[] alphas[i];
   delete[] betas[i];
   delete[] emissionProbabilities[i];
 }
 delete []iarray;
 delete []alphas;
 delete[]scales;
 delete[]betas;
 delete [] emissionProbabilities;
 
 
					   
 
 
 ofstream myfile;
 if(outputFile) {
   myfile.open(outputFile);
   myfile << numStates << "\n";
 
 for(int outerState = 0; outerState < numStates; outerState++) {
   myfile << startProbs[outerState] << "\t" << ratioMeans[outerState] << "\t" << ratioVariances[outerState]  ;
   for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
     myfile << "\t" << transitionProbabilities[outerState][innerStateIndex]  ;
     
   }
   myfile << endl;
 }
 
 oneStepTransitions = new double*[numStates];
 for(int i = 0; i < numStates; i++) {
   oneStepTransitions[i] = new double[numStates];
 }
 char hello[100];
 myfile.close();
 }
 /*sprintf(hello, "%s.1SNP", outputFile);
 myfile.precision(14);
 myfile.open(hello);	
 utility::rootOfAMatrix(transitionProbabilities, 1./3817., oneStepTransitions, numStates) ;
 
 
 for(int outerState = 0; outerState < numStates; outerState++) {
   myfile << startProbs[outerState] << "\t" << ratioMeans[outerState] << "\t" << ratioVariances[outerState] << "\t" ;
   for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {
     myfile << oneStepTransitions[outerState][innerStateIndex] << "\t";
    
     
   }
   myfile << endl;
 }  
 */
 
 




 
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
 



 cout << "DONE BAUM " << endl;

}





void HMMSNP::initializeTransitionValues( Individual** individuals, double **** tempTransitionValues) {

  vector <Chromosome*> tempChromosomes = individuals[0]->getChromosomes();
  if(!ltorTransitions) {
    ltorTransitions = new double ***[tempChromosomes.size()]();
    for(int i = 0; i < tempChromosomes.size(); i++) {
      ltorTransitions[i] = 0;
      tempTransitionValues[i] = 0;      
    }
  }

  for(unsigned int chromosomeIter = 0; chromosomeIter < tempChromosomes.size(); chromosomeIter++) {

    int size = individuals[0]->getChromosomes()[chromosomeIter]->getFullExpressions(1).size();
   
    if(!ltorTransitions[chromosomeIter]) {

      ltorTransitions[chromosomeIter] = new double **[size]();
      for(int i = 0; i < size; i++) {
        ltorTransitions[chromosomeIter][i] = 0;
      }
   
    }


    tempTransitionValues[chromosomeIter] = new double **[size];

    for(int SNPiter = 0; SNPiter < size; SNPiter++) {



      if(!ltorTransitions[chromosomeIter][SNPiter]) {
        ltorTransitions[chromosomeIter][SNPiter] = new double *[numStates]();

        for(int i = 0; i < numStates; i++) {

          ltorTransitions[chromosomeIter][SNPiter][i] = 0;
        }
      }
      tempTransitionValues[chromosomeIter][SNPiter] = new double*[numStates];




      for(int k = 0; k < numStates; k++) {
        if(!ltorTransitions[chromosomeIter][SNPiter][k]) {
          ltorTransitions[chromosomeIter][SNPiter][k] = new double [numStates]();
	  for(int l = 0; l < numStates;l++) {
	    ltorTransitions[chromosomeIter][SNPiter][k][l] = 0;

	  }
        }
	tempTransitionValues[chromosomeIter][SNPiter][k] = new double[numStates];
	

	
	

        for(int l = 0 ; l < numStates; l++) {
          if(transitionProbabilities[k][l]  == 0) {
            transitionProbabilities[k][l] = 1e-5;
          }
          if(!ltorTransitions[chromosomeIter][SNPiter][k][l]) {
            ltorTransitions[chromosomeIter][SNPiter][k][l] = transitionProbabilities[k][l];
            
          }
	  
          tempTransitionValues[chromosomeIter][SNPiter][k][l] = 0.0;

        }
      }

    }

  }
}






void HMMSNP::perChromosomeSet(Individual** individuals, int** maxStates, double**** tempTransitionValues,  float** iarray, int chromosomeIndex) {
 
  int size = individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1).size(); 

  for(int i = 0; i <numIndividuals; i++) {
    maxStates[i] = new int [size];
  }
  
  vector <ExpressionInfo*> expressions = individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1);
  
  

  
  
  for(int expressionIndex = 0; expressionIndex < size; expressionIndex++) {
    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
      
	iarray[expressionIndex][stateIndex] = 0.0;
      
    }    

    for(int innerIndex = 0; innerIndex < numStates; innerIndex++) {
      for(int outerIndex = 0; outerIndex < numStates; outerIndex++) {
	
	
	tempTransitionValues[chromosomeIndex][expressionIndex][innerIndex][outerIndex] = 0.;
      }
    
    }
  }
}

  



void HMMSNP::setEmissionProbabilities(  Individual** individuals, double **emissionProbabilities, vector <ExpressionInfo*> expressionSuperVector, int individualIndex, int chromosomeIndex,  float summationRatio) {

  runifValue = summationRatio;
  float* unexpressedTableLocal = unexpressedTable;
  float* expressedTableLocal = expressedTable;

  for(unsigned int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
    
    int logRatioIndex = (int)((expressionSuperVector[expressionIndex]->getExpressionLogRatio() + 10) *10);
    if(logRatioIndex < 0) {
      logRatioIndex = 0;
    }
    if(logRatioIndex > 139) {
      logRatioIndex = 139;
    }
    float emissionSummation  = 0.0;
    for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
      emissionProbabilities[expressionIndex][emitIndex] = 1.0;
      
      
      
      if(emitIndex < numStates-1 && (numStates %2 == 0)) {
	
	emissionProbabilities[expressionIndex][emitIndex]  *= expressedTableLocal[logRatioIndex];
      }
      else if(numStates %2 == 0) {
	emissionProbabilities[expressionIndex][emitIndex]  *= unexpressedTableLocal[logRatioIndex];
      }
      
      
      
      if(expressionSuperVector[expressionIndex]->isHeterozygote()) {
	emissionProbabilities[expressionIndex][emitIndex] *= computeProbability(ratioMeans[emitIndex], ratioVariances[emitIndex], expressionSuperVector[expressionIndex]->getRatioOfRatios() ) ;
      }
      
      
      
      
      
      

	
	emissionProbabilities[expressionIndex][emitIndex] *= 10;
	emissionSummation +=   emissionProbabilities[expressionIndex][emitIndex];
	emissionProbabilities[expressionIndex][emitIndex] += runifValue;

    }
  }
}
 
      
      
 



void HMMSNP::forward(Individual** individuals, double** alphas, double **emissionProbabilities, double* scales, double*** transitionsLocal, int **viterbiStates, int size, int individualIndex, vector<ExpressionInfo*> expressionSuperVector) {

  double sumAlphas = 0.0;

 
 
  for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {

    alphas[0][stateIndex] = startProbs[stateIndex] * emissionProbabilities[0][stateIndex];
    sumAlphas += alphas[0][stateIndex];
  }


  scales[0] = sumAlphas;
  for(int emitIndex = 0; emitIndex < numStates; emitIndex++) {
    alphas[0][emitIndex] /= sumAlphas;
  }





 
  for(unsigned int orderingIndex = 1; orderingIndex < size; orderingIndex++) {
    int alphaIndex = orderingIndex;
    int previousAlphaIndex = orderingIndex-1;
    

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








void HMMSNP::backward(Individual **individuals, double** betas, double* scales, double*** transitionsLocal, double** emissionProbabilities, int size, int individualIndex, vector<ExpressionInfo*> expressionSuperVector) {

  int individualToPrint = -6;
  int minToPrint = 0;
  int maxToPrint = 100000000;
  

  for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
    betas[size-1][stateIndex] = 1;
  }

  
  for( int orderingIndex = size-2; orderingIndex >= 0; orderingIndex--) {
    
    int betaIndex = orderingIndex;
    int nextBetaIndex = orderingIndex+1;
    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
      betas[betaIndex][stateIndex] = 0.0;


      for(int innerStateIndex = 0; innerStateIndex < numStates; innerStateIndex++) {


	betas[betaIndex][stateIndex] += transitionsLocal[betaIndex][stateIndex][innerStateIndex]*betas[nextBetaIndex][innerStateIndex]*emissionProbabilities[nextBetaIndex][innerStateIndex];
	if(individualIndex == individualToPrint && expressionSuperVector[orderingIndex]->getSNP()->getLocation() >= minToPrint && expressionSuperVector[orderingIndex]->getSNP()->getLocation() <= maxToPrint) {
	  printf("individual %d For  emit %d innerEmit %d betaIndex %d  and   from loc %d we now have ratio %f exp %f transition %f beta next %f \n", individualIndex, stateIndex, innerStateIndex, betaIndex,   expressionSuperVector[betaIndex]->getSNP()->getLocation(),  expressionSuperVector[nextBetaIndex]->getRatioOfRatios(), expressionSuperVector[nextBetaIndex]->getExpressionLogRatio(),  transitionsLocal[betaIndex][stateIndex][innerStateIndex], betas[nextBetaIndex][innerStateIndex]);
	}
      }



    }

    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
      betas[betaIndex][stateIndex] /= scales[betaIndex];

      if(isnan(betas[betaIndex][stateIndex])) {
	printf("This is a problem in beta %d %d\n", betaIndex, stateIndex);
      }
      else if (betas[betaIndex][stateIndex] == 0.0) {
	printf("we are at EXACTLY 0  %d %d from %f %f now at %f \n", betaIndex, stateIndex, betas[betaIndex+1][stateIndex], scales[betaIndex], betas[betaIndex][stateIndex]);
      }

    }
  }
}



void HMMSNP::viterbiTraceback(double** alphas, int** maxStates, int** viterbiStates, int size, int individualIndex) {
  int maxState = 1;
  float maxVal = -1;

  
  for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
    if(alphas[size-1][stateIndex] > maxVal) {
      maxState = stateIndex;
      maxVal = alphas[size-1][stateIndex];
    }

  }
  

  for(int orderingIndex =size -1; orderingIndex >= 1; orderingIndex--) {

    maxStates[individualIndex][orderingIndex] = maxState;
    
    
    maxState = viterbiStates[orderingIndex][maxState];
   
  }
  maxStates[individualIndex][0] = maxState;
 
}




void HMMSNP::fileOutput(Individual** individuals, ofstream &myViterbiFile,  ofstream &myfile, int chromosomeIndex, vector<ExpressionInfo*> expressionSuperVector, float** iarray, int** maxStates, int individualIndex, char* fileName) {
 
  myfile << "track type=bedGraph name=ind" << individuals[individualIndex]->getName() <<  "chr" << chromosomeIndex << endl;
 
  myViterbiFile << "track type=bedGraph name=ind" <<individuals[individualIndex]->getName()  << "viterbi" <<  "chr" << chromosomeIndex << endl;
  for(int expressionIndex = 0; expressionIndex < expressionSuperVector.size(); expressionIndex++) {
    
    myViterbiFile << "chr" ;
    myViterbiFile << chromosomeIndex << " " ;
    myViterbiFile << expressionSuperVector[expressionIndex]->getSNP()->getLocation() << " " ; 
    myViterbiFile << expressionSuperVector[expressionIndex]->getSNP()->getLocation()+1 <<   " ";
    //    myViterbiFile  << expressionSuperVector[expressionIndex]->getExpressionLogRatio() << " "; 
    //myViterbiFile << expressionSuperVector[expressionIndex]->getRatioOfRatios() << " ";
    myViterbiFile << maxStates[individualIndex][expressionIndex];
    
    myfile << "chr" << chromosomeIndex << " " << expressionSuperVector[expressionIndex]->getSNP()->getLocation() << " " << expressionSuperVector[expressionIndex]->getSNP()->getLocation()+1 <<  " " ;
    float eval = 0;
    for(int stateIndex = 0; stateIndex < numStates; stateIndex++) {
      
      eval += iarray[expressionIndex][stateIndex] * ratioMeans[stateIndex];
    }
  
    
    
    
    myfile << eval << "\n";
    myViterbiFile <<  "\n";
  }

}






void HMMSNP::addToTempTransitions (double**** tempTransitionValues,double** alphas, double**emissionProbabilities, double** betas, double***transitionsLocal, vector<ExpressionInfo*> expressionSuperVector, int chromosomeIndex) {
  
  for( int expressionIndex = 0; expressionIndex < expressionSuperVector.size()-1; expressionIndex++) {
    

    if(expressionSuperVector[expressionIndex+1]->isHeterozygote() || numStates % 2 == 0) {
      double transitionSums[numStates][numStates];
      double transitionFullSums[numStates];
      for(int firstState = 0 ; firstState < numStates; firstState++) {
	transitionFullSums[firstState] = 0.0;
	for(int secondState = 0; secondState < numStates; secondState++) {
	  
	transitionSums[firstState][secondState] = 0.0;
	
	}
	
    }
      
      
      for(int firstState = 0; firstState < numStates; firstState++) {
	
	for(int secondState = 0; secondState < numStates; secondState++) {
	  
	  transitionSums[firstState][secondState] = alphas[expressionIndex][firstState]*transitionsLocal[expressionIndex][firstState][secondState] *emissionProbabilities[expressionIndex+1][secondState]*betas[expressionIndex+1][secondState];

	  


	  transitionFullSums[firstState] += transitionSums[firstState][secondState];
	  tempTransitionValues[chromosomeIndex][expressionIndex][firstState][secondState] += transitionSums[firstState][secondState];
	}
      }
    }
  }
}





void HMMSNP::updateTransitions(double**** tempTransitionValues, float weight, int size, int chromosomeIndex) { 



  for( int expressionIndex = 0; expressionIndex < size; expressionIndex++) {

    double **distanceDependentTransitions = transitionProbabilities;

    double sums [numStates];

    int expressedStates = numStates;
    if(numStates % 2 == 0) {
      expressedStates = numStates-1;
    }

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
	double original  =   ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex];
	double mirrorOriginal =   ltorTransitions[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond];

	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] = tempTransitionValues[chromosomeIndex][expressionIndex][firstIndex][secondIndex] + tempTransitionValues[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond];




	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] += weight * (distanceDependentTransitions[firstIndex][secondIndex] + distanceDependentTransitions[mirrorFirst][mirrorSecond]);

	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] /= (sums[firstIndex] + sums[mirrorFirst] + weight*2);

	ltorTransitions[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond] = ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex];
	



      }


      if(numStates % 2 == 0) {
	int mirrorSecond = numStates-1;
	double original  =   ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][numStates-1];
	double mirrorOriginal =   ltorTransitions[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond];

	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][numStates-1] = tempTransitionValues[chromosomeIndex][expressionIndex][firstIndex][numStates-1] + tempTransitionValues[chromosomeIndex][expressionIndex][mirrorFirst][numStates-1];

	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][numStates-1]   += weight * (distanceDependentTransitions[firstIndex][numStates-1] + distanceDependentTransitions[mirrorFirst][numStates-1]);
	
	
	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][numStates-1] /= (sums[firstIndex] + sums[mirrorFirst] + weight*2);
	ltorTransitions[chromosomeIndex][expressionIndex][mirrorFirst][numStates-1] =  ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][numStates-1];


      }


    }

    if(numStates % 2 == 0) {
      int firstIndex = numStates-1;
      sums[firstIndex] = 0.;
      for(int secondIndex = 0; secondIndex < numStates; secondIndex++) {
	sums[firstIndex] +=   tempTransitionValues[chromosomeIndex][expressionIndex][firstIndex][secondIndex];

      }

      int mirrorFirst = numStates-1;
      for(int secondIndex = 0; secondIndex < numStates ; secondIndex++) {

	int mirrorSecond = expressedStates-secondIndex-1;
	if(secondIndex == numStates-1) {
	  mirrorSecond = secondIndex;

	}
	double original  =   ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex];
	double mirrorOriginal =   ltorTransitions[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond];

	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] = tempTransitionValues[chromosomeIndex][expressionIndex][firstIndex][secondIndex] + tempTransitionValues[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond];




	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] += weight * (distanceDependentTransitions[firstIndex][secondIndex] + distanceDependentTransitions[mirrorFirst][mirrorSecond]);

	ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex] /= (sums[firstIndex] + sums[mirrorFirst] + weight*2);

	ltorTransitions[chromosomeIndex][expressionIndex][mirrorFirst][mirrorSecond] = ltorTransitions[chromosomeIndex][expressionIndex][firstIndex][secondIndex];






      }

    }



    

  }
  
}









void HMMSNP::leftToRight (Individual ** individuals, char* outputFileName, double weight) {
 ofstream myfile;
 ofstream myViterbiFile;
  cout << endl;

  int maxOuter = 100;

  double maxChanged = 1e-02;
  

 


  
  
  
   bool print = 0;

   vector <Chromosome*> tempChromosomes = individuals[0]->getChromosomes();
   int *** minMaxIndices = new int**[tempChromosomes.size()];
   
   for(int i = 0; i < tempChromosomes.size(); i++) {
     minMaxIndices[i] = new int*[90000];
    
     for(int j = 0; j < 90000; j++) {
       minMaxIndices[i][j] = new int[2];
    
     } 
   }




   


  double **** tempTransitionValues = new double ***[tempChromosomes.size()];
  initializeTransitionValues( individuals,  tempTransitionValues) ;
  





  float **iarray = new float*[90000];
  double **alphas = new double*[90000]   ;
  double **betas = new double*[90000];
  int ***viterbiStates = new int**[numIndividuals];
  double **emissionProbabilities = new double*[90000];
  double *scales  = new double[90000];
  for(int individualIndex = 0; individualIndex < numIndividuals; individualIndex++) {
    viterbiStates[individualIndex] = new int*[90000]();
    for(int i = 0; i < 90000; i++) {
      viterbiStates[individualIndex][i] = new int[numStates]();
    }
  }
  

  for(int i = 0; i < 90000; i++) {
    
    iarray[i] = new float[numStates];
    alphas[i] = new double[numStates];
    betas[i] = new double[numStates];
    emissionProbabilities[i] = new double[numStates];
  }


  
  
  
  for(int outer = 0; outer < maxOuter; outer++) {
    cout << "ltor" << outer << endl;
  
    if(outputFileName && (outer  == maxOuter-1 || outer % 99 == 0)) {                                                                                                                                                                                                                                                                                                                           


    }
   

    printf("Left to right loop %d\n", outer);
    bool changed = 0;
    int totalChanged=0;

    for(unsigned int chromosomeIndex = 0; chromosomeIndex < tempChromosomes.size(); chromosomeIndex++) {
      double *** transitionsLocal = ltorTransitions[chromosomeIndex];
      int sizeOfArray = individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1).size();      
      if(sizeOfArray < 1) {
        continue;
      }
      
      int**  maxStates = new int* [numIndividuals];
    
    
      
      perChromosomeSet(individuals, maxStates, tempTransitionValues, iarray, chromosomeIndex);
 

      for(int individualIndex = 0; individualIndex < numIndividuals; individualIndex++) {
	
	vector<Chromosome*> chromosomes = individuals[individualIndex]->getChromosomes();
        vector <ExpressionInfo*> expressionSuperVector= chromosomes[chromosomeIndex]->getFullExpressions(1);
	
	
	for( int index = 0; index < expressionSuperVector.size(); index++) {
	  
	  scales[index] = 0.0;
	  for( int innerIndex = 0; innerIndex < numStates; innerIndex++) {
	    betas[index][innerIndex] = 0.0;
	    alphas[index][innerIndex] = 0.0;
	    emissionProbabilities[index][innerIndex] = 0.0;
	    viterbiStates[individualIndex][index][innerIndex] = 0;
	  }
	  
	}
     
	
	  setEmissionProbabilities(   individuals, emissionProbabilities,  expressionSuperVector,   individualIndex, chromosomeIndex,  runifValue) ;	
	
              
	
	  forward(individuals,  alphas, emissionProbabilities,  scales,  transitionsLocal, viterbiStates[individualIndex], expressionSuperVector.size(), individualIndex, expressionSuperVector) ;
	

	  backward( individuals,  betas,  scales,  transitionsLocal,  emissionProbabilities,  expressionSuperVector.size(), individualIndex, expressionSuperVector);
	  
	  
	  viterbiTraceback( alphas,  maxStates,  viterbiStates[individualIndex], expressionSuperVector.size(), individualIndex) ;
	
       
	  
	 
	    
	  
      

	  if(outputFileName && (outer  == maxOuter-1 || outer % 99 == 0)) {	  

	    
	    if(individualIndex == 0) {
	      char viterbiName[249];
	      char evalName[249];
	      sprintf(evalName, "%schr%diter%d.evals",  outputFileName, chromosomeIndex, outer);
	      myfile.open(evalName);
	      sprintf(viterbiName, "%schr%diter%d.viterbi",  outputFileName, chromosomeIndex, outer);
	      myViterbiFile.open(viterbiName);
	    }
	    fileOutput(individuals, myViterbiFile,  myfile,  chromosomeIndex,  expressionSuperVector, iarray, maxStates, individualIndex, outputFileName) ;
	    if(individualIndex == numIndividuals-1) {
	      myfile.close();
	      myViterbiFile.close();
	    }
	  }
	
	
      

      
	
		

	addToTempTransitions ( tempTransitionValues,alphas, emissionProbabilities,  betas, transitionsLocal,  expressionSuperVector, chromosomeIndex) ;

      }
    

      updateTransitions(tempTransitionValues, weight, individuals[0]->getChromosomes()[chromosomeIndex]->getFullExpressions(1).size(), chromosomeIndex);      
    
      for(int indIndex = 0 ; indIndex < numIndividuals; indIndex++) {
	delete[] maxStates[indIndex];
      }
      delete[] maxStates;
    }
    

    
  }
  
  

  for(int i = 0; i < 90000; i++) {
    
    delete[]iarray[i];
    delete[] alphas[i];
    for(int individualIndex = 0; individualIndex < numIndividuals; individualIndex++)
      delete[] viterbiStates[individualIndex][i];
    delete[] betas[i];
    delete[] emissionProbabilities[i];
  }
  for(int individualIndex = 0; individualIndex < numIndividuals; individualIndex++)
    delete[] viterbiStates[individualIndex];
  
  delete []iarray;
  delete []alphas;
  delete[]scales;
  delete[]viterbiStates;
  delete[]betas;
  delete [] emissionProbabilities;
  
  for(unsigned int chromosomeIter = 0; chromosomeIter < tempChromosomes.size(); chromosomeIter++) {
    
   int size = individuals[0]->getChromosomes()[chromosomeIter]->getFullExpressions(1).size();
   
   for(int SNPiter = 0; SNPiter < size; SNPiter++) {
     
     
     for(int k = 0; k < numStates; k++) {
      
       delete []tempTransitionValues[chromosomeIter][SNPiter][k];
       for(int l = 0 ; l < numStates; l++) {
	 
	    
       }
     }
     
     delete[]tempTransitionValues[chromosomeIter][SNPiter];
     
   }
   
   delete []tempTransitionValues[chromosomeIter];
    
 } 
 
 delete[]tempTransitionValues;
 deleteLTOR(individuals);
}


    


void HMMSNP::deleteLTOR(Individual **individuals) {
  for(unsigned int chromosomeIter = 0; chromosomeIter < individuals[0]->getChromosomes().size(); chromosomeIter++) {
	
    int size = individuals[0]->getChromosomes()[chromosomeIter]->getFullExpressions(1).size();
	
    for(int SNPiter = 0; SNPiter < size; SNPiter++) {
     
    
      for(int k = 0; k < numStates; k++) {
	delete []ltorTransitions[chromosomeIter][SNPiter][k];
      
      
          
        
      }
      delete []ltorTransitions[chromosomeIter][SNPiter];
      
      
    }
    delete []ltorTransitions[chromosomeIter];
    
    
  } 
  delete[]ltorTransitions;
  ltorTransitions = 0;
  
}



