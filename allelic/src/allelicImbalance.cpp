///============================================================================
// Name        : allelicImbalance.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : 
//============================================================================

#include <iostream>
#include "constants.h"
#include "FileReader.h"
#include "Individual.h"
#include "Chromosome.h"
#include "ExpressionInfo.h"

#include "HMMSNP.h"
#include <string.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <vector>
#include <string.h>
#include <time.h>


using namespace std;



extern int numChromosomes;
extern int numChromosomes2;




void bwLeftToRight(Individual** individuals,  char* startParametersFile, char* outputFile) {
  HMMSNP hmm(startParametersFile);
  //    hmm.baumWelch(individuals, outputFile);
    char* outputName = new char[250]();


    if(outputFile) {
      strcpy(outputName, outputFile);
      strcat(outputName, "LTOR");
    }
    else
      outputName = 0;
    hmm.leftToRight(individuals, outputName, 5.); 

}


void parseCommandLineArguments(int argc, char** argv) {

Individual* permutedIndividuals[numIndividuals];
int** foldNumbers;
 FileReader fileReader;
  map <string, char*> arguments;
  arguments["numChromosomes"] = "1";
  string string1;
  arguments["outputFile"] = 0;
  arguments["numIndividuals"] = "55";

  arguments["fileNamePrefix"] = "/scratch/james/ceu/ceuRevisedRatiosChrom";
  for(int argIndex = 1; argIndex < argc-1; argIndex+=2) {
   
    char* flag = argv[argIndex];
    if(flag[0] != '-') {
      cout <<"Problem with argument " << argv[argIndex] << endl;
      
    }
      printf("Inserting %s %s\n", argv[argIndex]+1, argv[argIndex+1]);
    
    string1 = argv[argIndex]+1;
    arguments[string1] = argv[argIndex+1];
  }







  
  numIndividuals = atoi(arguments["numIndividuals"]);
  
  Individual* individuals[numIndividuals];
  

  
int  numChromosomes = atoi(arguments["numChromosomes"]);
  
  char*  fileNamePrefix = arguments["fileNamePrefix"]; 
  
  fileReader.readFiles(individuals, 1, numChromosomes,  fileNamePrefix);


  bwLeftToRight( individuals,  arguments["startFile"], arguments["outputFile"]);
 

}




int main(int argc, char *argv[]) {
 
   
  parseCommandLineArguments(argc, argv);
 



  
}

