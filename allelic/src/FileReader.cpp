/*
 * FileReader.cpp
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#include "FileReader.h"
#include "SNP.h"
#include <string.h>
#include <fstream>
#include <iostream>
#include "ExpressionInfo.h"
#include "constants.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
FileReader::FileReader() {
	

}

FileReader::~FileReader() {
  
}




void  FileReader::readFiles(Individual **individuals,  int startChromosome, int endChromosome,  char* fileNamePrefix) {

  
  for(int individualIter = 0; individualIter < numIndividuals; individualIter++) {
    individuals[individualIter] = new Individual();
  }
    
  char line[100000]; 
   for(int indIter = 0; indIter < numIndividuals; indIter++) {
    for(int chromosomeIndex = 1; chromosomeIndex <= 23; chromosomeIndex++) {
      Chromosome *chromosome = new Chromosome;
      individuals[indIter]->addChromosome(chromosome);
    }
   }
   
  
   for(int chromosomeIter = startChromosome ; chromosomeIter <= endChromosome; chromosomeIter++) {
     
     char fileName[500];
     sprintf(fileName, "%s%d.txt", fileNamePrefix, chromosomeIter);
     ifstream infile;
    cout << endl;
    infile.open(fileName, ifstream::in);
    
    infile.getline(line, 1000);
    
    //sample id stuff;
    char* name = new char[40];
    strtok(line, " \t");
    
    for(int i = 0; i < numIndividuals; i++) {
      strcpy(name, strtok(NULL, " \t"));
      individuals[i]->setName(name);
    }

    while(infile.good()) {

      infile.getline(line, 100000);
      
      if(strlen(line) < 10) {
	continue;
	
     }
      
      char* token;
      
      char *name = new char[25];
      char *chromosome = new char[3];
      char *location = new char[25];
      
      strcpy(name, strtok(line, " \t"));
      strcpy(chromosome, strtok(NULL, " \t"));
      strcpy(location, strtok(NULL, " \t"));
      
      SNP *snp = new SNP(name, chromosome, atoi(location));

      char* tokens[numIndividuals*3];
      int index = 0;
      token = strtok(NULL, " \t");
      while(token != NULL && index < numIndividuals*3) {
	
       tokens[index++] = token;
       
       token = strtok(NULL, " \t");
                                                                              
     }
     
     for(int individualIndex = 0; individualIndex < numIndividuals; individualIndex++) {
 
       bool isHet = atoi(tokens[individualIndex*3]);
       float expressionLogRatio = atof(tokens[individualIndex*3+1]) ;;
       float ratioOfRatios = atof(tokens[individualIndex*3+2]);
       
       ExpressionInfo *expressionInfo  = new ExpressionInfo(snp, isHet, ratioOfRatios, expressionLogRatio);
       
       individuals[individualIndex]->getChromosomes()[atoi(chromosome)]->addExpression(expressionInfo);
       
     }
     delete []name;
     delete [] chromosome;
     delete [] location;
    }
  

    
   }



    for(int chromosomeIter = 1; chromosomeIter <=22; chromosomeIter++)  
      cout << "This is the size " << chromosomeIter << " " << individuals[0]->getChromosomes()[chromosomeIter]->getFullExpressions(1).size() << endl;
    
}

