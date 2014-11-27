/*
 * Individual.cpp
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#include "Individual.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <vector>


#include <string.h>

#include <stdio.h>

using namespace std;
Individual::Individual() {

  name = new char[40];
}

Individual::Individual(char* newname) {
  name = new char[40];
  strcpy(name, newname);
}
char* Individual::getName() {
  return name;
}
void Individual::setName(const char* newname) {
  strcpy(name, newname);
}
Individual::~Individual() {

  chromosomes.clear(); 
}
void Individual::setExpression(int chromosomeIndex, vector <ExpressionInfo*> expressionInfo) {

  chromosomes[chromosomeIndex]->setExpression(expressionInfo);

}

vector <Chromosome*> Individual::getChromosomes() {

  
	return chromosomes;
}

void Individual::addChromosome(Chromosome *chromosome) {
  chromosomes.push_back( chromosome);
}

