/*
 * Chromosome.cpp
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#include "Chromosome.h"
#include <stdlib.h>
#include <stdio.h>


#include <stdlib.h>
#include <math.h>

#include <vector>


#include <string.h>
#include <iostream>


using namespace std;
Chromosome::Chromosome() {
  

	

}

//Used to change the vector of ExpressionInfo's corresponding to the Chromosome
void Chromosome::setExpression(vector <ExpressionInfo*> ei) {
  expressions.clear();
  heterozygotes.clear();
  expressions = ei;
 
}


Chromosome::~Chromosome() {

  
  
  heterozygotes.resize(0);
  expressions.clear();  
}





//Call with useLow=1 to return all heterozygous SNPs in the chromosome, 0 to only include those whose total RNA 
//expression value (both alleles) is >= 1000
vector <ExpressionInfo*> Chromosome::getFullExpressions(bool useLow) {
  if(useLow) {
    return expressions;
  }
  else {
    //    return getHighOnlyExpressions();
  }
}









void Chromosome::addExpression(ExpressionInfo* ei) {
  
	expressions.push_back( ei);

}






