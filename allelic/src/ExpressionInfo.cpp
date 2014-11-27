/*
 * ExpressionInfo.cpp
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#include "ExpressionInfo.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
using namespace std;
ExpressionInfo::ExpressionInfo() {

}


ExpressionInfo::~ExpressionInfo() {

  delete snp;
  

}
void ExpressionInfo::setRatioOfRatios(float ratio) {
	ratioOfRatios = ratio;
}


ExpressionInfo::ExpressionInfo(SNP* new_snp, bool new_heterozygote, float new_ratioOfRatios, float new_expressionRatio) {

	snp = new_snp;
	ratioOfRatios = new_ratioOfRatios;


	if(ratioOfRatios < -5) {
		ratioOfRatios = -5;
	}
	else if(ratioOfRatios > 5) {
		ratioOfRatios = 5;
	}
	expressionRatio  = new_expressionRatio;
	heterozygote = new_heterozygote;


}

bool ExpressionInfo::isHeterozygote() {
  return heterozygote;
}



void ExpressionInfo::setHeterozygote(bool isHeterozygote) {
  heterozygote = isHeterozygote;
}

SNP* ExpressionInfo::getSNP() {
  return snp;
}

float ExpressionInfo::getExpressionRatio() {
  
  return expressionRatio;
}
void ExpressionInfo::setExpressionRatio(float new_ratio) {
  
  expressionRatio = new_ratio;
}

float ExpressionInfo::getExpressionLogRatio() {


  if(expressionRatio <= -10) {
    expressionRatio = -10;
    
  }
  else if(expressionRatio >= 4) {
    expressionRatio = 1;
  }
  return expressionRatio;
}

float ExpressionInfo::getRatioOfRatios() {
  return ratioOfRatios;
}
