/*
 * Chromosome.h
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_
#include <vector>
#include "ExpressionInfo.h"
using namespace std;
class Chromosome {
public:
	Chromosome();
	vector <ExpressionInfo*> getFullExpressions(bool useLow) ;
	virtual ~Chromosome();
	void addExpression(ExpressionInfo *ei);
	void setExpression(vector <ExpressionInfo*> ei);
	
private:
	vector <ExpressionInfo*> expressions;
	vector <ExpressionInfo*> heterozygotes;
};

#endif /* CHROMOSOME_H_ */
