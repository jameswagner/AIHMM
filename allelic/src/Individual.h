/*
 * Individual.h
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_
#include "Chromosome.h"

#include <vector>
class Individual {
public:
	Individual();
	Individual(char* name);
	virtual ~Individual();
	vector <Chromosome*> getChromosomes();
	void setExpression(int chromosomeIndex, vector <ExpressionInfo*> ei);
	void addChromosome(Chromosome *chromosome);

	void setName(const char* newName);
	char* getName();


private:

	char* name;

	vector <Chromosome*> chromosomes;

};

#endif /* INDIVIDUAL_H_ */
