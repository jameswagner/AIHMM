#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include "Chromosome.h"
#include <vector>
#include <string>
#include <iostream>

class Individual {
public:
    Individual();
    Individual(const std::string& name);
    ~Individual();

    std::vector<Chromosome*>& getChromosomes();
    void setExpression(unsigned int chromosomeIndex, const std::vector<ExpressionInfo*>& ei);
    void addChromosome(Chromosome* chromosome);

    void setName(const std::string& newName);
    std::string getName();

private:
    std::string name;
    std::vector<Chromosome*> chromosomes;
};

#endif /* INDIVIDUAL_H_ */