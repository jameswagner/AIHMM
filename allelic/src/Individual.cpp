#include "Individual.h"

Individual::Individual() : name("") {}

Individual::Individual(const std::string& newname) : name(newname) {}

std::string Individual::getName() {
    return name;
}

void Individual::setName(const std::string& newname) {
    name = newname;
}

Individual::~Individual() {
    for (auto chromosome : chromosomes) {
        delete chromosome;
    }
    chromosomes.clear();
}

#include <iostream>

void Individual::setExpression(unsigned int chromosomeIndex, const std::vector<ExpressionInfo*>& expressionInfo) {
    if (chromosomeIndex >= 0 && chromosomeIndex < chromosomes.size()) {
        chromosomes[chromosomeIndex]->setExpression(expressionInfo);
    } else {
        std::cerr << "Invalid chromosome index." << std::endl;
    }
}

std::vector<Chromosome*>& Individual::getChromosomes() {
    return chromosomes;
}

void Individual::addChromosome(Chromosome* chromosome) {
    chromosomes.push_back(chromosome);
}
