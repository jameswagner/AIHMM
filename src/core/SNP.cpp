#include "SNP.h"


SNP::SNP(const std::string& new_name, const std::string& new_chromosome, int new_location,  const std::string& new_GeneName)
  : name(new_name), chromosome(new_chromosome), location(new_location), geneName(new_GeneName) {}


SNP::~SNP() {}

int SNP::getLocation() {
    return location;
}

std::string SNP::getName() {
    return name;
}

std::string SNP::getChromosome() {
    return chromosome;
}
