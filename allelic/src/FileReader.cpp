#include "FileReader.h"
#include <cstring>



FileReader::FileReader() {}

FileReader::~FileReader() {}

  
void FileReader::readFiles(std::vector<Individual*>& individuals, int startChromosome, int endChromosome, const std::string& fileNamePrefix) {
    
    for (unsigned int individualIter = 0; individualIter < individuals.size(); individualIter++) {
        individuals[individualIter] = new Individual();
    }

    char line[100000];
    for (unsigned int individualIter = 0; individualIter < individuals.size(); individualIter++) {
        for (int chromosomeIndex = startChromosome; chromosomeIndex <= endChromosome; chromosomeIndex++) {
            Chromosome* chromosome = new Chromosome;
            individuals[individualIter]->addChromosome(chromosome);
        }
    }
   
  
    for (int chromosomeIter = startChromosome; chromosomeIter <= endChromosome; chromosomeIter++) {
        std::string fileName = fileNamePrefix + std::to_string(chromosomeIter) + ".txt";
        std::ifstream infile(fileName);
        
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            continue;
        }

        infile.getline(line, sizeof(line));
        
        // Sample id stuff
        std::vector<std::string> sampleNames;
        std::string sampleNamesLine(line);
        std::istringstream ssSampleNames(sampleNamesLine);
        std::string sampleName;

        while (std::getline(ssSampleNames, sampleName, '\t')) {
            sampleNames.push_back(sampleName);
        }

        for (size_t i = 0; i < sampleNames.size() && i < individuals.size(); ++i) {
            individuals[i]->setName(sampleNames[i]);
        }

        while (infile.getline(line, sizeof(line))) {
            if (strlen(line) < 10) {
                continue;
            }
      
            std::istringstream ss(line);
            std::string name, chromosome, location;
            std::getline(ss, name, '\t');
            std::getline(ss, chromosome, '\t');
            std::getline(ss, location, '\t');
          
            SNP* snp = new SNP(name, chromosome, std::stoi(location));
            std::vector<std::string> tokens;
            std::string token;
            while (std::getline(ss, token, '\t')) {
                tokens.push_back(token);
            }
        
            for (unsigned int i = 0; i < tokens.size(); i += 3) {
                bool isHet = std::stoi(tokens[i]);
                float expressionLogRatio = std::stof(tokens[i + 1]);
                float ratioOfRatios = std::stof(tokens[i + 2]);
                ExpressionInfo* expressionInfo = new ExpressionInfo(snp, isHet, ratioOfRatios, expressionLogRatio);
                individuals[i / 3]->getChromosomes()[std::stoi(chromosome)-1]->addExpression(expressionInfo);
            }
        }
   }

    for(int chromosomeIter = 1; chromosomeIter <=endChromosome; chromosomeIter++)  
      std::cout << "This is the size of chromosome" << chromosomeIter << " " << individuals[0]->getChromosomes()[chromosomeIter]->getFullExpressions(1).size() << std::endl;
    
}

