#ifndef FILEREADER_H_
#define FILEREADER_H_

#include "Individual.h"
#include <sstream>
#include <string>
#include <vector>
#include "SNP.h"
#include <fstream>
#include <iostream>
#include "ExpressionInfo.h"
#include "constants.h"

class FileReader {
public:
    FileReader();
    virtual ~FileReader();

    void readFiles(std::vector<Individual*>& individuals, int startChromosome, int endChromosome, const std::string& fileNamePrefix);

};

#endif /* FILEREADER_H_ */
