///============================================================================
// Name        : allelicImbalance.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : 
//============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <boost/program_options.hpp>
#include "constants.h"
#include "FileReader.h"
#include "Individual.h"
#include "Chromosome.h"
#include "ExpressionInfo.h"
#include "HMMSNP.h"

namespace po = boost::program_options;

void bwLeftToRight(std::vector<std::unique_ptr<Individual>>& individuals, const std::string& startParametersFile, const std::string& outputFile) {
    HMMSNP hmm(startParametersFile);
    std::string outputName;
    if (!outputFile.empty()) {
        outputName = outputFile + "BW";
    }
    std::cout << "Output file: " << outputFile << std::endl;
    hmm.baumWelch(individuals, outputName);
    //hmm.leftToRight(individuals, outputName, 5.0);
}


void parseCommandLineArguments(int argc, char** argv, int& numIndividuals, int& numChromosomes, std::string& fileNamePrefix, std::string& outputFile, std::string& startFile) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("outputFile,o", po::value<std::string>(), "Set the output file name")
        ("numIndividuals,n", po::value<int>()->default_value(55), "Set the number of individuals")
        ("numChromosomes,c", po::value<int>()->default_value(1), "Set the number of chromosomes")
        ("fileNamePrefix,p", po::value<std::string>()->default_value("./ceu/ceuRevisedRatiosChrom"), "Set the file name prefix")
        ("startFile,s", po::value<std::string>(), "Set the start file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return;
    }

    if (!vm.count("startFile")) {
        throw std::runtime_error("Error: startFile is required.");
    }

    numIndividuals = vm["numIndividuals"].as<int>();
    numChromosomes = vm["numChromosomes"].as<int>();
    fileNamePrefix = vm["fileNamePrefix"].as<std::string>();
    startFile = vm["startFile"].as<std::string>();

    if (vm.count("outputFile")) {
        outputFile = vm["outputFile"].as<std::string>();
    }
}

int main(int argc, char *argv[]) {
    try {
        int numIndividuals = 0;
        int numChromosomes = 0;
        std::string fileNamePrefix;
        std::string outputFile;
        std::string startFile;

        parseCommandLineArguments(argc, argv, numIndividuals, numChromosomes, fileNamePrefix, outputFile, startFile);

        std::vector<std::unique_ptr<Individual>> individuals(numIndividuals);

        FileReader fileReader;
        fileReader.readFiles(individuals, 1, numChromosomes, fileNamePrefix);

        bwLeftToRight(individuals, startFile, outputFile);


    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
