#ifndef SNP_H_
#define SNP_H_

#include <string>

class SNP {
public:
    SNP(const std::string& new_name, const std::string& new_chromosome, int new_location, const std::string& new_GeneName = "");
    void setGeneName(const std::string& geneName);
    virtual ~SNP();
    std::string getName();
    std::string getGeneName();
    bool getGene();
    int getLocation();
    std::string getChromosome();

private:
    std::string name;
    std::string chromosome;
    int location;
    std::string geneName;
};

#endif /* SNP_H_ */
