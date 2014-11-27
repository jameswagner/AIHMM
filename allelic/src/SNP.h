/*
 * SNP.h
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#ifndef SNP_H_
#define SNP_H_

class SNP {
public:
	SNP();
	SNP(char* new_name, char* new_chromosome, int new_location, char* new_type, char* new_GeneName, int new_exon);
	SNP(char* new_name, char* new_chromosome, int new_location);
	void setGeneName(char* geneName);
	virtual ~SNP();
	char* getName();
	char* getGeneName();
	bool getGene();
	
	int getLocation();
	char* getChromosome();
	
private:
	int location;
	char chromosome [3];
	char name[25];
	
};

#endif /* SNP_H_ */

