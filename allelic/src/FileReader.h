/*
 * FileReader.h
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#ifndef FILEREADER_H_
#define FILEREADER_H_
#include "Individual.h"
#include <string>
class FileReader {
public:
	FileReader();
	virtual ~FileReader();
	 
	void readFiles(Individual **individuals,  int startChromosome, int endChromosome,  char* fileNamePrefix );
	 

};

#endif /* FILEREADER_H_ */
