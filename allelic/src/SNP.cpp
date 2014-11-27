/*
 * SNP.cpp
 *
 *  Created on: Dec 3, 2008
 *      Author: jameswagner
 */

#include "SNP.h"
#include <stdio.h>
#include <string.h>
SNP::SNP() {


}


SNP::SNP(char* new_name, char* new_chromosome, int new_location) {
  location = new_location;
  strcpy(name, new_name);
  
  strcpy(chromosome,  new_chromosome);
	



}







SNP::~SNP() {
  delete name;
  delete chromosome;

}

int SNP::getLocation() {
	return location;
}

char* SNP::getName() {
	return name;
}
char* SNP::getChromosome() {
  return chromosome;
}
