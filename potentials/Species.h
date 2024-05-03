#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <algorithm>

#include <vector>

#include "Element.h"

class Species
{

private:
    int
	    /** the number of elements */
        numberOfElements; 

    
    std::vector<Element> 
	    /** vectors of elements containing names, mass and charge */
		eleData;         

public:
    // default constructor
    Species();
    // destructor
    virtual ~Species();

		/**
		 *  read elements from file
		 * @param inStream (ifstream&) : the input stream (already open
		 * @param numPotentials (int) : the expected number of potentials
		 */
		void loadSpecies(std::ifstream& inStream, int numPotentials);

		/**
		 *  print element data
		 * @param outStream (ofstream&) : output stream
		 */
		void printSpecies(std::ofstream& outStream);

		/**
		 * returns the data of an element
		 * @param i (int) : the number of the element to be returned
		 * @return (Element) : element i
		 */
		Element getSpecies(int i) const;

		/**
		 * returns the number of different elements
		 * @return int : number of elements
		 */
		int getNumSpecies(void) const
		{
			return numberOfElements;
		}

		// returns the maxximum number of different interactions between the elements - not used
		int getMaxNumInteractions(void) const
		{
			return numberOfElements * numberOfElements; //return (numberOfElements * (numberOfElements + 1)) / 2;
		}

		double getMass(int);  // returns mass


};

