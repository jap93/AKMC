// species.cpp: implementation of the Species class.
//
//////////////////////////////////////////////////////////////////////

#include "Species.h"

std::vector<std::string> split(std::string s);

//////////////////////////////////////////////////////////////////////
// construction/destruction
//////////////////////////////////////////////////////////////////////
/* default constructor */
Species::Species()
{
    numberOfElements = 0;
}

/* default destructor */
Species::~Species()
{

}

/** dimensions and loads elemental data from file
 * m_element, m_mass, m_charge assumed to be on the same line
 * receives the number of different elements present and
 * the file to read from.
 * returns integer for error checking
 */

void Species::loadSpecies(std::ifstream& inStream, int num)
{
    int
        i;

    std::string
        dummy,
        line;

    Element
        ele;

    std::vector<std::string> words;

    numberOfElements = num;

    for (i = 0; i < numberOfElements; i++)
    {

        std::getline(inStream, line);
        words = split(line);

        ele.name = words[0];

        dummy = words[1];
        ele.mass = std::stod(dummy);

        dummy = words[2];
        ele.charge = std::stod(dummy);

	    if (words.size() > 3)
            dummy = words[3];
            ele.atomicNumber = std::stoi(dummy);

        eleData.push_back(ele);
        // error checking

    }

}


/**
  function to print species data
 **/

void Species::printSpecies(std::ofstream& outstream)
{
    Element
        ele;

    outstream << "\n \n *********************************************************************" << std::endl;
    outstream << "                         atom data ";
    outstream << "\n \n *********************************************************************" << std::endl;

    outstream << "\n        element         charge           mass" << std::endl;

    for (int i = 0; i < numberOfElements; i++)
    {
        ele = eleData[i];

        outstream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
            << std::setprecision(3) << "\n"
            << std::setw(15) << ele.name
            << std::setw(15) << ele.charge
            << std::setw(15) << ele.mass << std::endl;
    }

    outstream.flush();

}


/************************************************************
 returns the data for element i
 ***********************************************************/
Element Species::getSpecies(int i) const
{
    Element
        dummy;

    dummy = eleData[i];
    return dummy;
}



/************************************************************
 returns the mass of an element
 ***********************************************************/
double Species::getMass(int i)
{
    return eleData[i].mass;
}
