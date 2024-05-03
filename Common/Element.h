#pragma once

#include <string>

struct Element
{

    int
        /** an integer type for the atom species */
        inttype = -1,
	    
	/** atomic number */
	atomicNumber = -1;

    std::string
        /** Element name */
        name;

        /** type ie core/shell */
        //type;

    double
        /** mass of this element */
        mass,

        /** atomic charge of this element */
        charge;

    Element() {};
    virtual ~Element() {};

};
