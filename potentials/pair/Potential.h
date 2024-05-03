#pragma once

#include <string>

struct Potential
{
public:
    int
        potentialType;                   // potential type

    std::string
        iona,                    // element for atom a
        ionb;                    // element for atom b

    double
        param1,                  // short- range interaction parameters
        param2,
        param3,
        param4,
        param5,
        param6;                    // maximum distance for potential

    virtual ~Potential()
    {
    }

    Potential()
    {
        param1 = 0.0;                  // short- range interaction parameters are set to zero
        param2 = 0.0;
        param3 = 0.0;
        param4 = 0.0;
        param5 = 0.0;
        param6 = 0.0;

    }


};