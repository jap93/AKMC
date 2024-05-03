#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "Species.h"
#include "Constants.h"
#include "Memory.h"




struct MetalPotential
{
    int
        potentialType = -1,
        numComponents = 0,
        numPtsRep = 0,
        numPtsDens = 0,
        numPtsEmbed = 0;

    double
        startRep = 0.0,
        finishRep = 0.0,
        startDens = 0.0,
        finishDens = 0.0,
        startEmbed = 0.0,
        finishEmbed = 0.0,
        stepEmbed = 0.0,
        stepDens = 0.0,
        stepRep = 0.0,
        potCutoff = 0.0;

    double
        param1 = 0.0,                  // short- range interaction parameters for parameterised potentials
        param2 = 0.0,
        param3 = 0.0,
        param4 = 0.0,
        param5 = 0.0,
        param6 = 0.0,
        param7 = 0.0, 
        param8 = 0.0,
        param9 = 0.0;

    double
        *parRep = nullptr,
        *parDens = nullptr,
        *parEmbed = nullptr,
        *derivRep = nullptr,
        *derivDens = nullptr,
        *derivEmbed = nullptr,
        *distRep = nullptr,
        *distDens = nullptr,
        *distEmbed = nullptr;         // max distance

    std::string
        atA,
        atB;

    void readEAM(std::string infile, std::ofstream& outStream);

    void readLEAM(std::string infile, std::ofstream& outStream);

    void readParameters(double* par, int numPts, std::ifstream& inStream);

    void calculateDerivativeMesh(double* pot, double* deriv, double dr, int npts);

    MetalPotential& operator=(const MetalPotential& src);

    MetalPotential(const MetalPotential& src);

    MetalPotential();

    virtual ~MetalPotential();
    
};

