#pragma once


#include <fstream>
#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>

#include "MPICommunicator.h"
#include "Constants.h"

struct RelaxCntrl
{
    int
        maxIter = 3000,
        print = 1;

    double
        normTol = 1.0e-5,
        forceTol = 1.0e-4;
        //energyTol = 1.0e-7;

    
    //FIRE defaults
    int
        N_delaySteps = 20,
        N_maxUphill = 50,
        N_min = 5;

    double
        initTimeStep = 0.0005,
        alphaDecrease = 0.99,
        stepIncrease = 1.1,
        stepDecrease = 0.5,
        alphaStart = 0.1,
        maxTimeStep = 0.00075,
        minTimeStep = 0.0001;

    bool
        initialDelay = true,
        debug = false;

    std::string
        method = "fire";

    bool
        verletIntegration = false;

    void readRelax(std::ifstream& inStream, std::ofstream& outStream);

    RelaxCntrl() {};
    virtual ~RelaxCntrl() {};

    RelaxCntrl& operator=(const RelaxCntrl& src);
};
