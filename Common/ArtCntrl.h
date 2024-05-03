#pragma once

#include <fstream>
#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>

#include "MPICommunicator.h"
#include "Constants.h"

struct ARTCntrl
{
    int
        /* the number of Lanczos vectors in j dimensional basis */
        numLanczosVectors = 20,
        
        /* the number of allowed steps in the activation phase */
        maxActivationSteps = 10;

    double
        /* the step size used for numerical derivatives */
        deltaX = 0.01,

        /* the maximum displacement in line search */
        maxDeltaX = 0.2,

        /* the size of the initial displacement for the move out of the basin */
        initialDisplacement = 0.01,

        /* the maximum allowed eigen value in activation stage */
        maxEigenvalue = 0.0,
        
        /* the tolereance to check whether Lanczos step has convereged */
        eigValTol = 0.001,
        
        /* damping factor for the parallel force */
        parDamp = 0.5;

    bool
        /* flag indicating that the perpendicular force should be minimised first */
        goPeprepedicular = false,

        /* turns on extra printing for diagnostic purposes */
        debug = true;

    //below are the parameters for the parameters used within the FIRE minimisation
    int
        maxParIter = 1000,
        maxPerpIter = 20,
        print = 1;

    double
        normTol = 1.0e-4,
        forceTol = 1.0e-3,
        //energyTol = 1.0e-7,
        maxStep = 0.01;

    std::string
        minMethod = "fire";

    //FIRE defaults
    int
        N_min = 5;

    double
        initTimeStep = 0.005,
        alphaDecrease = 0.99,
        stepIncrease = 1.1,
        stepDecrease = 0.5,
        alphaStart = 0.1,
        maxTimeStep = 0.01,
        initTimeStepPerp = 0.005,
        maxTimeStepPerp = 0.01;

    /**
     * @brief writes out the details from the ART input and checks consistency
     * 
     * @param outStream (std::ofstream&) : the output stream
     * @return numErrors (int) : the number of errors in the input
     */
    int writeARTDetails(std::ofstream& outStream);

    void readART(std::ifstream& inStream, std::ofstream& outStream);

    ARTCntrl() {};
    virtual ~ARTCntrl() {};
};

