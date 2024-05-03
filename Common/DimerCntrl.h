#pragma once

#include <fstream>
#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>

#include "MPICommunicator.h"
#include "Constants.h"
#include "RelaxCntrl.h"
#include "ArtCntrl.h"



struct DimerCntrl
{
    double
        /* the angle used for numerical derivatives (in radians) */
        deltaTheta = 0.005,

        /* the size of the displacement for the dimer */
        dimerDisplacement = 0.01,

        /* the size of the displacement towards the saddle point for the line search */
        lineDisplacement = 0.001,

        /* the maximum displacement in line search */
        maxDeltaX = 0.5;

    bool
        /* turns on extra printing for diagnostic purposes */
        debug = true;


    int
        /* the maximum number of rotations in Dimer optimisation */
        maxRotations = 20,
        
        /* the type of activation employed. 1 = global(default), 2 = local, 3 = globaltypes */
        regionStyle = 1;

    double
        /* the radius around an atom that is selected */
        regionRadius = 5.0;

    std::vector<std::string>
        /** vector of names of atom types to be targeted in type selection/displacement */
        dimerTypes1,
        dimerTypes2;

    std::vector<int>
        /** vector containing the ideal coordination number of an ion type */
        idealCoordNumber;

    std::vector<double>
        /** vector containing the distance used to determine the coordination number */
        coordDistance,

        /** vector holding minimum position to mask atom movement along x vector*/
        maskLow,

        /** vector holding minimum position to mask atom movement */
        maskHigh;

    //below are the parameters for the parameters used within the Dimer minimisation
    int
        maxParIter = 10,
        print = 1;

    double
        /* the maximum allowed eigenvalue */
        maxEigenvalue = 0.0,
        /* the rotation force tolerance to check whether dimer minimum mode step has converged */
        normTol = 0.1,
        forceTol = 1.0e-2,
        energyTol = 1.0e-4;


    //Dimer defaults
    int
        N_min = 5;

    double
        alphaDecrease = 0.99,
        stepIncrease = 1.1,
        stepDecrease = 0.5,
        alphaStart = 0.1;

    void readDimer(std::ifstream& inStream, std::ofstream& outStream);

    DimerCntrl() {};
    virtual ~DimerCntrl() {};
};

