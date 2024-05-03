#pragma once

#include <fstream>
#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm>

#include "MPICommunicator.h"
#include "Constants.h"




/**
 * Control parameters for the nudged elastic band method - also includes
 * FIRE relaxation parameters specific to NEB calcultion
 */
struct NudgeCntrl
{
    int
        numBasins = 1,
        maxIter = 3000,
        print = 10;

    double
        normTol = 1.0e-4,
        forceTol = 1.0e-4,
        energyTol = 1.0e-5;


    //FIRE defaults
    int
        N_min = 5;

    double
        initTimeStep = 0.0005,
        alphaDecrease = 0.99,
        stepIncrease = 1.1,
        stepDecrease = 0.5,
        alphaStart = 0.1,
        maxTimeStep = 0.001;

    int
        m_mode = 1,
        numImages = 9,
        m_midpoint = numImages / 2,
        m_maxiter = 1000,

        m_print = 1;

    double
        m_sprneb = 5.0,
        m_engtol = 0.0001,
        m_frctol = 0.01,
        m_gnorm_tol = 0.01;

    bool
        debug = false;

    void readNudgedElasticBand(std::ifstream& inStream, std::ofstream& outstream);

    NudgeCntrl() {};
    virtual ~NudgeCntrl() {};
};

