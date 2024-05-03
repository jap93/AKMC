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
#include "DimerCntrl.h"
#include "KMCCntrl.h"
#include "NudgeCntrl.h"
#include "MDCntrl.h"

struct JobControl
{
	int
        /** the number of cores per calculation */
        numGrps = 0,

        /** seed for random number generator */
        seed = 0,

        /** the number of MPI cores used for the parent process that coordinates work groups */
        parentCores = 1;

    bool
        /** pre-relax - relax the initial structure */
        preRelax = false,

        /** flag to indicate a restart of the program */
        restart = false;



    std::string 
        inFormat = "dlpoly",
        //geomUnit = "",
        //energyUnit = "",
        externalFile = "onetep.dat",
        pubRootName = "onetep";
        
    std::vector<std::string>
        /** vector of names of atom types to be frozen in simulations */
        frozenTypes;

    KMCCntrl
        kmcParameters;

    RelaxCntrl
        minParameters;

    NudgeCntrl
        nebParameters;

    ARTCntrl
        artParameters;

    DimerCntrl
        dimerParameters;

    MDCntrl
        mdParameters;

	void readJobControl(std::ifstream& instream, std::ofstream& outstream);

    void writeSimulationDetails(const MPIComms& mpi, std::ofstream& outStream);

    

    JobControl() {};
    virtual ~JobControl() {};
};

