#pragma once

#include "Memory.h"

struct NbrListPBC
{

    int
        /** the maximum number of neighbors */
        maxNbrs = 0;

    int
        /* the number of neighbours of a given atom */
        *numNbrs = nullptr;

    int
        /** stores the numbers of atoms that neighbour a given atom */
        *map = nullptr;

    double
        /** the reference positions between nbrlist searches */
        *rZeroX = nullptr,
        *rZeroY = nullptr,
        *rZeroZ = nullptr;

    NbrListPBC();
    virtual ~NbrListPBC();

    /**
     * allocates nbr list arrays
     * @param natoms (int) : the number of atoms
     * @param maxNbrs (int) : the maximum number of neighbors usually natoms/2
     */
    void allocNbrListArrays(int natoms, int maxNbrs);

    /**
     * Calculates the full neighbourlist of interacting atoms
     * @param mpi (const MPIComms&) : the mpi communicator
     * @param posX (double*) : position of atom along X
     * @param posY (double*) : position of atom along Y
     * @param posZ (double*) : position of atom along Z
     * @param latVector (double*) : lattice vector
     * @param rcpVector (double*) : reciprocal space vector
     * @param frozen (int*) : not used here but are frozen atoms for cluster
     * @param natoms   (int) : the number of atoms/particles
     * @param cutoff (double) : the real space cutoff
     */
    void calculateNbrList(double* posX, double* posY, double* posZ,  double* latVector, double* rcpVector, int* frozen, int natoms, double cutoff);

    /**
     * Calculates the full neighbourlist of interacting atoms
     * @param mpi (const MPIComms&) : the mpi communicator
     * @param i (int) : the atom number
     * @param posX (double*) : position of atom along X
     * @param posY (double*) : position of atom along Y
     * @param posZ (double*) : position of atom along Z
     * @param latVector (double*) : lattice vector
     * @param rcpVector (double*) : reciprocal space vector
     * @param frozen (int*) : not used here but are frozen atoms for cluster
     * @param natoms   (int) : the number of atoms/particles
     * @param cutoff (double) : the real space cutoff
     */
    void calculateAtomNbrList(int i, double* posX, double* posY, double* posZ,  double* latVector, double* rcpVector, int* frozen, int natoms, double cutoff);

    /**
     * Determines whether the neighbourlist needs to be recalculated
     * @param velocityX (double*) : the volicity along the x-axis
     * @param velocityY (double*) : the volicity along the y-axis
     * @param velocityZ (double*) : the volicity along the z-axis
     * @param verletShell (double) : the length of the buffer around cutoof
     * @return bool          : returns true if the nbrlist reqyuires recalculation
     */
    bool checkNbrList(int i, double deltaX, double deltaY, double deltaZ, double verletShell);
};

