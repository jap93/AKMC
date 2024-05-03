#pragma once

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <complex>
#include <math.h>

#ifdef __NVCC__

    #include <cuda_runtime_api.h>
    #include <cuda.h>
    #include <cuComplex.h>

    #include "cuComplexExtra.cuh"

#endif

#ifdef SYCL

        #include <CL/sycl.hpp>
        using namespace sycl;

#endif

#include "Memory.h"
#include "Constants.h"
#include "Basis.h"
#include "MPICommunicator.h"

class Ewald
{
private:

    double
        etaMad,
        rcpSpaceCut;                    // real space cutoff

    int
        gcellx = 1,
        gcelly = 1,
        gcellz = 1;                           // number of cells over which the real space summation

    int
        maxGVec = 0;

    int
        numGVec = 0;                         // number of reciprocal space vectors

    double
        *gx = nullptr,
        *gy = nullptr,
        *gz = nullptr,
        *g0 = nullptr,
        *g1 = nullptr;                   
/*
    #ifdef __NVCC__

        double
            *d_g0 = nullptr,
            *d_gx = nullptr,
            *d_gy = nullptr,
            *d_gz = nullptr,
            *d_g1 = nullptr;

    #endif
*/
    bool newJob = true;

    std::complex<double>
        /** atom rciprocal space sums */
        *rcpSum = nullptr,
        /** working store for the above */
        *tmpRcpSum = nullptr;


    #ifdef __NVCC__

        cuDoubleComplex
            *d_rcpSum,
            *d_tmpRcpSum;

    #endif

    #ifdef SYCL

        queue ewaldQueue;

    #endif


public:

    Ewald();
    virtual ~Ewald();

    void printNumGVec(const MPIComms& mpi, std::ofstream&);

    void setEwaldImage(double, double, double, double, double*latVector, double* rcpVector, std::ofstream&);

    void printEwald(const MPIComms& mpi, std::ofstream&);

    void estimateGvectorvoid(const MPIComms& mpi, double volume, double* rcpVector);

    void findGvector(const MPIComms& mpi, double, double* , bool);

    double recipSpaceEnergy(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* charge, int numAtoms);

    double recipSpaceForce(const MPIComms& mpi, double* __restrict__ posX, double*  __restrict__ posY, double*  __restrict__ posZ,
                            double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__, double* __restrict__ atomcharge, int numatoms);

    void atomMoveRecipSpaceEnergy(double* oldPos, double* newPos, double atomCharge, double& engrcpOld, double& engrcpNew);

    void atomSwapRecipSpaceEnergy(double* pos, double oldCharge, double newCharge, double& engrcpOld, double& engrcpNew);

    double removeFragmentRecip(Basis& basin, int* atomList, int numatoms);

    double insertFragmentRecip(Basis& basin, int numAtomsStart, int numAtomsFin);

    void updatePhaseSum(void);

    void updateInsertPhaseSum(void);

    void updateRemovePhaseSum(void);

    void allocateArrays(void);

};

