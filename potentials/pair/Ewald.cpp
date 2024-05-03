#include "Ewald.h"

void dcell(double*, double[10]);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GPU kernals
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// the C++ code starts here
//////////////////////////////////////////////////////////////////////

Ewald::Ewald()
{
    rcpSpaceCut = 0.0;

    gcellx = gcelly = gcellz = 0;

    gx = nullptr;
    gy = nullptr;
    gz = nullptr;

    g0 = nullptr;
    rcpSum = nullptr;
    tmpRcpSum = nullptr;

}

Ewald::~Ewald()
{

    free_dvector(gx);
    free_dvector(gy);
    free_dvector(gz);
    free_dvector(g0);

    free_cvector(rcpSum);
    free_cvector(tmpRcpSum);

}

void Ewald::allocateArrays(void)
{

#ifdef SYCL

    gx = malloc_shared<double>(maxGVec, ewaldQueue);
    gy = malloc_shared<double>(maxGVec, ewaldQueue);
    gz = malloc_shared<double>(maxGVec, ewaldQueue);

    g0 = malloc_shared<double>(maxGVec, ewaldQueue);
    rcpSum = malloc_shared<std::complex<double>>(maxGVec, ewaldQueue);
    tmpRcpSum = malloc_shared<std::complex<double>>(maxGVec, ewaldQueue);

#else
    gx = alloc_dvector(maxGVec, "ewald", 0.0);
    gy = alloc_dvector(maxGVec, "ewald", 0.0);
    gz = alloc_dvector(maxGVec, "ewald", 0.0);

    g0 = alloc_dvector(maxGVec, "ewald", 0.0);
    g1 = alloc_dvector(maxGVec, "ewald", 0.0);
    rcpSum = alloc_cvector(maxGVec, "ewald");
    tmpRcpSum = alloc_cvector(maxGVec, "ewald");

#endif

    
    
}

/* error function */

void Ewald::printNumGVec(const MPIComms& mpi, std::ofstream& outStream)
{
    
    outStream << "\n the number of reciprocal lattice vectors " << numGVec << std::endl;
}



/* *************************************************************
 * calculates the ewald parameters usiming the image convention
 *
 * input: accuracy required for madelung sum
 ***************************************************************/

void Ewald::setEwaldImage(double eta, double cutoff, double eps, double tol, double*latVector, double* rcpVector, std::ofstream& outStream)
{
    double
        cellx,                        // cell dimensions in x y + z
        celly,
        cellz,
        tol1,
        celprp[10];

    etaMad = eta;

    /* determines the maximum number of k vectors
     * in the ewald summation */
    dcell(latVector, celprp);

    tol1 = sqrt(-log(eps * cutoff * pow(2.0 * tol * eta, 2)));

    gcellx = rint(0.25 + celprp[0] * eta * tol1 / PI);
    gcelly = rint(0.25 + celprp[1] * eta * tol1 / PI);
    gcellz = rint(0.25 + celprp[2] * eta * tol1 / PI);


    dcell(rcpVector, celprp);
    cellx = double(gcellx) * celprp[6];
    celly = double(gcelly) * celprp[7];
    cellz = double(gcellz) * celprp[8];

    rcpSpaceCut = std::min(cellx, celly);
    rcpSpaceCut = std::min(rcpSpaceCut, cellz);

    return;

}


void Ewald::printEwald(const MPIComms& mpi, std::ofstream& outStream)
{
    
    outStream << "\n cutoffs for lattice summation:"
              << "\n eta                     = " << etaMad
              << "\n reciprocal space cutoff = " << rcpSpaceCut << " angstroms^-1" << std::endl;

    outStream << "\n maximum k-vectors "
              << std::setw(4) << gcellx
              << std::setw(4) << gcelly
              << std::setw(4) << gcellz << std::endl;

    outStream.flush();

}

void Ewald::estimateGvectorvoid(const MPIComms& mpi, double volume, double* rcpVector)
{
    int
        nx,
        ny,
        nz;

    double
        x,
        y,
        z,
        xx,
        yy,
        zz,
        gsq = 0.0,
        rcpcutsq = rcpSpaceCut * rcpSpaceCut;

    maxGVec = 0; // zero number of vectors

    for (nx = 0; nx <= gcellx; nx++)
    {

        for (ny = -gcelly; ny <= gcelly; ny++)
        {

            for (nz = -gcellz; nz <= gcellz; nz++)
            {

                x = nx * rcpVector[0] + ny * rcpVector[1] + nz * rcpVector[2];
                y = nx * rcpVector[3] + ny * rcpVector[4] + nz * rcpVector[5];
                z = nx * rcpVector[6] + ny * rcpVector[7] + nz * rcpVector[8];

                gsq = x * x + y * y + z * z;

                if (gsq < rcpcutsq && gsq > 1.0e-8)
                {
                    maxGVec++;
                }
            }
        }
    }

    maxGVec += 100; // add a bit of padding
}

/* calculates g-vectors within a given radius.
 * returns the number of g vectors.
 * nb update every time volume called */

void Ewald::findGvector(const MPIComms& mpi, double volume, double* rcpVector, bool flag)
{
    int
        nx,
        ny,
        nz;

    double
        gw,
        gzz,
        etasq = etaMad * etaMad,
        gsqfct = 1.0 / (4.0 * etasq),
        gfct0 = CTOEV * PI / (etasq * volume),
        gfct1 = gsqfct * gfct0;;

    double
        x,
        y,
        z,
        xx,
        yy,
        zz,
        gsq = 0.0,
        rcpcutsq = rcpSpaceCut * rcpSpaceCut;

    numGVec = 0; // zero number of vectors

    for (nx = 0; nx <= gcellx; nx++)
    {

        for (ny = -gcelly; ny <= gcelly; ny++)
        {

            for (nz = -gcellz; nz <= gcellz; nz++)
            {

                x = nx * rcpVector[0] + ny * rcpVector[1] + nz * rcpVector[2];
                y = nx * rcpVector[3] + ny * rcpVector[4] + nz * rcpVector[5];
                z = nx * rcpVector[6] + ny * rcpVector[7] + nz * rcpVector[8];

                gsq = x * x + y * y + z * z;

                if (gsq < rcpcutsq && gsq > 1.0e-8)
                {

                    //add the new vector and increase counter
                    gx[numGVec] = x;
                    gy[numGVec] = y;
                    gz[numGVec] = z;

                    numGVec++;

                    if (numGVec == maxGVec)
                    {

                        std::cout << "\n\n***too many g vectors";
                        std::cout.flush();
                        exit(EXIT_FAILURE);

                    }
                }
            }
        }
    }
    
#ifdef SYCL

    ewaldQueue.submit([&](handler &h) {

        double* lgx = gx;
        double* lgy = gy;
        double* lgz = gz;
        double* lg0 = g0;

        h.parallel_for(numGVec, [=](id<1> i) {

        lgx[i] *= TWOPI;
        lgy[i] *= TWOPI;
        lgz[i] *= TWOPI;

        double gsq = lgx[i] * lgx[i] + lgy[i] * lgy[i] + lgz[i] * lgz[i];
        double x = gsq * gsqfct;

        double gzz = 2.0 * exp(-x) / x;
        gzz = exp(-x) / x;

        if (xx < 1.0e-8)
        {
            lg0[i] = gzz * gfct0;
        }
        else
        {
            lg0[i] = gzz * gfct0 * 2.0;
        }

        });
    });
    ewaldQueue.wait();
#else 
    
    for (u_int32_t i = 0; i < numGVec; i++)
    {

        gx[i] *= TWOPI;
        gy[i] *= TWOPI;
        gz[i] *= TWOPI;

        xx = gx[i];
        yy = gy[i];
        zz = gz[i];
        gsq = xx * xx + yy * yy + zz * zz;
        x = gsq * gsqfct;

        gzz = 2.0 * exp(-x) / x;
        gzz = exp(-x) / x;

        if (xx < 1.0e-8)
        {
            g0[i] = gzz * gfct0;
        }
        else
        {
            g0[i] = gzz * gfct0 * 2.0;
        }
        gw = -2.0 * gzz * (1.0 + 1.0 / x);
        g1[i] = gw * gfct1;

    }

#endif

}

void Ewald::updatePhaseSum(void)
{
    u_int32_t
       ircp;

      
    #ifdef SYCL
        ewaldQueue.submit([&](handler &h) {

        std::complex<double>* locSum = rcpSum;
        std::complex<double>* locTmp = tmpRcpSum;
        h.parallel_for(numGVec, [=](id<1> ircp) {
            // access sharedArray
            locSum[ircp] = locTmp[ircp];
            });
        });
        ewaldQueue.wait();

    #else
        
        for(ircp = 0; ircp < numGVec; ircp++)
            rcpSum[ircp] = tmpRcpSum[ircp];

    #endif

}

void Ewald::updateInsertPhaseSum(void)
{
    int
       ircp;

    for(ircp = 0; ircp < numGVec; ircp++)
        rcpSum[ircp] += tmpRcpSum[ircp];

}

void Ewald::updateRemovePhaseSum(void)
{
    int
       ircp;

    for(ircp = 0; ircp < numGVec; ircp++)
        rcpSum[ircp] -= tmpRcpSum[ircp];

}

double Ewald::recipSpaceEnergy(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* charge, int numAtoms) //double* posX, double* posY, double* posZ, double* atomCharge, int numatoms)
{
    u_int32_t
        rcp;

    double
        pphase,   
        ercp = 0.0;

    std::complex<double> 
        tmpSum1,
        tmpSum2,
        cplxConj,
        COMPLEXI = {0.0, 1.0};

    #ifdef SYCL
        //ewaldQueue.submit([&](handler &h) {
        //h.parallel_for(numGVec, [=](id<1> rcp) {
        //    tmpRcpSum[rcp] = 0.0 + 0.0 * COMPLEXI;
        //    });
        //});
        //ewaldQueue.wait();

    #else
        for (rcp = 0; rcp < numGVec; rcp++)
            tmpRcpSum[rcp] = 0.0 + 0.0 * COMPLEXI;
    #endif

    #ifdef SYCL
        ewaldQueue.submit([&](handler &h) {
        
        double* lgx = gx;
        double* lgy = gy;
        double* lgz = gz;
        double* px = basin.posX;
        double* py = basin.posY;
        double* pz = basin.posZ;
        double* chg = basin.charge;
        std::complex<double>* ltmp = tmpRcpSum;

        h.parallel_for(numGVec, [=](id<1> rcp) {

            ltmp[rcp] = 0.0 + 0.0 * COMPLEXI;

            for (u_int32_t i = 0; i < natoms; i++)
            {
                double pphase = lgx[rcp] * px[i] + lgy[rcp] * py[i] + lgz[rcp] * pz[i];
                std::complex<double> tmpSum1 = COMPLEXI * pphase;
                std::complex<double> tmpSum2 = exp(tmpSum1);
                ltmp[rcp] = ltmp[rcp] + tmpSum2 * chg[i]; //* exp(COMPLEXI * pphase);
            }
        });
        });
        ewaldQueue.wait(); 
         
    #else 
        for (rcp = 0; rcp < numGVec; rcp++)
        {

            for (u_int32_t i = 0; i < numAtoms; i++)
            {
                pphase = gx[rcp] * posX[i] + gy[rcp] * posY[i] + gz[rcp] * posZ[i];
                tmpSum1 = COMPLEXI * pphase;
                tmpSum2 = exp(tmpSum1);
                tmpRcpSum[rcp] = tmpRcpSum[rcp] + tmpSum2 * charge[i]; //* exp(COMPLEXI * pphase);
            }
        
        }
    #endif
    

    //reduction on energy 
    for (rcp = 0; rcp < numGVec; rcp++)
    {
        
        cplxConj = tmpRcpSum[rcp] * conj(tmpRcpSum[rcp]);
        ercp += g0[rcp] * cplxConj.real();

    }

    ercp *= 0.5;

    return ercp;
}
#ifdef _OPENMP
double Ewald::recipSpaceForce(const MPIComms& mpi, double* __restrict__ posX, double*  __restrict__ posY, double*  __restrict__ posZ,
                              double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__ atomcharge, int numatoms)
{
    int
        rcp,
        atom;

    double
        xx, yy, zz,
        gx, gy, gz, g0,
        //dummy,
        cossum,
        sinsum,
        phase,
        phasefact,
        energyrcp = 0.0;


    double* sinphase; //[numatoms];  // cannot dynamically allocate for OMP
    double* cosphase; //[numatoms];  // you may be able to do this within the parallel region
    double* tmpForceX; //[numatoms];
    double* tmpForceY; //[numatoms];
    double* tmpForceZ; //[numatoms];
   
    #pragma omp parallel default(none)  \
        reduction(+: energyrcp) \
        shared(gx, gy, gz, g0, posX, posY, posZ, forceX, forceY, forceZ, atomcharge,  numGVec, numatoms, std::cout ) \
        private(rcp, cossum, sinsum, phase, cosphase, sinphase, phasefact, atom, gx, gy, gz, g0, xx, yy, zz, tmpForceX, tmpForceY, tmpForceZ)
    { 
        sinphase = (double*) calloc(numatoms, sizeof(double));
        cosphase = (double*) calloc(numatoms, sizeof(double));
        tmpForceX = (double*) calloc(numatoms, sizeof(double));
        tmpForceY = (double*) calloc(numatoms, sizeof(double));
        tmpForceZ = (double*) calloc(numatoms, sizeof(double));

        #pragma omp for schedule(static,4)
        for (rcp = 0; rcp < numGVec; rcp++)
        {
            cossum = 0.0;
            sinsum = 0.0;

            gx = gx[rcp];
            gy = gy[rcp];
            gz = gz[rcp];
            g0 = g0[rcp];

            for (atom = 0; atom < numatoms; atom++)
            {

                xx = posX[atom];
                yy = posY[atom];
                zz = posZ[atom];

                phase = gx * xx + gy * yy + gz * zz;

                cosphase[atom] = atomcharge[atom] * cos(phase);
                sinphase[atom] = atomcharge[atom] * sin(phase);

                cossum = cossum + cosphase[atom];
                sinsum = sinsum + sinphase[atom];

            }

            for (atom = 0; atom < numatoms; atom++)
            {
                energyrcp += g0 * (cosphase[atom] * cossum + sinphase[atom] * sinsum);

                phasefact = -g0 * (-sinphase[atom] * cossum + cosphase[atom] * sinsum);

                tmpForceX[atom] += phasefact * gx;
                tmpForceY[atom] += phasefact * gy;
                tmpForceZ[atom] += phasefact * gz;

            }

        } // end of loop over reciprocal lattice vectors
        for (atom = 0; atom < numatoms; atom++)
        {
            #pragma omp atomic
            forceX[atom] += tmpForceX[atom];
        }
        for (atom = 0; atom < numatoms; atom++)
        {
            #pragma omp atomic
            forceY[atom] += tmpForceY[atom];
        }
        for (atom = 0; atom < numatoms; atom++)
        {
            #pragma omp atomic
            forceZ[atom] += tmpForceZ[atom];
        }

        free(tmpForceX);
        free(tmpForceY);
        free(tmpForceZ);
        free(sinphase);
        free(cosphase);

    }
    // energy, forces and stresses must be converted to eV
    energyrcp *= 0.5;

    return energyrcp;
}

#else

double Ewald::recipSpaceForce(const MPIComms& mpi, double* __restrict__ posX, double*  __restrict__ posY, double*  __restrict__ posZ,
                              double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__ stress, 
                              double* __restrict__ atomcharge, int numatoms)
{
    int
        rank = mpi.getWorkGroupRank(),
        numProcs = mpi.getNumWorkGroupProcs();

    double
        xx, yy, zz,
        vx, vy, vz, v0, v1,
        //dummy,
        cossum,
        sinsum,
        phase,
        phasefact,
        phSumSq,
        phaseTotal,
        factor,
        tmp_energy = 0.0,
        energyrcp = 0.0;

    double* sinphase = alloc_dvector(numatoms, "ewald force", 0.0);
    double* cosphase = alloc_dvector(numatoms, "ewald force", 0.0);

    for (int rcp = rank; rcp < numGVec; rcp = rcp + numProcs)
    {
        cossum = 0.0;
        sinsum = 0.0;

        vx = gx[rcp];
        vy = gy[rcp];
        vz = gz[rcp];
        v0 = g0[rcp];
        v1 = g1[rcp];

        for (int atom = 0; atom < numatoms; atom++)
        {

            xx = posX[atom];
            yy = posY[atom];
            zz = posZ[atom];

            phase = vx * xx + vy * yy + vz * zz;

            cosphase[atom] = atomcharge[atom] * cos(phase);
            sinphase[atom] = atomcharge[atom] * sin(phase);

            cossum = cossum + cosphase[atom];
            sinsum = sinsum + sinphase[atom];

        }

        phSumSq = 0.5 * (cossum * cossum + sinsum * sinsum);

        if(phSumSq > 1.0e-16)
        {
            phaseTotal += phSumSq * v0;

            factor = 2.0 * phSumSq * v1;

            // add in stress 
            stress[0] += factor * vx * vx;
            stress[1] += factor * vy * vy;
            stress[2] += factor * vz * vz;
            stress[3] += factor * vy * vz;
            stress[4] += factor * vx * vz;
            stress[5] += factor * vx * vy;

        }


        //#pragma clang loop vectorize(enable)
        //#pragma omp simd aligned(forceX:64)
        //#pragma vector aligned //all vectors aligned

        for (int atom = 0; atom < numatoms; atom++)
        {
            tmp_energy += v0 * (cosphase[atom] * cossum + sinphase[atom] * sinsum);

            phasefact = -v0 * (-sinphase[atom] * cossum + cosphase[atom] * sinsum);

            forceX[atom] += phasefact * vx;
            forceY[atom] += phasefact * vy;
            forceZ[atom] += phasefact * vz;

        }

    } // end of loop over reciprocal lattice vectors

    phaseTotal = mpi.sumDoubleGroup(phaseTotal);
    stress[0] -= phaseTotal;
    stress[1] -= phaseTotal;
    stress[2] -= phaseTotal;
    
    tmp_energy *= 0.5;
    energyrcp = mpi.sumDoubleGroup(tmp_energy);

    free_dvector(sinphase);
    free_dvector(cosphase);

    return energyrcp;
}

#endif


 
 //calculates the reciprocal lattice part of ewald energy
 void Ewald::atomMoveRecipSpaceEnergy(double* oldPos, double* newPos, double atomCharge, double& engrcpOld, double& engrcpNew)
 {
    int
        pid = 0, //mpi.getWorkGroupRank(),
        numGrp = 1, //mpi.getNumWorkGroupProcs(),
        ircp;
 
    double
        
        pphaseOld,
        pphaseNew;
         
    std::complex<double> 
        cplxConj,
        COMPLEXI = {0.0,1.0};
 
    for (ircp = pid; ircp < numGVec; ircp += numGrp)    
        tmpRcpSum[ircp] = rcpSum[ircp];
 
    for (ircp = pid; ircp < numGVec; ircp += numGrp) 
    {
        cplxConj = tmpRcpSum[ircp] * conj(tmpRcpSum[ircp]);
        engrcpOld += g0[ircp] * cplxConj.real();
    }
 
    for (ircp = pid; ircp < numGVec; ircp += numGrp)
    {
 
        //old phase factor
        pphaseOld = gx[ircp] * oldPos[0] + gy[ircp] * oldPos[1] + gz[ircp] * oldPos[2];
 
        //new phase 
        pphaseNew = gx[ircp] * newPos[0] + gy[ircp] * newPos[1] + gz[ircp] * newPos[2];;
 
        //update the sum by calculating differences
        tmpRcpSum[ircp] = tmpRcpSum[ircp] + atomCharge * (exp(pphaseNew * COMPLEXI) - exp(pphaseOld * COMPLEXI));
        cplxConj = tmpRcpSum[ircp] * conj(tmpRcpSum[ircp]);
        engrcpNew += g0[ircp] * cplxConj.real();
 
    }
 
    engrcpOld *= 0.5;
    engrcpNew *= 0.5;
}
 
//calculates the reciprocal lattice part of ewald energy
void Ewald::atomSwapRecipSpaceEnergy(double* pos, double oldCharge, double newCharge, double& engrcpOld, double& engrcpNew)
 
{
    int
        ircp;
 
    double
        pphase,
        deltachg = newCharge - oldCharge;
         
 
    std::complex<double> 
        cplxConj,
        COMPLEXI = {0.0,1.0};
 
    for (ircp = 0; ircp < numGVec; ircp++)  
        tmpRcpSum[ircp] = rcpSum[ircp];
 
    for (ircp = 0; ircp < numGVec; ircp++)
    {
        cplxConj = tmpRcpSum[ircp] * conj(tmpRcpSum[ircp]);
        engrcpOld += g0[ircp] * cplxConj.real();
    }
 
    for (ircp = 0; ircp < numGVec; ircp++)
    {
 
        //new phase 
        pphase = gx[ircp] * pos[0] + gy[ircp] * pos[1] + gz[ircp] * pos[2];;
 
        //update the sum by calculating differences
        tmpRcpSum[ircp] = tmpRcpSum[ircp] + deltachg * exp(pphase * COMPLEXI);
        cplxConj = tmpRcpSum[ircp] * conj(tmpRcpSum[ircp]);
        engrcpNew += g0[ircp] * cplxConj.real();
 
    }
 
}

double Ewald::removeFragmentRecip(Basis& basin, int* atomList, int numatoms)
{
    int
        rcp,  
        j,
        i;

    double
        pphase,   
        ercp = 0.0;

    std::complex<double> 
        cplxConj,
        COMPLEXI = {0.0,1.0};

    for (rcp = 0; rcp < numGVec; rcp++)
        //tmpRcpSum[rcp] = 0.0 + 0.0 * COMPLEXI;
        tmpRcpSum[rcp] = rcpSum[rcp];

    for (rcp = 0; rcp < numGVec; rcp++)
    {

        for (i = 0; i < numatoms; i++)
        {
            j = atomList[i];
            pphase = gx[rcp] * basin.posX[j] + gy[rcp] * basin.posY[j] + gz[rcp] * basin.posZ[j];
            tmpRcpSum[rcp] = tmpRcpSum[rcp] - basin.charge[j] * exp(pphase * COMPLEXI);
        }

    }

    for (rcp = 0; rcp < numGVec; rcp++)
    {
        cplxConj = tmpRcpSum[rcp] * conj(tmpRcpSum[rcp]);
        ercp += g0[rcp] * cplxConj.real();
    }

    //if (numProcs > 1)  // sum the energies

    return ercp;
}

double Ewald::insertFragmentRecip(Basis& basin, int numAtomsStart, int numAtomsFin)
{
    int
        rcp,  
        i;

    double
        pphase,   
        ercp = 0.0;

    std::complex<double> 
        cplxConj,
        COMPLEXI = {0.0,1.0};

    for (rcp = 0; rcp < numGVec; rcp++)
    {
        tmpRcpSum[rcp] = 0.0 + 0.0 * COMPLEXI;
        //tmpRcpSum[rcp] = rcpSum[rcp];
    }

    for (rcp = 0; rcp < numGVec; rcp++)
    {

        for (i = numAtomsStart; i < numAtomsFin; i++)
        {
            pphase = gx[rcp] * basin.posX[i] + gy[rcp] * basin.posY[i] + gz[rcp] * basin.posZ[i];
            tmpRcpSum[rcp] = tmpRcpSum[rcp] + basin.charge[i] * exp(pphase * COMPLEXI);
        }

    }

    for (rcp = 0; rcp < numGVec; rcp++)
    {
        cplxConj = tmpRcpSum[rcp] * conj(tmpRcpSum[rcp]);
        ercp += g0[rcp] * cplxConj.real();
    }

    //if (numProcs > 1)  // sum the energies

    return ercp;

    /*   cuDoubleComplex
        tmp1, tmp2,
        COMPLEXI = make_cuDoubleComplex(0.0, 1.0);

    while (tid < numGVec)
    {
        d_tmpRcpSum[tid] = d_rcpSum[tid]; //make_cuDoubleComplex(0.0, 0.0); //

        for (unsigned int i = numAtomStart; i < numAtomFin; i++)
        {

            phase = gx[tid] * posX[i] + gy[tid] * posY[i] + gz[tid] * posZ[i];

            tmp1 = cuCmul(COMPLEXI, phase);
            tmp2 = cuCexp(tmp1);
            d_tmpRcpSum[tid] = cuCadd(d_tmpRcpSum[tid], cuCmul(tmp2, atomCharge[i])); 
            //printf("\n phase %10.4f, %10.4f, %10.4f, %10.4f",gx[tid],gy[tid],gz[tid],phase);
        }

        tid += blockDim.x * gridDim.x;
    }*/
}
 
