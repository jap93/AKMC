#include "Field.h"
void dcell(double*, double[10]);
std::vector<std::string> split(std::string s);

//int omp_get_thread_num(void);

Field::Field()
{
    
}

Field::~Field()
{
    free_dvector(tmpRho);
    free_dvector(wrkRho);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// utility functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Field::setup(const MPIComms& mpi, double* latVector, double* rcpVector, double minDimension, double volume, 
                  int maxAtoms, const Species& spec, std::ofstream& outStream)
{


    double
        radius,
        eps,
        tol,
        celprp[10];

    if (shortRangeCut < 1.0)
    {
        outStream << "\n\n\n *** short range cutoff is too small : " << shortRangeCut << std::endl;
        outStream.flush();
        mpi.commsAbortWorld();
        exit(EXIT_FAILURE);
    }

    if (useNbrList)
    {
        radius = shortRangeCut + verletShell;
    }
    else
    {
        radius = shortRangeCut;
    }

    // check the image convention is obeyed
    if (minimumImage && radius > (minDimension / 2))
    {
        outStream << "\n\n\n *** short range shortRangeCut to large for simulation box." << std::endl;
        outStream << " cutoff : " << shortRangeCut << " minimum dimension: " << minDimension << std::endl;
        outStream.flush();
        mpi.commsAbortWorld();
        exit(EXIT_FAILURE);
    }

    if (useNbrList && verletShell < 1e6)
    {
        outStream << "\n\n\n *** WARNING: the verlet does not appear to have been set." << std::endl;
        //outStream.flush();
        //mpi.commsAbortWorld();
        //exit(EXIT_FAILURE);
    }

    if (!minimumImage)   // setup for non-Image convention
    {
        dcell(latVector, celprp);
        cellX = int(shortRangeCut / celprp[0]) + 1;
        cellY = int(shortRangeCut / celprp[1]) + 1;
        cellZ = int(shortRangeCut / celprp[2]) + 1;  

        outStream << "\n\n Image convention is not being used!" << std::endl;
        outStream << "\n The number of cells used in the summation are " << cellX << " " << cellY << " " << cellZ << std::endl;       
    }
    else
    {
        outStream << "\n\n\n short range cutoffs for lattice summation (image):"
                  << "\n real space shortRangeCut       = " << shortRangeCut << " Angstroms" << std::endl;
    }
    
    //alocate the density work arrays
    outStream << "\n allocating temporary metal density arrays of size " << maxAtoms << std::endl;

    tmpRho = alloc_dvector(maxAtoms, "rho", 0.0);
    wrkRho = alloc_dvector(maxAtoms, "rho", 0.0);
    wrkSize = maxAtoms;

    // set up short-range potential mesh
    mbdy.calculateEnergyMesh(spec, shortRangeCut, outStream);


}

void Field::print(const MPIComms& mpi, std::ofstream& outStream)
{

    
    /**** key words for internal energy evaluation i.e. interatomic potentials */

	outStream << " short range shortRangeCut                                           " << shortRangeCut << std::endl;

    if (useNbrList)
	{
		outStream << " a neighbourlist will be used with a buffer of " << verletShell << " A" << std::endl;
		outStream << " the maximum number of neighbours is " << maximumNbrs << std::endl;
	}

    mbdy.printPotential(outStream);
}




/******************************************************
reads in description of the species and the potential model
*******************************************************/
void Field::readPotential(std::ifstream& inStream, std::ofstream& outStream, Species& spec)
{
    std::string
        line,
        subWord,
        keyWord,
        dummy,
        dummy2;

    int
        num;

    bool
       finishRead = false;

    // read in and print out the species and potentails
    std::vector<std::string> words;

    while (!inStream.eof() || finishRead == false)
    {                                             // or start directive

        std::getline(inStream, line);
        transform(line.begin(), line.end(), line.begin(), ::tolower);
	    words = split(line);

        keyWord = words[0];

        if (keyWord == "cutoff")
		{
			dummy = words[1];
			shortRangeCut = std::stod(dummy);
		}
		else if (keyWord == "nbrlist")
		{
			useNbrList = true;
			dummy = words[1];
			verletShell = std::stod(dummy);
		}
        else if (keyWord == "maxNbrs")
		{
			dummy = words[1];
			maximumNbrs = std::stoi(dummy);
		}
        
		else if (keyWord == "noimage")
		{
			minimumImage = false;
		}
        else if (keyWord == "species")
        {

            dummy = words[1];
            num = std::stoi(dummy);
            spec.loadSpecies(inStream, num);
            spec.printSpecies(outStream);
        }
        else if (keyWord == "many" || keyWord == "manybody")
        {

            dummy = words[1];
            num = std::stoi(dummy);

            // get the correct energyUnits
            dummy = words[2];

            mbdy.loadPotential(num, inStream, outStream);
            //mbdy.printPotential(outStream);

        }
        else if (keyWord == "close")
        {
            finishRead = true;
            break;
        }
        else
        {
            outStream << "\n unrecognised key word used in potentials file " << keyWord << std::endl;
            outStream.flush();
            throw std::invalid_argument("Invalid keyword!");
        }

    }
    outStream.flush();
}


/************************************************************************
 driver for the calculation of forces and energy of all ions
 ************************************************************************/
void Field::calculateForces(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ, 
                            double* stress, double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int* frozen, int numAtoms, Energy& eng)
{
    double
        manyEnergy = 0.0,
        rcpenergy = 0.0,
        twoenergy = 0.0;

    for (int i = 0; i < numAtoms; i++)
    {
        forceX[i] = 0.0;
        forceY[i] = 0.0;
        forceZ[i] = 0.0;
    }
    for (int i = 0; i < 6; i++)
        stress[i] = 0.0;


    
    if (minimumImage)
    {
        resetSimulationCellMinImage(mpi, posX, posY, posZ, latVector, rcpVector, numAtoms);

        if (useNbrList)
        {
            realSpaceForceImageNbrList(mpi, manyEnergy, twoenergy, posX, posY, posZ, forceX, forceY, forceZ, stress, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
        }
        else
        {
                realSpaceForceImage(mpi, manyEnergy, twoenergy, posX, posY, posZ, forceX, forceY, forceZ, stress, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
        }
               
    }
    else
    {
        resetSimulationCellNoImage(mpi, posX, posY, posZ, latVector, rcpVector, numAtoms);
        realSpaceForceNoImage(mpi, manyEnergy, twoenergy, posX, posY, posZ, forceX, forceY, forceZ, stress, latVector, rcpVector, 
                                atomCharge, atmLabel, numAtoms);

        
    }

    eng.manyEnergy = manyEnergy;
    eng.vdwEnergy = twoenergy;
    
    eng.totalEnergy = rcpenergy + manyEnergy + twoenergy;

    if (mpi.getNumWorkGroupProcs() > 1)
    {
        sumForces(mpi, forceX, forceY, forceZ, numAtoms);
        mpi.sumDoubleVectorGroup(stress, 6);
    }

    /*if (mpi.getWorkGroupRank() == 0){
    for (int i = 0; i < numAtoms; i++)
    {
        std::cout << forceX[i] << " " << forceY[i] << " " << forceZ[i] << std::endl;
    }
    }

    std::cout.flush();
    exit(0); */
    
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// routines using the nearest image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/**************************************************************
 * calculate energy and forces using image convention- this subroutine
 * uses explicit calculation of the energy and forces
 *
 * before doing the calculation it is necessary to
 * to reset the box so that any atom that has
 * migrated out of the box can be placed back in according
 * to the image convention
 **************************************************************/
#if defined  _OPENMP


void Field::realSpaceForceImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* __restrict__ posX, double* __restrict__ posY, double* __restrict__ posZ, 
                                double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__ stress, double* __restrict__ latVector,
                                double* __restrict__ rcpVector, double* __restrict__ atomCharge, int* __restrict__ atmLabel, int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_s,
        tmp_energy_s = 0.0,
        tmp_energy_m = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just squared

    double
        *tmpForceX = nullptr, 
        *tmpForceY = nullptr, 
        *tmpForceZ = nullptr,
        *tmpStress = nullptr,
        *tmpRho = nullptr;

    manyEnergy = 0.0;
    twoenergy = 0.0;

    #pragma omp parallel default(none)  \
        reduction(+: manyEnergy, tmp_energy_s) \
        shared(atmLabel, posX, posY, posZ, numAtoms, rcpVector, latVector, forceX, forceY, forceZ, wrkRho, stress, radius, shortRangeCut, std::cout) \
        private(ltypea, ltypeb, ax, ay, az, bx, by, bz, rx, ry, rz, r, rsq, xx, yy, zz, aatom, batom) \
        private(eng_s, forc_m, forc_s, forceTot, tmp_energy_m, tmpRho, tmpForceX, tmpForceY, tmpForceZ, tmpStress)
    {
        tmpForceX = (double*) calloc(numAtoms, sizeof(double));
        tmpForceY = (double*) calloc(numAtoms, sizeof(double));
        tmpForceZ = (double*) calloc(numAtoms, sizeof(double));
        tmpRho = (double*) calloc(numAtoms, sizeof(double));
        tmpStress = (double*) calloc(6, sizeof(double));

        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            ltypea = atmLabel[aatom];
            tmpRho[aatom] = calculateAtomRho(aatom, posX, posY, posZ, latVector, rcpVector, atmLabel, numAtoms);
        }

        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            //#pragma omp atomic
            wrkRho[aatom] = tmpRho[aatom];
        }

        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            tmp_energy_m = mbdy.embed(ltypea, ltypea, tmpRho[aatom]);
            #pragma omp atomic
            manyEnergy += tmp_energy_m;
        }


        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {

            ltypea = atmLabel[aatom];

            ax = posX[aatom];
            ay = posY[aatom];
            az = posZ[aatom];

            for (batom = 0; batom < numAtoms; batom++)
            {
                if (aatom == batom)
                    continue;

                ltypeb = atmLabel[batom];

                bx = posX[batom];
                by = posY[batom];
                bz = posZ[batom];

                rx = ax - bx;
                ry = ay - by;
                rz = az - bz;

                xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
                yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
                zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

                xx -= rint(xx);
                yy -= rint(yy);
                zz -= rint(zz);

                rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
                ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
                rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

                rsq = rx * rx + ry * ry + rz * rz;

                eng_s = 0.0;
                forc_m = 0.0;
                forc_s = 0.0;

                if (rsq <= radius)
                {

                    r = sqrt(rsq);

                    mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                    //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;
                
                    mbdy.calculateManyBodyForce(ltypea, ltypeb, r, rsq, wrkRho[aatom], wrkRho[batom], forc_m);
                    //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                    //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                    //exit(-1);
                
                    forceTot = (forc_s - forc_m) / (rsq * 2.0);
            
                    tmp_energy_s += 0.5 * eng_s;

                    tmpForceX[aatom] += rx * forceTot;
                    tmpForceY[aatom] += ry * forceTot;
                    tmpForceZ[aatom] += rz * forceTot;
                    tmpForceX[batom] -= rx * forceTot;
                    tmpForceY[batom] -= ry * forceTot;
                    tmpForceZ[batom] -= rz * forceTot;

                    tmpStress[0] -= forceTot * rx * rx;
                    tmpStress[1] -= forceTot * ry * ry;
                    tmpStress[2] -= forceTot * rz * rz;
                    tmpStress[3] -= forceTot * ry * rz;
                    tmpStress[4] -= forceTot * rx * rz;
                    tmpStress[5] -= forceTot * rx * ry;

                }

            } // end of loop over batom
 
        } // end of loop over aatom

        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            #pragma omp atomic
            forceX[aatom] += tmpForceX[aatom];
        }
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            #pragma omp atomic
            forceY[aatom] += tmpForceY[aatom];
        }
        for (aatom = 0; aatom < numAtoms; aatom++)
        {
            #pragma omp atomic
            forceZ[aatom] += tmpForceZ[aatom];
        }
        for (int i = 0; i < 6; i++)
        {
            #pragma omp atomic
            stress[i] += tmpStress[i];
        }

        //finally delete memory
        free(tmpForceX);
        free(tmpForceY);
        free(tmpForceZ);
        free(tmpRho);
        free(tmpStress);
    } //end of parallel

    
    twoenergy = tmp_energy_s;

    /*
    for (aatom = 0; aatom < numAtoms; aatom++)
        std::cout << "frc " << aatom << " " << forceX[aatom] << " " << forceY[aatom] << " " << forceZ[aatom] << std::endl;

    std::cout << "\n two, many mody energy " << twoenergy << " " << manyEnergy << std::endl;
    exit(0);
    */
}

#else

void Field::realSpaceForceImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* __restrict__ posX, double* __restrict__ posY, double* __restrict__ posZ, 
                                double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__ stress, double* __restrict__ latVector,
                                double* __restrict__ rcpVector, double* __restrict__ atomCharge, int* __restrict__ atmLabel, int numAtoms)
{
    int
        rank = mpi.getWorkGroupRank(),
        numProcs = mpi.getNumWorkGroupProcs(),
        batom,
        ltypea,
        ltypeb,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        fxi, fyi, fzi,
        eng_s,
        tmp_energy_s = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just squared

    manyEnergy = 0.0;
    twoenergy = 0.0;

    for (aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
    //for (aatom = 0; aatom < numAtoms; aatom++)
    {
        ltypea = atmLabel[aatom];
        //wrkRho[aatom] = 0.0;
        wrkRho[aatom] = calculateAtomRho(aatom, posX, posY, posZ, latVector, rcpVector, atmLabel, numAtoms);
        //std::cout << "rho " << rank << " " << numProcs << " " << aatom << " " << wrkRho[aatom] << std::endl;
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }

    if (numProcs > 1)
    {
        mpi.sumDoubleVectorGroup(wrkRho, numAtoms);
        manyEnergy = mpi.sumDoubleGroup(manyEnergy);
    }
    

    for (aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
    //for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        fxi = 0.0;
        fyi = 0.0;
        fzi = 0.0;

        for (batom = 0; batom < numAtoms; batom++)
        {
            if (aatom == batom)
                continue;

            ltypeb = atmLabel[batom];

            bx = posX[batom];
            by = posY[batom];
            bz = posZ[batom];

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
            yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
            zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            eng_s = 0.0;
            forc_m = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;
                
                mbdy.calculateManyBodyForce(ltypea, ltypeb, r, rsq, wrkRho[aatom], wrkRho[batom], forc_m);
                //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                //exit(-1);
                
                forceTot = (forc_s - forc_m) / (rsq * 2.0);
            
                tmp_energy_s += 0.5 * eng_s;

                fxi += rx * forceTot;
                fyi += ry * forceTot;
                fzi += rz * forceTot;
                forceX[batom] -= rx * forceTot;
                forceY[batom] -= ry * forceTot;
                forceZ[batom] -= rz * forceTot;

                stress[0] -= forceTot * rx * rx;
                stress[1] -= forceTot * ry * ry;
                stress[2] -= forceTot * rz * rz;
                stress[3] -= forceTot * ry * rz;
                stress[4] -= forceTot * rx * rz;
                stress[5] -= forceTot * rx * ry;

            }

        } // end of loop over batom

        forceX[aatom] += fxi;
        forceY[aatom] += fyi;
        forceZ[aatom] += fzi;

        //std::cout << "\n force " << aatom+1 << " " << fxi << " " << fyi << " " << fzi;
 
    } // end of loop over aatom

    //exit(-1);

    // sum the energies
    if (numProcs > 1)
    {
        //manyEnergy = mpi.sumDoubleGroup(tmp_energy_r);
        twoenergy = mpi.sumDoubleGroup(tmp_energy_s);
    }
    else
    {
        //manyEnergy = tmp_energy_r;
        twoenergy = tmp_energy_s;
    }

    //std::cout << "\n two, many mody energy " << twoenergy << " " << manyEnergy << std::endl;

    
}

#endif

void Field::realSpaceForceImageNbrList(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                double* forceX, double* forceY, double* forceZ, double* stress, double* latVector,
                                double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        rank = mpi.getWorkGroupRank(),
        numProcs = mpi.getNumWorkGroupProcs(),
        ltypea,
        ltypeb;

    int
        batom,
        idx,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        fxi, fyi, fzi,
        eng_m,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    double
        rho;


    manyEnergy = 0.0;
    twoenergy = 0.0;

    for (aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
    {
        ltypea = atmLabel[aatom];
        wrkRho[aatom] = calculateAtomRho(aatom, posX, posY, posZ, latVector, rcpVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }

    if (numProcs > 1)
    {
        mpi.sumDoubleVectorGroup(wrkRho, numAtoms);
        manyEnergy = mpi.sumDoubleGroup(manyEnergy);
    }

    for (aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
    {

        ltypea = atmLabel[aatom];
        idx = nbrList.maxNbrs * aatom;

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        fxi = 0.0;
        fyi = 0.0;
        fzi = 0.0;
        
        for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
        {
        
            batom = nbrList.map[idx+nbr];
            ltypeb = atmLabel[batom];

            bx = posX[batom];
            by = posY[batom];
            bz = posZ[batom];

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
            yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
            zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            eng_m = 0.0;
            eng_s = 0.0;
            forc_m = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;
                
                mbdy.calculateManyBodyForce(ltypea, ltypeb, r, rsq, wrkRho[aatom], wrkRho[batom], forc_m);
                //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                //exit(-1);
                
                forceTot = (forc_s - forc_m) / (rsq * 2.0);
            
                tmp_energy_s += 0.5 * eng_s;

                fxi += rx * forceTot;
                fyi += ry * forceTot;
                fzi += rz * forceTot;
                forceX[batom] -= rx * forceTot;
                forceY[batom] -= ry * forceTot;
                forceZ[batom] -= rz * forceTot;

                stress[0] -= forceTot * rx * rx;
                stress[1] -= forceTot * ry * ry;
                stress[2] -= forceTot * rz * rz;
                stress[3] -= forceTot * ry * rz;
                stress[4] -= forceTot * rx * rz;
                stress[5] -= forceTot * rx * ry;

            }

        } // end of loop over batom

        forceX[aatom] += fxi;
        forceY[aatom] += fyi;
        forceZ[aatom] += fzi;

    } // end of loop over aatom

    if (numProcs > 1)
    {
        twoenergy = mpi.sumDoubleGroup(tmp_energy_s);
    }
    else
    {
        twoenergy = tmp_energy_s;
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions to calculate energies using repeting cell rather than image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Field::realSpaceForceNoImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                double* forceX, double* forceY, double* forceZ, double* stress, double* latVector,
                                double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        rank = mpi.getWorkGroupRank(),
        numProcs = mpi.getNumWorkGroupProcs(),
        ltypea,
        ltypeb;

    int
        batom,
        nx, ny, nz,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_s,
        tmp_energy_s = 0.0,
        forc_m,
        forc_s,
        forceTot;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoenergy = 0.0;

    for (aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
    {
        ltypea = atmLabel[aatom];
        wrkRho[aatom] = calculateAtomRhoNoImage(aatom, posX, posY, posZ, latVector, rcpVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }

    if (numProcs > 1)
    {
        mpi.sumDoubleVectorGroup(wrkRho, numAtoms);
        manyEnergy = mpi.sumDoubleGroup(manyEnergy);
    }

    for (aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
    {

        ltypea = atmLabel[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        for (batom = 0; batom < numAtoms; batom++)
        {

            ltypeb = atmLabel[batom];

            bx = posX[batom];
            by = posY[batom];
            bz = posZ[batom];
            
            for(nx = -cellX; nx <= cellX; nx++)
            {
                for(ny = -cellY; ny <= cellY; ny++)
                {
                    for(nz = -cellZ; nz <= cellZ; nz++)
                    {

                        //calculate distance
                        xx = bx + nx * latVector[0] + ny * latVector[3] + nz * latVector[6];
                        yy = by + nx * latVector[1] + ny * latVector[4] + nz * latVector[7];
                        zz = bz + nx * latVector[2] + ny * latVector[5] + nz * latVector[8];

                        rx = ax - xx;
                        ry = ay - yy;
                        rz = az - zz;

                        rsq = rx *rx + ry * ry + rz * rz;

                        if (rsq <= radius && rsq > 1.0e-3)
                        {

                            r = sqrt(rsq);

                            mbdy.calculateManyBodyPairForce(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                            //std::cout << "\n pair force " << aatom << " " << batom << " " << r  << " " << eng_s << " " << forc_s;
                
                            mbdy.calculateManyBodyForce(ltypea, ltypeb, r, rsq, wrkRho[aatom], wrkRho[batom], forc_m);
                            //std::cout << "\n nany body " << wrkRho[aatom] << " " <<  wrkRho[batom] << " " <<  forc_m;
                            //std::cout << "\n total " << batom+1 << "  " << r << "  " << (forc_s - forc_m) / rsq;
                            //exit(-1);
                
                            forceTot = (forc_s - forc_m) / (rsq * 2.0);
            
                            tmp_energy_s += 0.5 * eng_s;

                            forceX[aatom] += rx * forceTot;
                            forceY[aatom] += ry * forceTot;
                            forceZ[aatom] += rz * forceTot;
                            forceX[batom] -= rx * forceTot;
                            forceY[batom] -= ry * forceTot;
                            forceZ[batom] -= rz * forceTot;

                            stress[0] -= forceTot * rx * rx;
                            stress[1] -= forceTot * ry * ry;
                            stress[2] -= forceTot * rz * rz;
                            stress[3] -= forceTot * ry * rz;
                            stress[4] -= forceTot * rx * rz;
                            stress[5] -= forceTot * rx * ry;

                        }

                    } // nz

                } // ny

            } // nz

        } // end of loop over batom

    } // end of loop over aatom

    // sum the forces over processors here - not the best way but easiest for now
    if (numProcs > 1)
    {
        twoenergy = mpi.sumDoubleGroup(tmp_energy_s);
    }
    else
    {
        twoenergy = tmp_energy_s;
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// routines to calculate energy only i.e. MC
///////////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************
 driver for the calculation of  energy of all ions
 ************************************************************************/
void Field::calculateEnergy(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, double* atomCharge, 
                            int* atmLabel, int* frozen, int numAtoms, Energy& eng)
{

    eng.manyEnergy = 0.0,
    eng.vdwEnergy = 0.0;

    //mbdy.checkMesh();
    if (useNbrList)
    {
        realSpaceEnergyImageNbrList(mpi, eng.manyEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceEnergyImage(mpi, eng.manyEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                                 atomCharge, atmLabel, numAtoms);
    }

}

void Field::calculateAtomEnergy(const MPIComms& mpi, int aatom, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, 
                                double* atomCharge, int* atmLabel, int* frozen, int numAtoms, Energy& eng)
{

    eng.manyEnergy = 0.0,
    eng.vdwEnergy = 0.0;

    if (useNbrList)
    {
        realSpaceAtomEnergyImageNbrList(mpi, aatom, eng.manyEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceAtomEnergyImage(mpi, aatom, eng.manyEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }

}

void Field::calculateAtomEnergyRemove(const MPIComms& mpi, int aatom, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, 
                                double* atomCharge, int* atmLabel, int* frozen, int numAtoms, Energy& eng)
{

    eng.manyEnergy = 0.0,
    eng.vdwEnergy = 0.0;

    if (useNbrList)
    {
        realSpaceAtomEnergyImageNbrList(mpi, aatom, eng.manyEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceAtomEnergyRemoveImage(mpi, aatom, eng.manyEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }

}
void Field::calculateAtomEnergyDiff(const MPIComms& mpi, int aatom, double* oldPos, double* newPos, double* posX, double* posY, double* posZ, 
                            double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int* frozen, int numAtoms, 
                            Energy& engDiff)
{
    double
        oldMbdy = 0.0,
        newMbdy = 0.0,
        oldVdW = 0.0,
        newVdW = 0.0;

    engDiff.manyEnergy = 0.0;
    engDiff.vdwEnergy = 0.0;  

    //calculate original energy
    if (useNbrList)
    {
        realSpaceAtomEnergyImageNbrList(mpi, aatom, oldMbdy, oldVdW, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceAtomEnergyImage(mpi, aatom, oldMbdy, oldVdW, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }

    //make the move here and calculate new energy
    posX[aatom] = newPos[0];
    posY[aatom] = newPos[1];
    posZ[aatom] = newPos[2];

    if (useNbrList)
    {
        realSpaceAtomEnergyImageNbrList(mpi, aatom, newMbdy, newVdW, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceAtomEnergyImage(mpi, aatom, newMbdy, newVdW, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }

    //return energy difference
    engDiff.manyEnergy = newMbdy - oldMbdy;
    engDiff.vdwEnergy = newVdW - oldVdW;   
    
    // put old position back incase move is rejected
    posX[aatom] = oldPos[0];
    posY[aatom] = oldPos[1];
    posZ[aatom] = oldPos[2];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// routines using the nearest image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**************************************************************
 * calculate energy using image convention- this subroutine
 * uses explicit calculation of the energy and forces
 *
 * before doing the calculation it is necessary to
 * to reset the box so that any atom that has
 * migrated out of the box can be placed back in according
 * to the image convention
 **************************************************************/
#if defined  _OPENMP

void Field::realSpaceEnergyImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                 double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb,
        //i, j,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_p,
        tmp_energy_s = 0.0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoenergy = 0.0;

    for (aatom = 0; aatom < numAtoms; aatom++)
        wrkRho[aatom] = 0.0;

    //for (int aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        double rho = 0.0;
        
        for (batom = 0; batom < numAtoms; batom++)
        {

            if (aatom == batom)
                continue;

            ltypeb = atmLabel[batom];

            bx = posX[batom];
            by = posY[batom];
            bz = posZ[batom];

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
            yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
            zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                mbdy.calculateManyBodyPairEnergy(ltypea, ltypeb, r, rsq, eng_p);
                rho += mbdy.calculateManyBodyDensEnergy(ltypea, ltypeb, r, rsq);
                tmp_energy_s += eng_p * 0.5;
            }

        } // end of loop over batom

        wrkRho[aatom] = rho;

    } // end of loop over aatom

    twoenergy = tmp_energy_s;

    for (aatom = 0; aatom < numAtoms; aatom++)
    {
        ltypea = atmLabel[aatom];
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }      

    
    return;

}

#else
void Field::realSpaceEnergyImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                 double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb,
        //i, j,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_p,
        tmp_energy_s = 0.0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoenergy = 0.0;

    for (aatom = 0; aatom < numAtoms; aatom++)
        wrkRho[aatom] = 0.0;

    //for (int aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        double rho = 0.0;
        
        for (batom = 0; batom < numAtoms; batom++)
        {

            if (aatom == batom)
                continue;

            ltypeb = atmLabel[batom];

            bx = posX[batom];
            by = posY[batom];
            bz = posZ[batom];

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
            yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
            zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                mbdy.calculateManyBodyPairEnergy(ltypea, ltypeb, r, rsq, eng_p);
                rho += mbdy.calculateManyBodyDensEnergy(ltypea, ltypeb, r, rsq);
                tmp_energy_s += eng_p * 0.5;
            }

        } // end of loop over batom

        wrkRho[aatom] = rho;

    } // end of loop over aatom

    twoenergy = tmp_energy_s;

    for (aatom = 0; aatom < numAtoms; aatom++)
    {
        ltypea = atmLabel[aatom];
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    }      

    
    return;

}

#endif


void Field::realSpaceEnergyImageNbrList(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                        double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        batom,
        idx,
        ltypea,
        ltypeb,
        //i, j,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_p,
        tmp_energy_s = 0.0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoenergy = 0.0;

    for (aatom = 0; aatom < numAtoms; aatom++)
        wrkRho[aatom] = 0.0;

    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
    
        idx = nbrList.maxNbrs * aatom;

        double rho = 0.0;
        
        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];
        
        for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
        {

            batom = nbrList.map[idx+nbr];
            
            ltypeb = atmLabel[batom];

            bx = posX[batom];
            by = posY[batom];
            bz = posZ[batom];

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
            yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
            zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            if (rsq <= radius)
            {
                r = sqrt(rsq);
                mbdy.calculateManyBodyPairEnergy(ltypea, ltypeb, r, rsq, eng_p);
                tmp_energy_s += eng_p * 0.5;
                rho += mbdy.calculateManyBodyDensEnergy(ltypea, ltypeb, r, rsq);
            }

        } // end of loop over batom

        wrkRho[aatom] = rho;

    } // end of loop over aatom

    for (aatom = 0; aatom < numAtoms; aatom++)
    {
        ltypea = atmLabel[aatom];
        manyEnergy += mbdy.embed(ltypea, ltypea, wrkRho[aatom]);
    } 
    
    // sum the forces over processors here - not the best way but easiest for now
    twoenergy = tmp_energy_s;

    return;

}

void Field::realSpaceAtomEnergyImage(const MPIComms& mpi, int aatom, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                    double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_p,
        tmp_energy_s = 0.0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoenergy = 0.0;

    
    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (int i = 0; i < numAtoms; i++)
    {
        ltypea = atmLabel[i];
        wrkRho[i] = 0.0;
        wrkRho[i] = calculateAtomRho(i, posX, posY, posZ, latVector, rcpVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypeb, wrkRho[i]);
    }

    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    ltypea = atmLabel[aatom];

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];


    for (batom = 0; batom < numAtoms; batom++)
    {

        if (batom == aatom)
            continue;

        ltypeb = atmLabel[batom];

        bx = posX[batom];
        by = posY[batom];
        bz = posZ[batom];

        rx = ax - bx;
        ry = ay - by;
        rz = az - bz;

        xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
        yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
        zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        
        if (rsq <= radius)
        {
            r = sqrt(rsq);

            eng_p = 0.0;
            mbdy.calculateManyBodyPairEnergy(ltypea, ltypeb, r, rsq, eng_p);
            tmp_energy_s += eng_p;
        }

    } // end of loop over batom
    
    twoenergy = tmp_energy_s;

}

void Field::realSpaceAtomEnergyImageNbrList(const MPIComms& mpi, int aatom,  double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                            double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        batom,
        idx,
        ltypea,
        ltypeb;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_p,
        tmp_energy_s = 0.0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoenergy = 0.0;

   for (int i = 0; i < numAtoms; i++)
    {
        ltypea = atmLabel[i];
        wrkRho[i] = 0.0;
        wrkRho[i] = calculateAtomRho(i, posX, posY, posZ, latVector, rcpVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypeb, wrkRho[i]);
    }
    
    ltypea = atmLabel[aatom];
    idx = nbrList.maxNbrs * aatom;

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];

    for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
    {
    
        batom = nbrList.map[idx+nbr];
        ltypeb = atmLabel[batom];

        bx = posX[batom];
        by = posY[batom];
        bz = posZ[batom];

        rx = ax - bx;
        ry = ay - by;
        rz = az - bz;

        xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
        yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
        zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        eng_p = 0.0;
    
        if (rsq <= radius)
        {

            r = sqrt(rsq);

            mbdy.calculateManyBodyPairEnergy(ltypea, ltypeb, r, rsq, eng_p);

        }

        tmp_energy_s += eng_p;

    } // end of loop over batom

    twoenergy = tmp_energy_s;

}

void Field::realSpaceAtomEnergyRemoveImage(const MPIComms& mpi, int aatom, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                    double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        batom,
        ltypea,
        ltypeb;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        eng_p,
        tmp_energy_s = 0.0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    manyEnergy = 0.0;
    twoenergy = 0.0;

    
    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (int i = 0; i < numAtoms; i++)
    {
        if (i == aatom)
            continue;
            
        ltypea = atmLabel[i];
        wrkRho[i] = 0.0;
        wrkRho[i] = calculateAtomRho(i, posX, posY, posZ, latVector, rcpVector, atmLabel, numAtoms);
        manyEnergy += mbdy.embed(ltypea, ltypeb, wrkRho[i]);
    }

    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    ltypea = atmLabel[aatom];

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];


    for (batom = 0; batom < numAtoms; batom++)
    {

        if (batom == aatom)
            continue;

        ltypeb = atmLabel[batom];

        bx = posX[batom];
        by = posY[batom];
        bz = posZ[batom];

        rx = ax - bx;
        ry = ay - by;
        rz = az - bz;

        xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
        yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
        zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        
        if (rsq <= radius)
        {
            r = sqrt(rsq);

            eng_p = 0.0;
            mbdy.calculateManyBodyPairEnergy(ltypea, ltypeb, r, rsq, eng_p);
            tmp_energy_s += eng_p;
        }

    } // end of loop over batom
    
    twoenergy = tmp_energy_s;

}
////////////////////////////////////////////////////////////////
// dummy routines not used by potentials
/////////////////////////////////////////////////////////////////
void Field::setExternalCommunicator(const MPIComms& mpi, std::ofstream& outstream)
{

}

void Field::createFileSystem(const MPIComms& mpi, std::string rootName, std::string externalFile, bool restart, std::ofstream& outStream)
{

}
void Field::calculateConfig(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, double* atomCharge, 
                         int* atmLabel, int* frozen, int numAtoms, Energy& eng)
{
    
    eng.totalEnergy = 0.0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// utility functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Field::sumForces(const MPIComms& mpi, double* forceX, double* forceY, double* forceZ, int n)
{

    mpi.sumDoubleVectorGroup(forceX, n);
    mpi.sumDoubleVectorGroup(forceY, n);
    mpi.sumDoubleVectorGroup(forceZ, n);
}

double Field::calculateAtomRhoNoImage(int aatom, double* posX, double* posY, double* posZ, 
                                 double* latVector, double* rcpVector, int* atmLabel, int numAtoms)

{
    int
        batom,
        ltypea,
        ltypeb,
        nx, ny, nz;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        engMany = 0.0;

    double
        rho = 0.0,
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    ltypea = atmLabel[aatom];

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];

    for (batom = 0; batom < numAtoms; batom++)
    {
    
        ltypeb = atmLabel[batom];

        bx = posX[batom];
        by = posY[batom];
        bz = posZ[batom];

        for(nx = -cellX; nx <= cellX; nx++)
        {
            for(ny = -cellY; ny <= cellY; ny++)
            {
                for(nz = -cellZ; nz <= cellZ; nz++)
                {
                    
                    //calculate distance
                    xx = bx + nx * latVector[0] + ny * latVector[3] + nz * latVector[6];
                    yy = by + nx * latVector[1] + ny * latVector[4] + nz * latVector[7];
                    zz = bz + nx * latVector[2] + ny * latVector[5] + nz * latVector[8];

                    rx = ax - xx;
                    ry = ay - yy;
                    rz = az - zz;

                    rsq = rx *rx + ry * ry + rz * rz;

                    if (rsq <= radius && rsq > 1.0e-3)
                    {
                        r = sqrt(rsq);

                        rho += mbdy.calculateManyBodyDensEnergy(ltypea, ltypeb, r, rsq);
                    }

                }

            }
        }

    }
    
    return rho;
}
double Field::calculateAtomRho(int aatom, double* posX, double* posY, double* posZ, 
                                 double* latVector, double* rcpVector, int* atmLabel, int numAtoms)

{
    int
        batom,
        ltypea,
        ltypeb;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        engMany = 0.0;

    double
        rho = 0.0,
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    ltypea = atmLabel[aatom];

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];

    for (batom = 0; batom < numAtoms; batom++)
    {
        if (aatom == batom)
            continue;

        ltypeb = atmLabel[batom];

        bx = posX[batom];
        by = posY[batom];
        bz = posZ[batom];

        rx = ax - bx;
        ry = ay - by;
        rz = az - bz;

        xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
        yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
        zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        if (rsq <= radius)
        {
            r = sqrt(rsq);

            rho += mbdy.calculateManyBodyDensEnergy(ltypea, ltypeb, r, rsq);
        } 

    } // end of loop over aatom

    return rho;

}

void Field::resetSimulationCellNoImage(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int natoms)
{

    double
        xx, yy, zz,
        rx, ry, rz;

    for (int aatom = 0; aatom < natoms; aatom++)
    {

        rx = posX[aatom];
        ry = posY[aatom];
        rz = posZ[aatom];

        xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
        yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
        zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

        //reduced coordinates should lie between 0 and 1
        if (xx < 0.0)
            xx += 1.0;
        if (yy < 0.0)
            yy += 1.0;
        if (zz < 0.0)
            zz += 1.0;
        if (xx >= 1.0)
            xx -= 1.0;
        if (yy >= 1.0)
            yy -= 1.0;
        if (zz >= 1.0)
            zz -= 1.0;

        posX[aatom] = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        posY[aatom] = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        posZ[aatom] = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

    }
}

void Field::resetSimulationCellMinImage(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int numAtoms)
{
    double
        rx, ry, rz,
        xx, yy, zz;

    for (int i = 0; i < numAtoms; i++)
    {

        rx = posX[i];
        ry = posY[i];
        rz = posZ[i];

        xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
        yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
        zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        posX[i] = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        posY[i] = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        posZ[i] = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

    }
}

void Field::resetAtomSimulationCell(int i, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int numAtoms)
{
    double
        rx, ry, rz,
        xx, yy, zz;

    rx = posX[i];
    ry = posY[i];
    rz = posZ[i];

    xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
    yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
    zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

    xx -= rint(xx);
    yy -= rint(yy);
    zz -= rint(zz);

    posX[i] = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
    posY[i] = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
    posZ[i] = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

}

// dummy routines
void Field::tidy(const MPIComms& mpi)
{

}
void Field::copyWorkFiles(const MPIComms& mpi)
{

}

void Field::updatePhaseSum(void)
{
    for (int i = 0; i <  wrkSize; i++)
        tmpRho[i] = wrkRho[i];
}

void Field::resetAtomMove(int atm)
{
    //update th neighborlist if required
}
