#include "Field.h"
void dcell(double*, double[10]);
std::vector<std::string> split(std::string s);

Field::Field()
{
    newjob = true;
}

Field::~Field()
{
    free_dvector(erfc_e);
    free_dvector(erfc_f);
}

/****************************************************************************************
 create error function lookup table
****************************************************************************************/
void Field::fieldErfcGen(double alpha)
{

    int
        i;

    double
        tt, exp1,
        rrr, rsq;

    double
        //      epsq = job.dielec,
        alph = alpha,
        a1 = 0.254829592,
        a2 = -0.284496736,
        a3 = 1.421413741,
        a4 = -1.453152027,
        a5 = 1.061405429,
        pp = 0.3275911;

    if (alpha == 0.0)
        return;

    erfc_e = alloc_dvector(MAXMESH, "nonbonded erfc", 0.0);
    erfc_f = alloc_dvector(MAXMESH, "nonbonded erfc", 0.0);

    delMesh = shortRangeCut / (MAXMESH - 4);      // spacing between mesh points
    rMesh = (MAXMESH - 4) / shortRangeCut;

    for (i = 1; i < MAXMESH; i++)
    {

        rrr = double(i) * delMesh;
        rsq = rrr * rrr;

        tt = 1.0 / (1.0 + pp * alph * rrr);
        exp1 = exp(-pow((alph * rrr), 2.0));

        erfc_e[i] = tt * (a1 + tt * (a2 + tt * (a3 + tt * (a4 + tt * a5)))) * exp1 / rrr;
        erfc_f[i] = (erfc_e[i] + 2.0 * (alph / ROOTPI) * exp1) / rsq;
    }

    aa = erfc_f[MAXMESH - 4] * shortRangeCut;
    bb = -(erfc_e[MAXMESH - 4] + aa * shortRangeCut);

    //  b0 = 2.0 * (epsq - 1.0) / (2.0 * epsq + 1.0);
    //  rfld0 = b0 / pow(job.shortrangecut,3);
    //  rfld1 = (1.0 + 0.5 * b0) / job.shortrangecut;
    rfld0 = erfc_e[MAXMESH - 4];
    rfld1 = shortRangeCut * erfc_f[MAXMESH - 4];

    //cout << setiosflags(ios::showpoint | ios::fixed | ios::right)
    //   << setprecision(15) << "\n gen erf mesh " << job.shortrangecut << " " << rfld0  << " " << rfld1;
    //cout.flush();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// utility functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Field::setup(const MPIComms& mpi, double* latVector, double* rcpVector, double minDimension, double volume, 
                  int numAtoms, const Species& spec, std::ofstream& outStream)
{


    double
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

    // check the image convention is obeyed
    if (minimumImage && (shortRangeCut + verletShell) > (minDimension / 2))
    {
        outStream << "\n\n\n *** short range shortRangeCut to large for simulation box." << std::endl;
        outStream << " cutoff : " << shortRangeCut << " minimum dimension: " << minDimension << std::endl;
        outStream.flush();
        mpi.commsAbortWorld();
        exit(EXIT_FAILURE);
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

    if (doEwald > 0)
    {

        if (doEwald == 2)
        {
            eta = alpha;
            fieldErfcGen(eta);
        }
        else
        {
            eps = fmin(fabs(madelungAcc), 0.5);
            tol = sqrt(fabs(log(eps * shortRangeCut)));

            eta = sqrt(fabs(log(eps * shortRangeCut * tol))) / shortRangeCut;

        }

    }
    

    if (newjob == true)
    {
        outStream << "\n cutoffs for lattice summation :"
                  << "\n real space shortRangeCut       = " << shortRangeCut << " Angstroms" << std::endl;

        outStream << "\n cutoffs for lattice summation eta      = " << eta << std::endl;
    }

    if (newjob == true && doEwald  == 1)
    {

        coul.setEwaldImage(eta, shortRangeCut, eps, tol, latVector, rcpVector, outStream);
        coul.estimateGvectorvoid(mpi, volume, rcpVector);
        coul.allocateArrays();
        // calculate g vectors
        coul.findGvector(mpi, volume, rcpVector, false);

        if (newjob == true)
        {
            coul.printEwald(mpi, outStream);
            coul.printNumGVec(mpi, outStream);
        }

    }

    // set up short-range potential mesh
    if (newjob == true) vdw.calculateEnergyMesh(spec, shortRangeCut, outStream);

    newjob = false;

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

	if (doEwald == 1)
	{
		outStream << " the accuracy for the Ewald sum is " << madelungAcc << std::endl;
	}
	else if (doEwald == 2)
	{
		outStream << " the shifted and damped Coulomb sum will be used with an alpha of " << alpha << std::endl;
	}
    else if (doEwald == 0) 
	{
		outStream << " no Coulomb interactions will be evaluated" << std::endl;
	}
    
    if (doEwald == 1)
    {

        coul.printEwald(mpi, outStream);
        coul.printNumGVec(mpi, outStream);
        
    }

    vdw.printPotential(outStream);
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

    // read in and print out the species and potentails
    std::vector<std::string> words;

    while (!inStream.eof())
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
        else if (keyWord == "ewald")
		{
			subWord = words[1];

			if (subWord == "precis") 
			{
		        doEwald = 1;
				dummy = words[2];
                madelungAcc = std::stod(dummy);
			}
			else if (subWord == "damp")
			{
		        doEwald = 2;
				dummy = words[2];
                alpha = std::stod(dummy);
			}
			else if (subWord == "none") 
			{
				doEwald = 0;
			}

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
        else if (keyWord == "twobody")
        {

            dummy = words[1];
            num = std::stoi(dummy);
            cFact = CTOEV;
            

            vdw.loadPotential(inStream, outStream, num);
            vdw.convertPotential(energyUnit);
            vdw.printPotential(outStream);

        }
        else if (keyWord == "close")
        {
            break;
        }
        else
        {
            outStream << "\n unrecognised key word used in potentials file " << keyWord << "\n from line " << line << std::endl;
            outStream.flush();
            exit(EXIT_FAILURE);
            //throw std::invalid_argument("Invalid keyword!");
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
        realenergy = 0.0,
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


    if (doEwald == 1)   //calculate reciprocal space part of ewald
    {

        rcpenergy = coul.recipSpaceForce(mpi, posX, posY, posZ, forceX, forceY, forceZ, stress, atomCharge, numAtoms);

    }
    if (minimumImage)
    {
        resetSimulationCellMinImage(mpi, posX, posY, posZ, latVector, rcpVector, numAtoms);

        if (useNbrList)
        {
            realSpaceForceImageNbrList(mpi, realenergy, twoenergy, posX, posY, posZ, forceX, forceY, forceZ, stress, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
        }
        else
        {
            realSpaceForceImage(mpi, realenergy, twoenergy, posX, posY, posZ, forceX, forceY, forceZ, stress, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
        }
               
    }
    else
    {
            resetSimulationCellNoImage(mpi, posX, posY, posZ, latVector, rcpVector, numAtoms);
            realSpaceForceNoImage(mpi, realenergy, twoenergy, posX, posY, posZ, forceX, forceY, forceZ, stress, latVector, rcpVector, 
                                  atomCharge, atmLabel, numAtoms);
    }

    
    eng.rcpEnergy = rcpenergy;
    eng.realEnergy = realenergy;
    eng.vdwEnergy = twoenergy;
    
    eng.totalEnergy = rcpenergy + realenergy + twoenergy;

    if (mpi.getNumWorkGroupProcs() > 1)
        sumForces(mpi, forceX, forceY, forceZ, numAtoms);

    
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

void Field::realSpaceForceImage(const MPIComms& mpi, double& realenergy, double& twoenergy, double* __restrict__ posX, double* __restrict__ posY, double* __restrict__ posZ, 
                                double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__ latVector,
                                double* __restrict__ rcpVector, double* __restrict__ atomCharge, int* __restrict__ atmLabel, int numAtoms)
{
    int
        ltypea,
        ltypeb;

    int
        batom,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        h0,
        h1,
        eng,
        eng_s,
        tmp_self = 0.0,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = hfct0 * hfct0 * hfct0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_self += -h0 * chargeprod * hfct0;

        }

    }

    double* tmpForceX; //[numAtoms];
    double* tmpForceY; //[numAtoms];
    double* tmpForceZ; //[numAtoms];

    //for (aatom = mpi.getAtomLoopStart(); aatom < mpi.getAtomLoopFinish(); aatom++)
    #pragma omp parallel default(none)  \
        reduction(+: tmp_energy_r, tmp_energy_s) \
        shared(atmLabel, atomCharge, posX, posY, posZ, numAtoms, rcpVector, latVector, forceX, forceY, forceZ, radius, shortRangeCut, hfct0, hfct1, doEwald) \
        private(ltypea, ltypeb, chargea, chargeb, chargeprod, ax, ay, az, bx, by, bz, rx, ry, rz, r, rsq, xx, yy, zz, aatom, batom) \
        private(h0, h1, eng, eng_s, forc, forc_s, tmpForceX, tmpForceY, tmpForceZ)
    {

        tmpForceX = (double*) calloc(numAtoms, sizeof(double));
        tmpForceY = (double*) calloc(numAtoms, sizeof(double));
        tmpForceZ = (double*) calloc(numAtoms, sizeof(double));

        #pragma omp for schedule(static,4)
        for (aatom = 0; aatom < numAtoms; aatom++)
        {

            ltypea = atmLabel[aatom];
            chargea = cFact * atomCharge[aatom];

            ax = posX[aatom];
            ay = posY[aatom];
            az = posZ[aatom];


            for (batom = aatom + 1; batom < numAtoms; batom++)
            {

                chargeb = atomCharge[batom];
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

                eng = 0.0;
                eng_s = 0.0;
                forc = 0.0;
                forc_s = 0.0;

                if (rsq <= radius)
                {

                    r = sqrt(rsq);

                    h0 = 0.0;
                    h1 = 0.0;

                    chargeprod = chargea * chargeb;

                    if (doEwald == 1)
                    {
                        hfunc1(hfct0 * r, h0, h1);

                        eng = h0 * chargeprod * hfct0;
                        forc = h1 * chargeprod * hfct1;
                    }
                    else if (doEwald == 2)
                    {
                        //evaluate damped and shifted
                        interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                    }

                    if (r <= shortRangeCut && r > 1.0e-4)
                    {
                        vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                    }

                }

                tmp_energy_r += eng;
                tmp_energy_s += eng_s;

                forc -= forc_s;

                tmpForceX[aatom] -= rx * forc;
                tmpForceY[aatom] -= ry * forc;
                tmpForceZ[aatom] -= rz * forc;
                tmpForceX[batom] += rx * forc;
                tmpForceY[batom] += ry * forc;
                tmpForceZ[batom] += rz * forc;

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
        free(tmpForceX);
        free(tmpForceY);
        free(tmpForceZ);

    }  // end of parallel

    //for (aatom = 0; aatom < numAtoms; aatom++)
    //{
    //    std::cout << "\n force rreal " << aatom << " " << forceX[aatom] << " " << forceY[aatom] << " " << forceZ[aatom];
    //}
    realenergy = tmp_energy_r + tmp_self;
    twoenergy = tmp_energy_s;

}

#elif defined THREADS
void Field::realSpaceForceImage(const MPIComms& mpi, double& realenergy, double& twoenergy, double* __restrict__ posX, double* __restrict__ posY, double* __restrict__ posZ, 
                                double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__ latVector,
                                double* __restrict__ rcpVector, double* __restrict__ atomCharge, int* __restrict__ atmLabel, int numAtoms)
{
    int
        ltypea,
        ltypeb;

    int
        batom,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        h0,
        h1,
        eng,
        eng_s,
        tmp_self = 0.0,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = hfct0 * hfct0 * hfct0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_self += -h0 * chargeprod * hfct0;

        }

    }

}
#else
void Field::realSpaceForceImage(const MPIComms& mpi, double& realenergy, double& twoenergy, double* __restrict__ posX, double* __restrict__ posY, double* __restrict__ posZ, 
                                double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ, double* __restrict__ stress, double* __restrict__ latVector,
                                double* __restrict__ rcpVector, double* __restrict__ atomCharge, int* __restrict__ atmLabel, int numAtoms)
{
    int
        rank = mpi.getWorkGroupRank(),
        numProcs = mpi.getNumWorkGroupProcs(),
        batom,
        ltypea,
        ltypeb,
        //i, j,
        aatom;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        fxi, fyi, fzi,
        h0,
        h1,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = hfct0 * hfct0 * hfct0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        //for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            //chargeprod = cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    for (aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
    //for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        fxi = 0.0;
        fyi = 0.0;
        fzi = 0.0;

        for (batom = aatom; batom < numAtoms; batom++)
        {
            if (aatom == batom)
                continue;

            chargeb = atomCharge[batom];
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

            eng = 0.0;
            eng_s = 0.0;
            forc = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                h0 = 0.0;
                h1 = 0.0;

                chargeprod = chargea * chargeb;

                if (doEwald == 1)
                {
                    hfunc1(hfct0 * r, h0, h1);

                    eng = h0 * chargeprod * hfct0;
                    forc = h1 * chargeprod * hfct1;
                }
                else if (doEwald == 2)
                {
                    //evaluate damped and shifted
                    interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                }

                //if (r <= shortRangeCut)
                //{
                    vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                //}

            }

            tmp_energy_r += eng;
            tmp_energy_s += eng_s;

            forc -= forc_s;

            fxi -= rx * forc;
            fyi -= ry * forc;
            fzi -= rz * forc;
            forceX[batom] += rx * forc;
            forceY[batom] += ry * forc;
            forceZ[batom] += rz * forc;

            stress[0] += forc * rx * rx;
            stress[1] += forc * ry * ry;
            stress[2] += forc * rz * rz;
            stress[3] += forc * ry * rz;
            stress[4] += forc * rx * rz;
            stress[5] += forc * rx * ry;

        } // end of loop over batom

        forceX[aatom] += fxi;
        forceY[aatom] += fyi;
        forceZ[aatom] += fzi;

    } // end of loop over aatom

    // sum the energies
    if (numProcs > 1)
    {
        realenergy = mpi.sumDoubleGroup(tmp_energy_r);
        twoenergy = mpi.sumDoubleGroup(tmp_energy_s);
    }
    else
    {
        realenergy = tmp_energy_r;
        twoenergy = tmp_energy_s;
    }

}

#endif


void Field::realSpaceForceImageNbrList(const MPIComms& mpi, double& realenergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                double* forceX, double* forceY, double* forceZ, double* stress, double* latVector,
                                double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
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
        h0,
        h1,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = pow(hfct0, 3.0);

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;


    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    for (int aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];
        idx = nbrList.maxNbrs * aatom;

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];
        
        for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
        {
        
            batom = nbrList.map[idx+nbr];
            chargeb = atomCharge[batom];
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

            eng = 0.0;
            eng_s = 0.0;
            forc = 0.0;
            forc_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                h0 = 0.0;
                h1 = 0.0;

                chargeprod = chargea * chargeb;

                if (doEwald == 1)
                {
                    hfunc1(hfct0 * r, h0, h1);

                    eng = h0 * chargeprod * hfct0;
                    forc = h1 * chargeprod * hfct1;
                }
                else if (doEwald == 2)
                {
                    //evaluate damped and shifted
                    interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                }

                if (r <= shortRangeCut && r > 1.0e-4)
                {
                    vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);
                }

            }

            tmp_energy_r += eng;
            tmp_energy_s += eng_s;

            forc -= forc_s;

            forceX[aatom] -= rx * forc;
            forceY[aatom] -= ry * forc;
            forceZ[aatom] -= rz * forc;
            forceX[batom] += rx * forc;
            forceY[batom] += ry * forc;
            forceZ[batom] += rz * forc;

            stress[0] += forc * rx * rx;
            stress[1] += forc * ry * ry;
            stress[2] += forc * rz * rz;
            stress[3] += forc * ry * rz;
            stress[4] += forc * rx * rz;
            stress[5] += forc * rx * ry;

        } // end of loop over batom

    } // end of loop over aatom

    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions to calculate energies using repeting cell rather than image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Field::realSpaceForceNoImage(const MPIComms& mpi, double& realenergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                double* forceX, double* forceY, double* forceZ, double* stress, double* latVector,
                                double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms)
{
    int
        numProcs = 1,
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
        h0,
        h1,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        forc,
        forc_s,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = pow(hfct0, 3.0);

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    bool
        sameCell = false;

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0, h1);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        for (batom = aatom; batom < numAtoms; batom++)
        {

            chargeb = atomCharge[batom];
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
                        sameCell = false;
                        if (nx == 0 && ny == 0 && nz == 0)
                            sameCell = true;

                        if (aatom == batom && sameCell == true)
                            continue;

                        //calculate distance
                        xx = bx + nx * latVector[0] + ny * latVector[3] + nz * latVector[6];
                        yy = by + nx * latVector[1] + ny * latVector[4] + nz * latVector[7];
                        zz = bz + nx * latVector[2] + ny * latVector[5] + nz * latVector[8];

                        rx = ax - xx;
                        ry = ay - yy;
                        rz = az - zz;

                        rsq = rx *rx + ry * ry + rz * rz;

                        eng = 0.0;
                        eng_s = 0.0;
                        forc = 0.0;
                        forc_s = 0.0;

                        if (rsq <= radius)
                        {

                            r = sqrt(rsq);

                            h0 = 0.0;
                            h1 = 0.0;

                            chargeprod = chargea * chargeb;

                            if (doEwald == 1)
                            {
                                hfunc1(hfct0 * r, h0, h1);

                                eng = h0 * chargeprod * hfct0;
                                forc = h1 * chargeprod * hfct1;
                            }
                            else if (doEwald == 2)
                            {
                                //evaluate damped and shifted
                                interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng, forc);

                            }

                            vdw.shortrangeone(ltypea, ltypeb, r, rsq, eng_s, forc_s);

                            tmp_energy_r += eng;
                            tmp_energy_s += eng_s;

                            forc -= forc_s;

                            forceX[aatom] -= rx * forc;
                            forceY[aatom] -= ry * forc;
                            forceZ[aatom] -= rz * forc;
                            forceX[batom] += rx * forc;
                            forceY[batom] += ry * forc;
                            forceZ[batom] += rz * forc;

                            stress[0] += forc * rx * rx;
                            stress[1] += forc * ry * ry;
                            stress[2] += forc * rz * rz;
                            stress[3] += forc * ry * rz;
                            stress[4] += forc * rx * rz;
                            stress[5] += forc * rx * ry;
                        }

                    } // nz

                } // ny

            } // nz

        } // end of loop over batom

    } // end of loop over aatom

    // sum the forces over processors here - not the best way but easiest for now
    if (numProcs > 1)
    {
        realenergy = mpi.sumDoubleGroup(tmp_energy_r);
        twoenergy = mpi.sumDoubleGroup(tmp_energy_s);
    }
    else
    {
        realenergy = tmp_energy_r;
        twoenergy = tmp_energy_s;
    }
    //MPI_Allreduce(&tmp_energy_r, &realenergy, 1, MPI_DOUBLE, MPI_SUM, mpi.getCommunicator());
    //MPI_Allreduce(&tmp_energy_s, &twoenergy, 1, MPI_DOUBLE, MPI_SUM, mpi.getCommunicator());


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

    if (doEwald == 1)
        eng.rcpEnergy = coul.recipSpaceEnergy(mpi, posX, posY, posZ, atomCharge, numAtoms);

    if (useNbrList)
    {
        realSpaceEnergyImageNbrList(mpi, eng.realEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                                atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceEnergyImage(mpi, eng.realEnergy, eng.vdwEnergy, posX, posY, posZ, latVector, rcpVector, 
                             atomCharge, atmLabel, numAtoms);
    }
        


}

void Field::calculateAtomEnergy(const MPIComms& mpi, int aatom, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, 
                                double* atomCharge, int* atmLabel, int* frozen, int numAtoms, Energy& eng)
{

    eng.rcpEnergy = 0.0,
    eng.realEnergy = 0.0,
    eng.vdwEnergy = 0.0;
    eng.selfEnergy = 0.0;

    //Coulomb interaction 

    if (useNbrList)
    {
        realSpaceAtomEnergyImageNbrList(mpi, aatom, eng.realEnergy, eng.vdwEnergy, eng.selfEnergy, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceAtomEnergyImage(mpi, aatom, eng.realEnergy, eng.vdwEnergy, eng.selfEnergy, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
        
}

void Field::calculateAtomEnergyDiff(const MPIComms& mpi, int aatom, double* oldPos, double* newPos, double* posX, double* posY, double* posZ, 
                            double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int* frozen, int numAtoms, 
                            Energy& engDiff)
{
    double
        oldReal = 0.0,
        newReal = 0.0,
        oldVdW = 0.0,
        newVdW = 0.0, 
        oldSelf = 0.0,
        newSelf = 0.0;

    //calculate original energy
    if (useNbrList)
    {
        realSpaceAtomEnergyImageNbrList(mpi, aatom, oldReal, oldVdW, oldSelf, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceAtomEnergyImage(mpi, aatom, oldReal, oldVdW, oldSelf, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }

    //make the move here and calculate new energy
    posX[aatom] = newPos[0];
    posY[aatom] = newPos[1];
    posZ[aatom] = newPos[2];

    if (useNbrList)
    {
        realSpaceAtomEnergyImageNbrList(mpi, aatom, newReal, newVdW, newSelf, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }
    else
    {
        realSpaceAtomEnergyImage(mpi, aatom, newReal, newVdW, newSelf, posX, posY, posZ, latVector, rcpVector, 
                                    atomCharge, atmLabel, numAtoms);
    }

    //return energy difference
    engDiff.realEnergy = newReal - oldReal;
    engDiff.vdwEnergy = newVdW - oldVdW;   
    engDiff.selfEnergy = newSelf - oldSelf;    
    
    // put old position back incase move is rejected
    posX[aatom] = oldPos[0];
    posY[aatom] = oldPos[1];
    posZ[aatom] = oldPos[2];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// routines using the nearest image convention
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#if defined THREADS
void Field::realSpaceEnergyImage(double& realenergy, double& twoenergy, double* posX, double* posY, double* posZ, 
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
        h0,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = hfct0 * hfct0 * hfct0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    //for (int aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];


        for (batom = aatom + 1; batom < numAtoms; batom++)
        {

            chargeb = atomCharge[batom];
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

            eng = 0.0;
            eng_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                h0 = 0.0;

                chargeprod = chargea * chargeb;

                if (doEwald == 1)
                {
                    hfunc1(hfct0 * r, h0);

                    eng = h0 * chargeprod * hfct0;

                }
                else if (doEwald == 2)
                {
                    //evaluate damped and shifted
                    interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng);

                }

                if (r <= shortRangeCut && r > 1.0e-4)
                {
                    vdw.shortRange(ltypea, ltypeb, r, rsq, eng_s);
                }

            }

            tmp_energy_r += eng;
            tmp_energy_s += eng_s;

        } // end of loop over batom

    } // end of loop over aatom

    
    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;
    
    return;

}

#else
/**************************************************************
 * calculate energy using image convention- this subroutine
 * uses explicit calculation of the energy and forces
 *
 * before doing the calculation it is necessary to
 * to reset the box so that any atom that has
 * migrated out of the box can be placed back in according
 * to the image convention
 **************************************************************/
void Field::realSpaceEnergyImage(const MPIComms& mpi, double& realenergy, double& twoenergy, double* posX, double* posY, double* posZ, 
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
        h0,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = hfct0 * hfct0 * hfct0;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    //for (int aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];


        for (batom = aatom + 1; batom < numAtoms; batom++)
        {

            chargeb = atomCharge[batom];
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

            eng = 0.0;
            eng_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                h0 = 0.0;

                chargeprod = chargea * chargeb;

                if (doEwald == 1)
                {
                    hfunc1(hfct0 * r, h0);

                    eng = h0 * chargeprod * hfct0;

                }
                else if (doEwald == 2)
                {
                    //evaluate damped and shifted
                    interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng);

                }

                if (r <= shortRangeCut && r > 1.0e-4)
                {
                    vdw.shortRange(ltypea, ltypeb, r, rsq, eng_s);
                }

            }

            tmp_energy_r += eng;
            tmp_energy_s += eng_s;

        } // end of loop over batom

    } // end of loop over aatom

    
    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;
    
    return;

}

#endif
void Field::realSpaceEnergyImageNbrList(const MPIComms& mpi, double& realenergy, double& twoenergy, double* posX, double* posY, double* posZ, 
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
        h0,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = pow(hfct0, 3.0);

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;

    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            tmp_energy_r += -h0 * chargeprod * hfct0;

        }

    }

    //std::cout << "\n FIELD PBC nbrlist ewald: " << doEwald;

    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < numAtoms; aatom++)
    {

        ltypea = atmLabel[aatom];
        chargea = cFact * atomCharge[aatom];
        idx = nbrList.maxNbrs * aatom;

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];
        //std::cout << "\n FIELD PBC nbrlist  " << nbrList.numnbrs[aatom];
        for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
        {

            batom = nbrList.map[idx+nbr];
            chargeb = atomCharge[batom];
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

            eng = 0.0;
            eng_s = 0.0;

            if (rsq <= radius)
            {

                r = sqrt(rsq);

                h0 = 0.0;

                chargeprod = chargea * chargeb;

                if (doEwald == 1)
                {
                    hfunc1(hfct0 * r, h0);

                    eng = h0 * chargeprod * hfct0;

                }
                else if (doEwald == 2)
                {
                    //evaluate damped and shifted
                    interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng);

                }

                if (r <= shortRangeCut && r > 1.0e-4)
                {
                    vdw.shortRange(ltypea, ltypeb, r, rsq, eng_s);
                }

            }

            tmp_energy_r += eng;
            tmp_energy_s += eng_s;

        } // end of loop over batom

    } // end of loop over aatom

    // sum the forces over processors here - not the best way but easiest for now
    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;

    return;

}

void Field::realSpaceAtomEnergyImage(const MPIComms& mpi, int aatom, double& realenergy, double& twoenergy, double& selfEnergy, double* posX, double* posY, double* posZ, 
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
        h0,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta;

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;
    selfEnergy = 0.0;

    if (doEwald == 1)
    {

        hfunc2(0.0, h0);

        chargea = atomCharge[aatom];
        chargeprod = 0.5 * cFact * chargea * chargea;
        selfEnergy += -h0 * chargeprod * hfct0;

    }

    //for (aatom = rank; aatom < numAtoms - 1; aatom = aatom + numProcs)
    ltypea = atmLabel[aatom];
    chargea = cFact * atomCharge[aatom];

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];


    for (batom = 0; batom < numAtoms; batom++)
    {

        if (batom == aatom)
            continue;

        chargeb = atomCharge[batom];
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

        eng = 0.0;
        eng_s = 0.0;

        if (rsq <= radius)
        {

            r = sqrt(rsq);

            h0 = 0.0;

            chargeprod = chargea * chargeb;

            if (doEwald == 1)
            {
                hfunc1(hfct0 * r, h0);

                eng = h0 * chargeprod * hfct0;

            }
            else if (doEwald == 2)
            {
                //evaluate damped and shifted
                interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng);

            }

            if (r <= shortRangeCut && r > 1.0e-4)
            {
                vdw.shortRange(ltypea, ltypeb, r, rsq, eng_s);
            }

        }

        tmp_energy_r += eng;
        tmp_energy_s += eng_s;

    } // end of loop over batom
    
    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;

}

void Field::realSpaceAtomEnergyImageNbrList(const MPIComms& mpi, int aatom,  double& realenergy, double& twoenergy, double& selfEnergy, double* posX, double* posY, double* posZ, 
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
        h0,
        eng,
        eng_s,
        tmp_energy_r = 0.0,
        tmp_energy_s = 0.0,
        chargea,
        chargeb,
        chargeprod,
        hfct0 = eta,
        hfct1 = pow(hfct0, 3.0);

    double
        r,
        rsq,
        rx,  // dummy variables for distances
        ry,
        rz,
        radius = pow(shortRangeCut, 2.0);  // real space madelung and short range cut are same
                                           // for image convention so just square

    realenergy = 0.0;
    twoenergy = 0.0;
    selfEnergy = 0.0;

    if (doEwald == 1)
    {

        //for (int aatom = rank; aatom < numAtoms; aatom = aatom + numProcs)
        for (int aatom = 0; aatom < numAtoms; aatom++)
        {

            hfunc2(0.0, h0);

            chargea = atomCharge[aatom];
            chargeprod = 0.5 * cFact * chargea * chargea;
            selfEnergy += -h0 * chargeprod * hfct0;

        }

    }

    //std::cout << "\n FIELD PBC nbrlist ewald: " << doEwald;

    ltypea = atmLabel[aatom];
    chargea = cFact * atomCharge[aatom];
    idx = nbrList.maxNbrs * aatom;

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];

    for (int nbr = 0; nbr < nbrList.numNbrs[aatom]; nbr++)
    {
    
        batom = nbrList.map[idx+nbr];
        chargeb = atomCharge[batom];
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

        eng = 0.0;
        eng_s = 0.0;
    
        if (rsq <= radius)
        {

            r = sqrt(rsq);

            h0 = 0.0;

            chargeprod = chargea * chargeb;

            if (doEwald == 1)
            {
                hfunc1(hfct0 * r, h0);

                eng = h0 * chargeprod * hfct0;
            }
            else if (doEwald == 2)
            {
                //evaluate damped and shifted
                interpolate_shifted_damped_potential(r, chargeprod, shortRangeCut, eng);

            }

            if (r <= shortRangeCut && r > 1.0e-4)
            {
                vdw.shortRange(ltypea, ltypeb, r, rsq, eng_s);
            }

        }

        tmp_energy_r += eng;
        tmp_energy_s += eng_s;

    } // end of loop over batom

    realenergy = tmp_energy_r;
    twoenergy = tmp_energy_s;

}

void inline Field::interpolate_shifted_damped_potential(double r, double chg, double cutoff, double& e0)
{
    double
        tmp1,
        tmp2,
        tmp3,
        d1,
        d2,
        dlt,
        omega,
        gamma,
        rstep = 1.0 / delMesh;

    int
        point = int(r * rstep);

    dlt = r * rstep - double(point);

    e0 = 0.0;

    //interpolate energy
    tmp1 = erfc_e[point];
    tmp2 = erfc_e[point + 1];
    tmp3 = erfc_e[point + 2];

    d1 = tmp1 + (tmp2 - tmp1) * dlt;
    d2 = tmp2 + (tmp3 - tmp2) * (dlt - 1.0);

    omega = d1 + (d2 - d1) * dlt * 0.5;

    e0 = chg * (omega - rfld0 + rfld1 * (r - cutoff));

}

void inline Field::interpolate_shifted_damped_potential(double r, double chg, double shortRangeCut, double& e0, double& f1)
{
    double
        tmp1,
        tmp2,
        tmp3,
        d1,
        d2,
        dlt,
        omega,
        gamma,
        rstep = 1.0 / delMesh;

    int
        point = int(r * rstep);

    dlt = r * rstep - double(point);

    e0 = f1 = 0.0;

    //interpolate energy
    tmp1 = erfc_e[point];
    tmp2 = erfc_e[point + 1];
    tmp3 = erfc_e[point + 2];

    d1 = tmp1 + (tmp2 - tmp1) * dlt;
    d2 = tmp2 + (tmp3 - tmp2) * (dlt - 1.0);

    omega = d1 + (d2 - d1) * dlt * 0.5;

    e0 = chg * (omega - rfld0 + rfld1 * (r - shortRangeCut));

    //interpolate force
    tmp1 = erfc_f[point];
    tmp2 = erfc_f[point + 1];
    tmp3 = erfc_f[point + 2];

    d1 = tmp1 + (tmp2 - tmp1) * dlt;
    d2 = tmp2 + (tmp3 - tmp2) * (dlt - 1.0);

    gamma = d1 + (d2 - d1) * dlt * 0.5;

    f1 = -chg * (gamma - rfld1 / r);

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
/*************************************************************
 * function for ewald sum
 * **********************************************************/
void inline Field::hfunc1(double x, double& h0, double& h1)
{
    double
        fact = 1.128379167095512,
        expon,
        xsq;

    if (x == 0.0)
        return;

    h0 = erfc(x);
    h0 = h0 / x;

    xsq = x * x;
    expon = fact * exp(-xsq);
    h1 = (-h0 - expon) / xsq;

    return;
}

void inline Field::hfunc2(double x, double& h0, double& h1)
{
    double
        fact = 1.128379167095512,
        expon,
        xsq;
    xsq = x * x;

    if (xsq < 0.01)
    {
        h0 = fact * ((0.1 * xsq - 0.33333333333) * xsq + 1.0);
        h1 = fact * ((-xsq * 0.1428571429 + 0.4) * xsq - 0.6666666666667);
    }
    else
    {
        h0 = erf(x);
        h0 = h0 / x;

        expon = fact * exp(-xsq);
        h1 = (-h0 + expon) / xsq;

    }

    return;
}

void inline Field::hfunc1(double x, double& h0)
{

    if (x == 0.0)
        return;

    h0 = erfc(x);
    h0 = h0 / x;

}

void inline Field::hfunc2(double x, double& h0)

{
    double
        fact = 1.128379167095512,
        xsq = x * x;
       
    
    if (xsq < 0.01)
    {
        h0 = fact * ((0.1 * xsq - 0.33333333333) * xsq + 1.0);

    }
    else
    {
        h0 = erf(x);
        h0 = h0 / x;

    }

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
    if (doEwald == 1)
        coul.updatePhaseSum();
}

void Field::resetAtomMove(int atm)
{
    //update th neighborlist if required
}
