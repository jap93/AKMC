#include "SaddleSearch.h"

double getRandomNumber(void);

SaddleSearch::SaddleSearch()
{

}

SaddleSearch::~SaddleSearch()
{
        
}

void SaddleSearch::constructActiveRegion(const MPIComms& mpi, Basis& bas, const Species& spec, const JobControl& job, double* centre, std::ofstream& outStream)
{
    int
        rank = mpi.getWorkGroupRank(),
        i = 0, j = 0;

    int
        *coordNum = nullptr;

    vectorSize = 0;

    //first unfreeze all atoms and then reset atoms frozen for all routines
    bas.unFreezeAtoms();
    bas.freezeAtomTypes(spec, job.frozenTypes);

    if (job.kmcParameters.maskHigh.size() > 0)
    {
        for (unsigned i = 0; i < job.kmcParameters.maskHigh.size(); i++)
        {
            outStream << "\n masking region between " << job.kmcParameters.maskLow[i] << " and " << job.kmcParameters.maskHigh[i];
            bas.freezeRegion(job.kmcParameters.maskLow[i], job.kmcParameters.maskHigh[i]);
        }
    }

    //calculate the size of the arrays
    vectorSize = (bas.numberOfAtoms - bas.numFrozenAtoms) * 3;
    if (job.kmcParameters.regionStyle == 1)  // all atoms are displaced except those frozen by default
    {
        vectorSize = (bas.numberOfAtoms - bas.numFrozenAtoms) * 3;   //get the number of atoms being "moved"

        if (job.kmcParameters.useGauss) // if a Gaussian is used a centre is required
        {
        
            if (rank == 0)
            {
                bool found = false;
                while (!found)  // select an atom at random
                {
                    j = getRandomNumber() * bas.numberOfAtoms;
                    if (bas.frozen[j] == 0)  //check that it is allowed to move
                        found = true;
                }
            }

            j = mpi.sumIntGroup(j);

            centre[0] = bas.posX[j];
            centre[1] = bas.posY[j];
            centre[2] = bas.posZ[j];

            if (rank == 0)
                outStream << "\n Gaussian wieghting is used the centre " << centre[0] << " " << centre[1] << " " << centre[2] << std::endl;
        }
        else
        {
            centre[0] = centre[1] = centre[2] = 0.0;
        }

        outStream << "\n no of atoms, number frozen " << bas.numberOfAtoms << " " << bas.numFrozenAtoms << std::endl;
        
    }
    else if (job.kmcParameters.regionStyle == 2)  // all atoms selected at random
    {
        bool found = false;
        double ax, ay, az, rx, ry, rz, xx, yy, zz, rsq;
        double radius = job.kmcParameters.regionRadius * job.kmcParameters.regionRadius;

        j = 0;
        if (rank == 0)
        {
            while (!found)  // select an atom at random
            {
                j = getRandomNumber() * bas.numberOfAtoms;
                if (bas.frozen[j] == 0)  //check that it is allowed to move
                    found = true;
            }
        }

        j = mpi.sumIntGroup(j);

        //work out if atoms are within a given radius allowing for pbc
        ax = bas.posX[j];
        ay = bas.posY[j];
        az = bas.posZ[j];

        centre[0] = ax;
        centre[1] = ay;
        centre[2] = az;

        outStream << "\n atom selected at random with centre " << ax << " " << ay << " " << az << std::endl;

        for (i = 0; i < bas.numberOfAtoms; i++)
        {
            if (i == j || bas.frozen[i] == 1)
                continue;

            rx = bas.posX[i] - ax;
            ry = bas.posY[i] - ay;
            rz = bas.posZ[i] - az;

            xx = rx * bas.rcpVector[0] + ry * bas.rcpVector[3] + rz * bas.rcpVector[6];
            yy = rx * bas.rcpVector[1] + ry * bas.rcpVector[4] + rz * bas.rcpVector[7];
            zz = rx * bas.rcpVector[2] + ry * bas.rcpVector[5] + rz * bas.rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * bas.latVector[0] + yy * bas.latVector[3] + zz * bas.latVector[6];
            ry = xx * bas.latVector[1] + yy * bas.latVector[4] + zz * bas.latVector[7];
            rz = xx * bas.latVector[2] + yy * bas.latVector[5] + zz * bas.latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;
            if (rsq > radius)
            {
                //outStream << "\n distance " << i << " " << ax << " " << ay << " " << az << " " << bas.posX[i] << " " << bas.posY[i] << " " << bas.posZ[i] << " " << rsq;
                bas.frozen[i] = 1; //freeze the atom and increase counter
                bas.numFrozenAtoms++;
            }
        }

        outStream << "\n no of atoms, number frozen " << bas.numberOfAtoms << " " << bas.numFrozenAtoms << std::endl;
        
        vectorSize = (bas.numberOfAtoms - bas.numFrozenAtoms) * 3;   // finally get the number of atoms being "moved"
    }
    else if (job.kmcParameters.regionStyle == 4) // local types
    {
        outStream << "\n local type being found " << std::endl;
        outStream.flush();
        int num = spec.getNumSpecies();
        int typ = -1;
        std::string atomName;
        bool found = false;
        double ax, ay, az, rx, ry, rz, xx, yy, zz, rsq;
        double radius = job.kmcParameters.regionRadius * job.kmcParameters.regionRadius;
        Element ele;

        //get the integer type label of the atom species
        i = int(job.kmcParameters.artTypes1.size() * getRandomNumber());

        atomName = job.kmcParameters.artTypes1[i];
        
        for(j = 0; j < num; j++)
        {
            ele = spec.getSpecies(j);
            
            if (ele.name == atomName)
            {
                typ = j;
                break;
            }
        }
        
        //check the type has been found
        if (typ < 0)
        {
            outStream << "\n the atom type has not been found in the species list" << std::endl;
            outStream.flush();
            std::cerr << "\n the atom type has not been found in the species list" << std::endl;
            std::cerr.flush();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }

        j = 0;
        if (rank == 0)
        {
            while (!found)  // select an atom at random
            {
                j = getRandomNumber() * bas.numberOfAtoms;
                if (bas.atmLabel[j] == typ && bas.frozen[j] == 0)  //check that it is allowed to move
                    found = true;
            }
        }

        j = mpi.sumIntGroup(j);

        //work out if atoms are within a given radius allowing for pbc - it does not depend type
        ax = bas.posX[j];
        ay = bas.posY[j];
        az = bas.posZ[j];

        outStream << "\n atom " << ele.name << " " << j << " selected with centre " << ax << " " << ay << " " << az << std::endl;

        for (i = 0; i < bas.numberOfAtoms; i++)
        {
            if (i == j || bas.frozen[i] == 1)
                continue;

            rx = bas.posX[i] - ax;
            ry = bas.posY[i] - ay;
            rz = bas.posZ[i] - az;

            xx = rx * bas.rcpVector[0] + ry * bas.rcpVector[3] + rz * bas.rcpVector[6];
            yy = rx * bas.rcpVector[1] + ry * bas.rcpVector[4] + rz * bas.rcpVector[7];
            zz = rx * bas.rcpVector[2] + ry * bas.rcpVector[5] + rz * bas.rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * bas.latVector[0] + yy * bas.latVector[3] + zz * bas.latVector[6];
            ry = xx * bas.latVector[2] + yy * bas.latVector[4] + zz * bas.latVector[7];
            rz = xx * bas.latVector[3] + yy * bas.latVector[5] + zz * bas.latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;
            if (rsq > radius)
            {
                //outStream << "\n distance " << i << " " << ax << " " << ay << " " << az << " " << bas.posX[i] << " " << bas.posY[i] << " " << bas.posZ[i] << " " << rsq;
                bas.frozen[i] = 1; //freeze the atom and increase counter
                bas.numFrozenAtoms++;
            }
        }

        outStream << "\n no of atoms, number frozen " << bas.numberOfAtoms << " " << bas.numFrozenAtoms << std::endl;
        
        vectorSize = (bas.numberOfAtoms - bas.numFrozenAtoms) * 3;   // finally get the number of atoms being "moved"
    }
    else if (job.kmcParameters.regionStyle == 5)  // check whether coordination is different from that expected
    {
        int num = spec.getNumSpecies();
        int choice = -1;
        int typ1 = -1;
        int typ2 = -1;
        std::string atomName;
        bool found = false;
        double ax, ay, az, rx, ry, rz, xx, yy, zz, rsq;
        double radius = job.kmcParameters.regionRadius * job.kmcParameters.regionRadius;
        Element ele;

        //get the integer type label of the first atom species
        choice = int(job.kmcParameters.artTypes1.size() * getRandomNumber());

        atomName = job.kmcParameters.artTypes1[choice];

        for(j = 0; j < num; j++)
        {
            ele = spec.getSpecies(j);
            if (ele.name == atomName)
                typ1 = j;
        }

        //check the type has been found
        if (typ1 < 0)
        {
            outStream << "\n the atom type has not been found in the species list" << std::endl;
            outStream.flush();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }

        //get the integer type label of the second atom species
        atomName = job.kmcParameters.artTypes1[choice];

        for(j = 0; j < num; j++)
        {
            ele = spec.getSpecies(j);
            if (ele.name == atomName)
                typ2 = j;
        }

        //check the type has been found
        if (typ2 < 0)
        {
            outStream << "\n the atom type has not been found in the species list" << std::endl;
            outStream.flush();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }

        //calculate the coordination nos of 
        coordNum = bas.calculateCoordNumber(typ1, typ2, job.kmcParameters.coordDistance[choice]);

        j = 0;
        if (rank == 0)
        {
            while (!found)  // select an atom at random
            {
                j = getRandomNumber() * bas.numberOfAtoms;
                if (bas.atmLabel[j] == typ1 && bas.frozen[j] == 0 && (coordNum[j] - job.kmcParameters.idealCoordNumber[choice]) != 0)  //check that it is allowed to move
                    found = true;
            }
        }

        j = mpi.sumIntGroup(j);

        //work out if atoms are within a given radius allowing for pbc - it does not depend type
        ax = bas.posX[j];
        ay = bas.posY[j];
        az = bas.posZ[j];

        outStream << "\n central atom chosen " << j << " " << ax << " " << ay << " " << az << std::endl;;

        for (i = 0; i < bas.numberOfAtoms; i++)
        {
            //outStream << "\n " << i << " " << bas.atmLabel[i] << " " << coordNum[i] << " " << coordNum[i] - job.idealCoordNumber[choice];
            if (i == j || bas.frozen[i] == 1)
                continue;

            rx = bas.posX[i] - ax;
            ry = bas.posY[i] - ay;
            rz = bas.posZ[i] - az;

            xx = rx * bas.rcpVector[0] + ry * bas.rcpVector[3] + rz * bas.rcpVector[6];
            yy = rx * bas.rcpVector[1] + ry * bas.rcpVector[4] + rz * bas.rcpVector[7];
            zz = rx * bas.rcpVector[2] + ry * bas.rcpVector[5] + rz * bas.rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * bas.latVector[0] + yy * bas.latVector[3] + zz * bas.latVector[6];
            ry = xx * bas.latVector[1] + yy * bas.latVector[4] + zz * bas.latVector[7];
            rz = xx * bas.latVector[2] + yy * bas.latVector[5] + zz * bas.latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;
            if (rsq > radius)
            {
                bas.frozen[i] = 1; //freeze the atom and increase counter
                bas.numFrozenAtoms++;
            }
            //else
            //{
            //    outStream << "\n distance " << i << " " << ax << " " << ay << " " << az << " " << bas.posX[i] << " " << bas.posY[i] << " " << bas.posZ[i] 
            //              << " " << sqrt(rsq);
            //}
            
        }

        free_ivector(coordNum);

        outStream << "\n no of atoms, number frozen " << bas.numberOfAtoms << " " << bas.numFrozenAtoms << std::endl;
        
        vectorSize = (bas.numberOfAtoms - bas.numFrozenAtoms) * 3;   // finally get the number of atoms being "moved"
    }
    
    outStream << "\n finished constructing atoms " << rank;
    outStream.flush();

}

void SaddleSearch::makeRandomDisplacement(const MPIComms& mpi, Basis& bas, double initialDisplacement, std::ofstream& outStream)
{
    int
        rank = mpi.getWorkGroupRank(),
        grpSize = mpi.getNumWorkGroupProcs();

    double
        *deltaPos = nullptr;

    bool
        problem = false;

    if (vectorSize <= 0)
    {
        outStream << "\n vectorSize is zero - something has gone wrong with atom selection" << std::endl;
        outStream.flush();
        std::cerr << "\n vectorSize is zero - something has gone wrong with atom selection on rank " << rank << std::endl;
        std::cerr.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    //do a sanit check that vectorsize is the same on all cores within a group - it may indicate a problem
    //start the random displacement
    rVector = alloc_dvector(vectorSize, "ART::rVector", 0.0);
      
    deltaPos = alloc_dvector(vectorSize, "initial displacement of atoms", 0.0);

    if (rank == 0)
    {
        int j = 0;
        for (int i = 0; i < bas.numberOfAtoms; i++)
        {
            if (bas.frozen[i] == 0)
            {
                deltaPos[j] = getRandomNumber() - 0.5;
                deltaPos[j+1] = getRandomNumber() - 0.5;
                deltaPos[j+2] = getRandomNumber() - 0.5;
                j += 3;
            }
        }

    }

    //make sure random number id the same on all threads in the group
    mpi.sumDoubleVectorGroup(deltaPos, vectorSize);


    centreOfMass(deltaPos, vectorSize);  //prevent drift
        
    //the displacement has to be normalised
    double vNorm = vectorNorm(deltaPos, vectorSize);
    double invNorm = (1.0 / vNorm);

    for (int i = 0; i < vectorSize; i++)
    {
        rVector[i] = invNorm * deltaPos[i]; // unit vector to move out of basin
        deltaPos[i] *= initialDisplacement;
    }

    //displace the atom positions
    bas.displaceAtoms(deltaPos);

    free_dvector(deltaPos);

}

void SaddleSearch::makeGaussDisplacement(const MPIComms& mpi, Basis& bas, double* centre, double initialDisplacement, double gaussWidth, std::ofstream& outStream)
{
    int
        rank = mpi.getRank(),
        i, j;

    double
        dist = 0.0,
        *deltaPos = nullptr,
        *gaussWeight = nullptr;

    if (vectorSize <= 0)
    {
        outStream << "\n vectorSize is zero - something has gone wrong with atom selection" << std::endl;
        outStream.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    //start the random displacement
    rVector = alloc_dvector(vectorSize, "ART::rVector", 0.0);
        
    deltaPos = alloc_dvector(vectorSize, "initial displacement of atoms", 0.0);
    gaussWeight = alloc_dvector(vectorSize, "initial displacement of atoms", 0.0);
    

    if (rank == 0)
    {
        j = 0;
        for (i = 0; i < bas.numberOfAtoms; i++)
        {
            if (bas.frozen[i] == 0)
            {
                deltaPos[j] = getRandomNumber() - 0.5;
                deltaPos[j+1] = getRandomNumber() - 0.5;
                deltaPos[j+2] = getRandomNumber() - 0.5;
    
                dist = bas.getDistance(centre, i);
                dist = dist * dist;
                gaussWeight[j] = exp(-dist / (2.0 * pow(gaussWidth,2.0)));
                gaussWeight[j+1] = exp(-dist / (2.0 * pow(gaussWidth,2.0)));
                gaussWeight[j+2] = exp(-dist / (2.0 * pow(gaussWidth,2.0)));
                j += 3;
            }
        }
    }

    mpi.sumDoubleVectorGroup(deltaPos, vectorSize);
    mpi.sumDoubleVectorGroup(gaussWeight, vectorSize);
    
    centreOfMass(deltaPos, vectorSize);  //prevent drift
        
    //the displacement has to be normalised
    double vNorm = vectorNorm(deltaPos, vectorSize);
    double invNorm = (1.0 / vNorm);

    for (i = 0; i < vectorSize; i++)
    {
        rVector[i] = invNorm * deltaPos[i]; // unit vector to move out of basin
        deltaPos[i] *= initialDisplacement * gaussWeight[i];
    }

    //displace the atom positions
    bas.displaceAtoms(deltaPos);

    free_dvector(deltaPos);
    free_dvector(gaussWeight);

}
