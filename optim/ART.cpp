#include "ART.h"

double getRandomNumber(void);
std::vector<std::string> split(std::string s);

extern "C" int dstev_(char *jobz, int *n, double *d__, double *e, double *z__, int *ldz, double *work, int *info);

ART::ART()
{

}

ART::~ART()
{
        
}
void ART::findSaddlePoint(const MPIComms& mpi, const Species& spec, const JobControl& job, Field* fld, Basis& basin1, Basis& saddle,
                          Status& artStatus, Energy& basin1Energy, Energy& saddleEnergy, double* centre, bool continueSearch,
                          std::ofstream& outStream)
{
    int
        rank = mpi.getRank();

    

    // do some error checking first
    numLanczosVectors = job.artParameters.numLanczosVectors;

    artStatus.setStatusFailed();


    //make a copy of the initial basin for the saddle
    saddle = basin1;
    bool restart = false;
    if (restart)
    {  // this is where restart is read in
       // at this point read in the saddle
        double restartEnergy = 0.0;
        std::ifstream resStream;
        resStream.open("saddle.res", std::ios::in);
        saddle.inputBasis(itART, restartEnergy, spec, job, resStream, outStream);

        double tstep;
        loadARTDetails(saddle.frozen, tstep, saddle.numberOfAtoms, resStream);
        resStream.close();

        saddleEnergy.setTotalEnergy(restartEnergy);

        if (rank == 0)
        {
            outStream << "\n ART restarting on iteration " << itART << std::endl;
            outStream << " the energy of the saddle restart is " << saddleEnergy.getTotalEnergy() << " " << saddleEnergy.getEnergyUnit() << std::endl;
            outStream << " centre " << centre[0] << " " << centre[1] << " " << centre[2] << std::endl;
            outStream << " vectorSize " << vectorSize << std::endl;
            for (int i = 0; i < vectorSize; i++)
            {
                outStream << "\n rVector " << i << " " << rVector[i];
                outStream.flush();
            }
        }

	    centre[0] = ARTCentre[0];
	    centre[1] = ARTCentre[1];
	    centre[2] = ARTCentre[2];

    }
    else
    {
        //make a frozen list in the atom configuration dependent on the requirements
        constructActiveRegion(mpi, saddle, spec, job, centre, outStream);

        //randomly displaced unfrozen ions
        if (job.kmcParameters.useGauss)
        {
            makeGaussDisplacement(mpi, saddle, centre, job.artParameters.initialDisplacement, job.kmcParameters.gaussWidth, outStream);
        }
        else
        {
            makeRandomDisplacement(mpi, saddle, job.artParameters.initialDisplacement, outStream);
        }

        checkVectorSize();
        itART = 0;

	    ARTCentre[0] = centre[0];
	    ARTCentre[1] = centre[1];
	    ARTCentre[2] = centre[2];
    }

    if (job.artParameters.goPeprepedicular)
    {

        FIRESaddlePoint(mpi, saddle, spec, rVector, saddleEnergy, fld, job.artParameters, artStatus, true, outStream);

        //check that a transition point was found and that there was no interrupt
        if (artStatus.getStatus() == INTERRUPTED)
        {
            outStream << "\n interrupted about to return";
            outStream.flush();
        }
        
        if (artStatus.getStatus() != SUCCESS)
        {
            free_dvector(rVector);
            return;
        } 

    }
                
    FIRESaddlePoint(mpi, saddle, spec, rVector, saddleEnergy, fld, job.artParameters, artStatus, false, outStream);
     

    if (artStatus.getStatus() == SUCCESS)
    {
        basin1.reportBasisDifference(mpi, saddle, spec, job.kmcParameters.kmcBasinRadius, outStream);
    
        if (rank == 0)
        {
            outStream << "\n basin difference at saddle point : " << std::endl;
            outStream.flush();
        }
    }

    //free up the memory
    free_dvector(rVector);

   //mpi.commsAbortWorld();
   //exit(EXIT_FAILURE)  ;   
}



bool ART::calculateLanczosForce(const MPIComms& mpi, const ARTCntrl& artParameters, Basis& bas, double& lowestEigenValue, double* forcePar, 
                                double* forcePerp, Field* fld, Energy& totalEnergy, std::ofstream& outStream)
{
    int
        i,
        j,
        idx = 0,
        idxm1 = 0;

    double
        alpha,
        beta,
        delta,
        eigValNew,
        eigValOld;    // beta is the norm of the "r" vector

    double
        *qx = nullptr,
        *qVector = nullptr,
        *uVector = nullptr,
        *initForce = nullptr,
        *eigVecNew = nullptr,
        *deltaForce = nullptr,
        *Q = nullptr, // ttridiagonal matrix Q
        *diagonalT = nullptr, // tridiagonal matrix T 
        *subDiagT = nullptr;
        
    bool
        converged = false;  //converged covers both failure to converge and convergence to positive eigenvalue

    Energy
        dummyEnergy;

    lowestEigenValue = 1e6; // set to high positive number


    qx = alloc_dvector(vectorSize, " ART method ");
    qVector = alloc_dvector(vectorSize, " ART method ");
    uVector = alloc_dvector(vectorSize, " ART method ");
    initForce = alloc_dvector(vectorSize, " ART method ");
    deltaForce = alloc_dvector(vectorSize, " ART method ");
    eigVecNew = alloc_dvector(numLanczosVectors, " ART method ");
    Q = alloc_dvector(numLanczosVectors * vectorSize, "ART tridiagonal matrix Q", 0.0);
    diagonalT = alloc_dvector(numLanczosVectors, "ART tridiagonal matrix T", 0.0);
    subDiagT = alloc_dvector(numLanczosVectors, "ART tridiagonal matrix T", 0.0);
    
    //calculate beta for iteration i = 0
    beta = vectorNorm(rVector, vectorSize);
    
    //calculate the forces for the initial positions
    calculateForces(mpi, initForce, bas, fld, totalEnergy, outStream);
    //outStream << "\n init force " << std::endl;
    //totalEnergy.printEnergy(0, outStream);


    i = 0;
    do   //do equations 2 - 7
    {

        idx = vectorSize * i;
        idxm1 = vectorSize * (i - 1);

        for (j = 0; j < vectorSize; j++)
        {
            qVector[j] = rVector[j] / beta;      //eqn 2
            Q[idx+j] = qVector[j]; //Q(i,j) = qVector(j)
        }
        
        //do the finite difference part - eqns 3, 9 and 10
        //first create the vector q.r*
        for (j = 0; j < vectorSize; j++)
        {
            qx[j] = qVector[j] * artParameters.deltaX;
            //outStream << "\n qs " << j << " " << qx[j] << " " << qVector[j] << " " << artParameters.deltaX;
        }

        bas.displaceAtoms(qx);

        //outStream << " force vector iteration " << i;
        calculateForces(mpi, deltaForce, bas, fld, dummyEnergy, outStream);
        //outStream << "\n displaced  force " << std::endl;
        //totalEnergy.printEnergy(0, outStream);

        //put atoms back
        for (j = 0; j < vectorSize; j++)
             qx[j] *= -1.0;
        bas.displaceAtoms(qx);

        //numerical derivative
        for (j = 0; j < vectorSize; j++)
            uVector[j] = -(deltaForce[j] - initForce[j]) / artParameters.deltaX;


        if (i == 0) // equation 4
        {
            for (j = 0; j < vectorSize; j++)
                rVector[j] = uVector[j];
        }
        else
        {
            for (j = 0; j < vectorSize; j++)
            {
                rVector[j] = uVector[j] - beta * Q[idxm1+j]; //rVector(j) = uVector(j) - beta * Q(i-1,j); 
            }
        }

        alpha = dotProduct(qVector, rVector, vectorSize);  //equation 5
        
        for (j = 0; j < vectorSize; j++)
        {
            rVector[j] = rVector[j] - alpha * Q[idx+j];     //equation 6 i.e. rVector(j) = rVector(j) - alpha * Q(i,j); 
        }
            

        //diagonal part of matrix T
        diagonalT[i] = alpha;

        //sub diagonal parts of matrix T
        //ah blast must be greater than 0!
        if ( i > 0) 
            subDiagT[i-1] =  beta;
            

        //update beta for the next iteration - eqn 7
        beta = vectorNorm(rVector, vectorSize);
    
        //calculate eigenvalue and check Eigenvalues have not converged
        if (i == 0) 
        {
            eigValNew = eigValOld = alpha;
        }
        else
        {
            
            try
            {
                diagonalise(diagonalT, subDiagT, eigValNew, eigVecNew, i + 1); //get the lowest eigenvalue and corresponding eigenvectors
            }
            catch(const std::invalid_argument& e)
            {
                outStream << e.what() << std::endl;

                //delete memory
                free_dvector(qx);
                free_dvector(qVector);
                free_dvector(uVector);
                free_dvector(deltaForce);
                free_dvector(initForce);
                free_dvector(eigVecNew);

                free_dvector(Q);
                free_dvector(diagonalT);
                free_dvector(subDiagT);
                
                //rethrow exception
                throw;
            }
            
            
            //get the lowest eigen value tolerance - eqn 8
            delta = abs((eigValNew - eigValOld) / eigValOld);

            /*if (mpi.getWorkGroupRank() == 0)
            {
                outStream << "\n eigenval for vector " << i << " " << eigValNew << " " << eigValOld << " " << delta;
                outStream.flush();
            }*/
            //set to the old value for the next iteration
            eigValOld = eigValNew;
            lowestEigenValue = eigValNew;
            if (delta < artParameters.eigValTol) 
            {
                
                converged = true;
                //iterConverged = i;

            } 

        }
        
        i++;

    } while (i < numLanczosVectors && converged == false);


    // v = Qv^t i.e. calculate the eigenvector of full Hessian the bit on page 9778 - bottom left para
    double* eigenVector = alloc_dvector(vectorSize, " lanczos eigenvector ", 0.0);
    double* fvv = alloc_dvector(vectorSize, " lanczos eigenvector ", 0.0);

    for (i = 0; i < numLanczosVectors; i++)
    {
        idx = i * vectorSize;
        for (j = 0; j < vectorSize; j++)
        {
            eigenVector[j] += eigVecNew[i] * Q[idx+j];   ////eigenVector(j) += eigVecNew[i] * triDiagQ[i][j];  //numLanczosVectors, vectorSize
        }
    }

    double norm = vectorNorm(eigenVector, vectorSize);
    for (j = 0; j < vectorSize; j++)
        eigenVector[j] = eigenVector[j] / norm;

    //calculate forces using equation 18
    double dot;
    for (i = 0; i < vectorSize; i++)
    {
        dot = dotProduct(initForce, eigenVector, vectorSize);
        for (j = 0; j < vectorSize; j++)
            fvv[j] = -dot * eigenVector[j];
    }
    

    for (i = 0; i < vectorSize; i++)
        forcePar[i] = initForce[i] + 2.0 * fvv[i];

    for (i = 0; i < vectorSize; i++)
            forcePerp[i] = fvv[i];
  

    for (i = 0; i < vectorSize; i++)
        rVector[i] = eigenVector[i];

    
    //destroy rest of temporary arrays
    free_dvector(eigenVector);
    free_dvector(fvv);
    free_dvector(qx);
    free_dvector(qVector);
    free_dvector(uVector);
    free_dvector(deltaForce);
    free_dvector(initForce);
    free_dvector(eigVecNew);

    free_dvector(Q);
    free_dvector(diagonalT);
    free_dvector(subDiagT);

    return converged;
}
void ART::calculateForces(const MPIComms& mpi, double* forceVector, Basis& bas, Field* fld, Energy& totalEnergy, std::ofstream& outStream)
{
    

    double
        *forceX, 
        *forceY, 
        *forceZ,
        stress[6] = {0.0};

    Energy
        eng;

    forceX = alloc_dvector(bas.numberOfAtoms, "MD: force", 0.0);
    forceY = alloc_dvector(bas.numberOfAtoms, "MD: force", 0.0);
    forceZ = alloc_dvector(bas.numberOfAtoms, "MD: force", 0.0);

    fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceX, forceY, forceZ, stress, bas.latVector, bas.rcpVector, bas.charge, 
                         bas.atmLabel, bas.frozen, bas.numberOfAtoms, eng);

    totalEnergy = eng;
    
    //map the forces onto the single vector model
    int j = 0;
    //int grpID = mpi.getWorkGroupRank();
    for (int i = 0; i < bas.numberOfAtoms; i++)
    {
        //outStream << "\n ART " <<  i << " " << forceX[i] << " " << forceY[i] << " " << forceZ[i];
        if (bas.frozen[i] == 0)  //i.e. it is free to move
        {
            forceVector[j] = forceX[i]; 
            forceVector[j+1] = forceY[i]; 
            forceVector[j+2] = forceZ[i]; 
            //if (grpID == 0) std::cout << "\n ART " << grpID << " " << i << " " << forceVector[j] << " " << forceVector[j+1] << " " << forceVector[j+2];
            j += 3;
        }
    }

    //free memeory
    free_dvector(forceX);
    free_dvector(forceY);
    free_dvector(forceZ);

//std::cout.flush();
//exit(-1);
}



void ART::FIRESaddlePoint(const MPIComms& mpi, Basis& bas, const Species& spec, double* rVector, Energy& finalenergy, Field* fld, const ARTCntrl& artParameters, Status& status, 
                          bool goPerpendicular, std::ofstream& outStream)
{
    int
        rank = mpi.getWorkGroupRank(),
        i, j,
        maxIterations = 0,
        nStep = 1;

    double
        eigenValue,
        timeStep,
        timeStepMax,
        halfStep,
        halfStepSquared,
        alpha = 0.01, //alphaStart,
        OneLessAlpha,
        deltaE,
        fNormPar,
        fNormPerp,
        vNorm,
        fMaxPar,
        fMaxPerp,
        vMax;

    double
        P;

    double
        *velocity,
        *forceOld,
        *forcePar,
        *forcePerp;

    bool 
        convergedFire = false,
        convergedLanzos = false;

    Energy
        totalEnergyOld,
        totalEnergyNew;

    std::ofstream resStream;
    
    //Set FIRE parameters
    alpha = artParameters.alphaStart;

    if (goPerpendicular)
    {
        maxIterations = artParameters.maxPerpIter;
        timeStep =  artParameters.initTimeStepPerp;
        timeStepMax = artParameters.maxTimeStepPerp;
    }
    else
    {
        maxIterations = artParameters.maxParIter;
        timeStep =  artParameters.initTimeStep;
        timeStepMax = artParameters.maxTimeStep;
    }

    // allocate work arrays
    velocity = alloc_dvector(vectorSize, "FIRE", 0.0);
    forceOld = alloc_dvector(vectorSize, "FIRE", 0.0);
    forcePar = alloc_dvector(vectorSize, "FIRE", 0.0);
    forcePerp = alloc_dvector(vectorSize, "FIRE", 0.0);


    //calulate new forces
    try
    {
        convergedLanzos = calculateLanczosForce(mpi, artParameters, bas, eigenValue, forcePar, forcePerp, fld, totalEnergyOld, outStream);
    }
    catch(const std::invalid_argument& e)
    {
        status.setStatusFailed();
        free_dvector(velocity);
        free_dvector(forceOld);
        free_dvector(forcePar);
        free_dvector(forcePerp);
        return;
    }
    catch(const std::runtime_error& e)
    {
        status.setStatusFailed();
        free_dvector(velocity);
        free_dvector(forceOld);
        free_dvector(forcePar);
        free_dvector(forcePerp);
        return;
    }

    /* it may may not converge on first pass but appears to rectify itself on the next step
    if (!converged)                     
    {
        status.setStatusFailed();
        free_dvector(velocityHalf);
        free_dvector(velocityFull);
        free_dvector(forceOld);
        free_dvector(forcePar);
        free_dvector(forcePerp);
        return;
    } */

    if (artParameters.debug && rank == 0)
    {
        outStream << "\n FIRE Saddle Search: initial energy " << totalEnergyOld.getTotalEnergy() << " " << totalEnergyOld.getEnergyUnit() << std::endl;
        outStream.flush();
    }

    halfStep = 0.5 * timeStep * totalEnergyOld.getTimeStepUnit();
    halfStepSquared = 0.5 * timeStep * timeStep;

    while (!convergedFire && itART <= maxIterations)
    {
        itART++;  // increase counter
	/*if (rank == 0)
	{
            resStream.open("saddle.res", std::ofstream::out);
            bas.title = " ART restart file";
            bas.dumpBasis(mpi, spec, 0.0, totalEnergyOld.getTotalEnergy(), itART, resStream);
            //saveARTDetails(mpi, bas.frozen, timeStep, bas.numberOfAtoms, resStream);
            resStream.flush();
            resStream.close();
	}*/

        //the MD part of FIRE - note mass is set to 1.0 for FIRE

        if (eigenValue >= artParameters.maxEigenvalue  || goPerpendicular == true)   //goPerpendicular == true)
        {
            for (i = 0; i < vectorSize; i++)
            {
                forceOld[i] = forcePerp[i];
            }
        }
        else
        {
            for (i = 0; i < vectorSize; i++)
            {
                forceOld[i] = forcePar[i] * artParameters.parDamp;
            }
        }
        
        //outStream << "before " << bas.rpos[i].x << " " << bas.rpos[i].y  << " " << bas.rpos[i].z << std::endl;
        for (i = 0; i < vectorSize; i++)
        {
            velocity[i] += halfStep * forceOld[i];   
        } 

        //centreOfMass(velocityHalf, vectorSize);
        
        
        //!!FIRE Update
        vectorNormMax(forcePar, fNormPar, fMaxPar, vectorSize);
        vectorNormMax(forcePerp, fNormPerp, fMaxPerp, vectorSize);
        vectorNormMax(velocity, vNorm, vMax, vectorSize);



        // F1 of program
        if (eigenValue >= artParameters.maxEigenvalue || goPerpendicular == true)
        {
            P = dotProduct(forcePerp, velocity, vectorSize);

            OneLessAlpha = 1.0 - alpha;
            // this is part F2 of the paper
            for (i = 0; i < vectorSize; i++)
            {
                velocity[i] = OneLessAlpha * velocity[i] + (alpha * forcePerp[i] * vNorm / fNormPerp);
            }
        }
        else
        {
            P = dotProduct(forcePar, velocity, vectorSize);

            OneLessAlpha = 1.0 - alpha;
            // this is part F2 of the paper
            for (i = 0; i < vectorSize; i++)
                velocity[i] = OneLessAlpha * velocity[i] + (alpha * forcePar[i] * vNorm / fNormPar);
        }
        
        if(P > 0)  // downhill move
        {
            if (nStep > artParameters.N_min)
            {   //F3

                timeStep = fmin(timeStep * artParameters.stepIncrease, timeStepMax); // increase the time step
                alpha *= artParameters.alphaDecrease;                  // decrease alpha
            }

            nStep += 1;

        }
        else // uphill move (P <= 0)
        {   //F4

            nStep = 0;                        // reset counter for the number of steps between N_min
            timeStep *= artParameters.stepDecrease;         // decrease timestep
            alpha = artParameters.alphaStart;               // reset alpha to start

            for (i = 0; i < vectorSize; i++)
                velocity[i] = 0.0;              // freeze system

        }

        if (artParameters.debug)
            outStream << "\n ************************* iteration " << itART << " *************************";
        //3: x (t + t ) <- x (t ) + t ·v (t + 1/2 t )
        j = 0;
        double dx, dy, dz;
        double maxDisplacement = -1.0e10;
        double delta = 0.0;
        int atmMax = -1;
        for (i = 0; i < bas.numberOfAtoms; i++)
        {
            // advance coordinates
            if(bas.frozen[i] == 0)
            {
                dx = timeStep * velocity[j];
                dy = timeStep * velocity[j+1];
                dz = timeStep * velocity[j+2];
                bas.posX[i] += dx;
                bas.posY[i] += dy;
                bas.posZ[i] += dz;

                j += 3;

                delta = sqrt(dx * dx + dy * dy + dz * dz);
                if (delta > maxDisplacement)
                {
                    maxDisplacement = delta;
                    atmMax = i;
                }
                //outStream << "\n" << i << " " << bas.posX[i]  << " " << bas.posY[i]  << " " << bas.posZ[i];
            }
            
        }

        convergedLanzos = calculateLanczosForce(mpi, artParameters, bas, eigenValue, forcePar, forcePerp, fld, totalEnergyNew, outStream);

        if (eigenValue >= artParameters.maxEigenvalue || goPerpendicular == true)
        {
            for (i = 0; i < vectorSize; i++)
                velocity[i] += halfStep * forcePerp[i];
        }
        else
        {
            for (i = 0; i < vectorSize; i++)
                velocity[i] += halfStep * forcePar[i];
        }

        //centreOfMass(velocityFull, vectorSize);

        if (artParameters.debug && itART % artParameters.print == 0)
        {
            outStream << "\n FIRE: total energy " << " " << totalEnergyNew.getTotalEnergy() << " eigenvalue " << eigenValue << " " 
                      << totalEnergyNew.getEnergyUnit() << " timestep " << timeStep << std::endl;

            outStream << " FIRE: gnorm perp     " << fNormPerp << " force max " << fMaxPerp << " energy difference " 
                          << fabs(totalEnergyNew.getTotalEnergy() - totalEnergyOld.getTotalEnergy()) << std::endl;
            outStream << " FIRE: gnorm parallel " << fNormPar << " force max " << fMaxPar << " energy difference " 
                          << fabs(totalEnergyNew.getTotalEnergy() - totalEnergyOld.getTotalEnergy()) << std::endl;
            outStream << " FIRE: the maximum displacemt for atom " << atmMax << " is " << maxDisplacement << std::endl;
                
        }
        outStream.flush();

        // it is recomended that convergence is tested here
        deltaE = fabs(totalEnergyOld.getTotalEnergy() - totalEnergyNew.getTotalEnergy());

        if (itART > 5) //only allow if a reasonable number of iterations complete
        {
            
            if ((fNormPar + fNormPerp) < artParameters.normTol || (fMaxPar + fMaxPerp) < artParameters.forceTol) // || deltaE < artParameters.energyTol)
            {

                if (artParameters.debug)
                {
                    outStream << "\n FIRE converged : delta E = " << deltaE << std::endl;
                }

                convergedFire = true;
                if (eigenValue <= artParameters.maxEigenvalue)
                    status.setStatusSuccess();
                else
                    status.setStatusFailed(); // converged with eigenvalue above minimum

                finalenergy = totalEnergyNew;
                free_dvector(velocity);
                free_dvector(forceOld);
                free_dvector(forcePar);
                free_dvector(forcePerp);
                return;
            }
        }

        //the MD part of FIRE - note mass is set to 1.0 for FIRE
        halfStep = 0.5 * timeStep * totalEnergyOld.getTimeStepUnit();
        halfStepSquared = 0.5 * timeStep * timeStep;
  
        //can now update energies
        totalEnergyOld = totalEnergyNew;

        if (abs(halfStep) < 1.0e-10) // sometimes time step goes almost to zero so cut off quickly
        {
            status.setStatusFailed();
            free_dvector(velocity);
            free_dvector(forceOld);
            free_dvector(forcePar);
            free_dvector(forcePerp);
            return; 
        }

        mpi.checkForInterrupt(status, outStream);

        if (status.getStatus() == TERMINATED)
        {
            outStream << "\n interrupted whilst in ART!" << std::endl;
            outStream.flush();
            status.setStatusTerminated();
            free_dvector(velocity);
            free_dvector(forceOld);
            free_dvector(forcePar);
            free_dvector(forcePerp);
            return;
        } 

    }
    
    //to many iterations
    status.setStatusFailed();

    free_dvector(velocity);
    free_dvector(forceOld);
    free_dvector(forcePar);
    free_dvector(forcePerp);
}

void ART::checkVectorSize(void)
{
    if (vectorSize < numLanczosVectors) // the lanczos force routine requires this condition
        numLanczosVectors = vectorSize;
}

void ART::diagonalise(double* diagonalT, double* subDiagT, double& eigVal, double* eigVec, int n)
{
    int info = 0;

    char jobz = 'V';
    double dt[n];
    double ot[n];
    double z[n*n];
    double work[2*n-2];

    for (int i = 0; i < n; i++)
    {
        dt[i] = diagonalT[i];
        ot[i] = subDiagT[i];
    }

    dstev_(&jobz, &n, dt, ot, z, &n, work, &info);

    if (info != 0)
    {
        std::cout << "non-negative value returned from dstev" << std::endl;
        std::cout.flush();
        throw std::invalid_argument("non-negative value returned from dstev");
    }
    
    eigVal = dt[0];
    
    for (int i = 0; i < n; i++)
        eigVec[i] = z[i];

}
    /* Subroutine  int dstev_(char *jobz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info)
{*/
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEV computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric tridiagonal matrix A.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix A, stored in elements 1 to N-1 of E; E(N) need not   
            be set, but is used by the routine.   
            On exit, the contents of E are destroyed.   

    Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)   
            If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal   
            eigenvectors of the matrix A, with the i-th column of Z   
            holding the eigenvector associated with D(i).   
            If JOBZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= max(1,N).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))   
            If JOBZ = 'N', WORK is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of E did not converge to zero.   
}
*/

void ART::loadARTDetails(int* frozen, double& timeStep, int natom, std::ifstream& inStream)
{
    std::string
        line,
        dummy;

    std::vector<std::string> words;

    std::getline(inStream, line);
    words = split(line);

    vectorSize = std::stoi(words[0]);

    rVector = alloc_dvector(vectorSize, "loadARTDetails rVector", 0.0);

    for (int i = 0; i < vectorSize; i++)
    {
        std::getline(inStream, line);
        words = split(line);
        rVector[i] = std::stod(words[0]);
    }

    std::getline(inStream, line);
    words = split(line);
    ARTCentre[0] = std::stod(words[0]);
    ARTCentre[1] = std::stod(words[1]);
    ARTCentre[2] = std::stod(words[2]);

    for (int i = 0; i < natom; i++)
    {
        std::getline(inStream, line);
        words = split(line);
        frozen[i] = std::stoi(words[0]);
    }

    std::getline(inStream, line);
    words = split(line);
    timeStep = std::stod(words[0]);
}


void ART::saveARTDetails(const MPIComms& mpi, int* frozen, double timeStep, int natom, std::ofstream& outStream)
{
    if (mpi.getRank() != 0)
        return;

    outStream << vectorSize << std::endl;

    for (int i = 0; i < vectorSize; i++)
        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                  << std::setprecision(10) << rVector[i] << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
              << std::setprecision(10) << ARTCentre[0] << " " << ARTCentre[1] << " " << ARTCentre[2] << std::endl;

    for (int i = 0; i < natom; i++)
        outStream << frozen[i] << std::endl;

    outStream << timeStep;
}

