#include "Dimer.h"

double getRandomNumber(void);
std::vector<std::string> split(std::string s);
Dimer::Dimer()
{

}

Dimer::~Dimer()
{
        
}

void Dimer::findSaddlePoint(const MPIComms& mpi, const Species& spec, const JobControl& job, Field* fld, Basis& basin1, Basis& saddle,
                          Status& status, Energy& basin1Energy, Energy& saddleEnergy, double* centre, bool continueSearch,
                          std::ofstream& outStream)
{
    int
        grpRank = mpi.getWorkGroupRank();

    double
        *rVector = nullptr,
        deltaPosX[basin1.numberOfAtoms],
        deltaPosY[basin1.numberOfAtoms],
        deltaPosZ[basin1.numberOfAtoms];

    bool
        foundSaddle = false;

    Status
        relStatus,
        dimerStatus;

    Relax
        rel;


    
    //free_dvector(rVector);
    
}

void Dimer::FIRESaddlePoint(const MPIComms& mpi, Basis& bas, double* rVector, double& finalenergy, Field* fld, const DimerCntrl& dimerParameters, 
                            Status& status, bool goPerpendicular, std::ofstream& outStream)
{
    int
        worldRank = mpi.getRank(),
        i, j,
        it = 0,
        nStep = 1;

    double
        totalEnergyOld,
        totalEnergyNew,
        eigenValue,
        rcpEnergy, realEnergy, vdwEnergy,
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
        Norm,
        fMaxPar,
        fMaxPerp,
        vMax;

    double
        P;

    double
        *velocityOld,
        *velocityNew,
        *forceOld,
        *forcePar,
        *forcePerp;

    bool 
        recalcNbrList = false,
        convergedFire = false,
        convergedLanzos = false;
    
    //Set FIRE parameters
    timeStepMax = 2 * PI * sqrt(alpha) * 1.2e-1;
    timeStep = timeStepMax * 0.5;

    // allocate work arrays
    velocityOld = alloc_dvector(vectorSize, "FIRE", 0.0);
    velocityNew = alloc_dvector(vectorSize, "FIRE", 0.0);
    forceOld = alloc_dvector(vectorSize, "FIRE", 0.0);
    forcePar = alloc_dvector(vectorSize, "FIRE", 0.0);
    forcePerp = alloc_dvector(vectorSize, "FIRE", 0.0);

    //useNbrList = false;
    //calulate new forces
    try
    {
        
        //convergeDimer(mpi, bas, rVector, eigenValue, forcePar, forcePerp, totalEnergyOld, shortRangeCut, verletShell,
        //              fld, nbrList, useNbrList, dimerParameters, status, doEwald, outStream);
    }
    catch(const std::invalid_argument& e)
    {
        status.setStatusFailed();
        free_dvector(velocityOld);
        free_dvector(velocityNew);
        free_dvector(forceOld);
        free_dvector(forcePar);
        free_dvector(forcePerp);
        return;
    }
    catch(const std::runtime_error& e)
    {
        status.setStatusFailed();
        free_dvector(velocityOld);
        free_dvector(velocityNew);
        free_dvector(forceOld);
        free_dvector(forcePar);
        free_dvector(forcePerp);
        return;
    }
    
   /* if(!convergedLanzos)  // if the initial Lanczos force fails then bail out
    {
        status.setStatusFailed();
        free_dvector(velocityOld);
        free_dvector(velocityNew);
        free_dvector(forceOld);
        free_dvector(forcePar);
        free_dvector(forcePerp);
        return;
    } */

    if (dimerParameters.debug)
    {
        outStream << "\n FIRE Saddle Search: initial energy " << totalEnergyOld << std::endl;
        outStream.flush();
    }
    while (!convergedFire && it <= dimerParameters.maxParIter)
    {
        it++;  // increase counter

        if (dimerParameters.debug)
            outStream << "\n fire iteration " << it;
        
        //the MD part of FIRE - note mass is set to 1.0 for FIRE
        halfStep = 0.5 * timeStep;
        halfStepSquared = 0.5 * timeStep * timeStep;

        if (eigenValue >= dimerParameters.maxEigenvalue || goPerpendicular == true)
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
                forceOld[i] = forcePar[i];
            }
        }
        
        j = 0;
        //outStream << "before " << bas.rpos[i].x << " " << bas.rpos[i].y  << " " << bas.rpos[i].z << std::endl;
        for (i = 0; i < bas.numberOfAtoms; i++)
        {
                    
            // advance coordinates and make sure it is only the non-frozen atoms
            if(bas.frozen[i] == 0)
            {
                bas.posX[i] += timeStep * velocityOld[j] + halfStepSquared * forceOld[j];
                bas.posY[i] += timeStep * velocityOld[j+1] + halfStepSquared * forceOld[j+1];
                bas.posZ[i] += timeStep * velocityOld[j+2] + halfStepSquared * forceOld[j+2];

                j += 3;
            }

        }
        
        //bas.resetSimulationCell();

        try
        {
            
            //convergeDimer(mpi, bas, rVector, eigenValue, forcePar, forcePerp, totalEnergyNew, shortRangeCut, verletShell,
            //          fld, nbrList, useNbrList, dimerParameters, status, doEwald, outStream);
        }
        catch(const std::invalid_argument& e)
        {
            status.setStatusFailed();
            free_dvector(velocityOld);
            free_dvector(velocityNew);
            free_dvector(forceOld);
            free_dvector(forcePar);
            free_dvector(forcePerp);
            return;
        }
        catch(const std::runtime_error& e)
        {
            status.setStatusFailed();
            free_dvector(velocityOld);
            free_dvector(velocityNew);
            free_dvector(forceOld);
            free_dvector(forcePar);
            free_dvector(forcePerp);
            return;
        }
        /*if (!convergedLanzos)
        {
            status.setStatusFailed();
            free_dvector(velocityOld);
            free_dvector(velocityNew);
            free_dvector(forceOld);
            free_dvector(forcePar);
            free_dvector(forcePerp);
            return;
        } */
            
        if (eigenValue >= dimerParameters.maxEigenvalue || goPerpendicular == true)
        {     
            for (int i = 0; i < vectorSize; i++)
                velocityNew[i] = velocityOld[i] + halfStep * forcePerp[i] + halfStep * forceOld[i];   
            
        }
        else
        {
            for (int i = 0; i < vectorSize; i++)
                velocityNew[i] = velocityOld[i] + halfStep * forcePar[i] + halfStep * forceOld[i];
        }
        

        //!!FIRE Update
        vectorNormMax(forcePar, fNormPar, fMaxPar, vectorSize);
        vectorNormMax(forcePerp, fNormPerp, fMaxPerp, vectorSize);
        vectorNormMax(velocityNew, vNorm, vMax, vectorSize);



        if (dimerParameters.debug && it % dimerParameters.print == 0)
        {
            outStream << "\n FIRE: iteration " << it << " " << totalEnergyNew << " " << eigenValue << std::endl;
            if (eigenValue < 0.0)
                outStream << " FIRE: gnorm parallel " << fNormPar << " force max " << fMaxPar << " energy difference " 
                          << fabs(totalEnergyNew - totalEnergyOld) << std::endl;
            else
                outStream << " FIRE: gnorm perp " << fNormPerp << " force max " << fMaxPerp << " energy difference " 
                          << fabs(totalEnergyNew - totalEnergyOld) << std::endl;
        }

        // it is recomended that convergence is tested here
        deltaE = fabs(totalEnergyOld - totalEnergyNew);
        //if (deltaE < dimerParameters.energyTol || (fNorm / vectorSize) < dimerParameters.normTol || fMax < dimerParameters.forceTol)

        if (it > 5) //only allow if a reasonable number of iterations complete
        {
            if ((fNormPar + fNormPerp) < dimerParameters.normTol || (fMaxPar + fMaxPerp) < dimerParameters.forceTol || deltaE < dimerParameters.energyTol)
            {

                if (dimerParameters.debug)
                {

                    outStream << "\n FIRE converged : delta E (eV) = " << deltaE << std::endl;
                    outStream << "\n FIRE converged : delta E (internal) " << deltaE << " tolerance " 
                              << dimerParameters.energyTol << std::endl;                

                }

                convergedFire = true;
                if (eigenValue <= dimerParameters.maxEigenvalue)
                    status.setStatusSuccess();

                finalenergy = totalEnergyNew;
                free_dvector(velocityOld);
                free_dvector(velocityNew);
                free_dvector(forceOld);
                free_dvector(forcePar);
                free_dvector(forcePerp);
                return;
            }
        }

        // F1 of program
        if (eigenValue >= dimerParameters.maxEigenvalue || goPerpendicular == true)
        {
            P = dotProduct(forcePerp, velocityNew, vectorSize);


            OneLessAlpha = 1.0 - alpha;
            Norm = alpha / (fNormPerp * vNorm);

            // this is part F2 of the paper
            for (i = 0; i < vectorSize; i++)
                velocityOld[i] = OneLessAlpha * velocityNew[i] + forcePerp[i] * Norm;
        }
        else
        {
            P = dotProduct(forcePar, velocityNew, vectorSize);


            OneLessAlpha = 1.0 - alpha;
            Norm = alpha / (fNormPar * vNorm);

            // this is part F2 of the paper
            for (i = 0; i < vectorSize; i++)
                velocityOld[i] = OneLessAlpha * velocityNew[i] + forcePar[i] * Norm;
        }
        
        if(P > 0)  // downhill move
        {
            if (nStep > dimerParameters.N_min)
            {   //F3

                timeStep = fmin(timeStep * dimerParameters.stepIncrease, timeStepMax); // increase the time step
                alpha *= dimerParameters.alphaDecrease;                  // decrease alpha
            }

            nStep += 1;

        }
        else // uphill move (P <= 0)
        {   //F4

            nStep = 0;                        // reset counter for the number of steps between N_min
            timeStep *= dimerParameters.stepDecrease;         // decrease timestep
            alpha = dimerParameters.alphaStart;               // reset alpha to start

            for (i = 0; i < vectorSize; i++)
                velocityOld[i] = 0.0;              // freeze system

        }

        //can now update energies
        totalEnergyOld = totalEnergyNew;

	    #ifndef ARCHER2
        /*if (mpi.checkForInterrupt(outStream))
        {
            status.setStatusInterrupt();

            #ifdef DEBUG
	        std::cerr << "\n interrupt in FIRESaddle " << mpi.getRank() << " " << mpi.getWorkGroupRank() << " " << mpi.getGroupID();
            std::cerr.flush();
            #endif

            free_dvector(velocityOld);
            free_dvector(velocityNew);
            free_dvector(forceOld);
            free_dvector(forcePar);
            free_dvector(forcePerp);
            return;
        } */
	    #endif

    }
    
    

    free_dvector(velocityOld);
    free_dvector(velocityNew);
    free_dvector(forceOld);
    free_dvector(forcePar);
    free_dvector(forcePerp);
}


void Dimer::calculateForces(const MPIComms& mpi, double* forceVector, Basis& bas, Field* fld, double& totalEnergy)
{
    double
        *forceX, 
        *forceY, 
        *forceZ,
        stress[6] = {0.0};

    Energy
        eng;

    forceX = alloc_dvector(bas.numberOfAtoms, "Dimer: force", 0.0);
    forceY = alloc_dvector(bas.numberOfAtoms, "Dimer: force", 0.0);
    forceZ = alloc_dvector(bas.numberOfAtoms, "Dimer: force", 0.0);

    //jap added in false just to get it to compile with my version
    fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceX, forceY, forceZ, stress, bas.latVector, bas.rcpVector, bas.charge, 
                            bas.atmLabel, bas.frozen, bas.numberOfAtoms, eng);


    totalEnergy = eng.getTotalEnergy();

    // sum the forces on free-moving particles to find net force at centre/end of dimer
    int j = 0;
    for (int i = 0; i < bas.numberOfAtoms; i++)
    {
        if (bas.frozen[i] == 0)  //i.e. it is free to move
        {
            forceVector[j] = forceX[i];
            forceVector[j+1] = forceY[i];
            forceVector[j+2] = forceZ[i];
            j += 3;
        }
    }

    //free memory
    free_dvector(forceX);
    free_dvector(forceY);
    free_dvector(forceZ);

}


            //rotForce1 = dotProduct(rotForce, tVector, 3);
            
            // obtain new unit and rotational vectors by rotating original vectors
            // about plane formed by original translational and rotational vectors
/*
            rVector2[0] = rVector[0] + (rVector[0] * cosdeltheta + tVector[0] * sindeltheta);
            rVector2[1] = rVector[1] + (rVector[1] * cosdeltheta + tVector[1] * sindeltheta);
            rVector2[2] = rVector[2] + (rVector[2] * cosdeltheta + tVector[2] * sindeltheta);
        
            // vectorLen = vectorNorm(rVector2, 3);
            // rVectorLen = 1.0/vectorLen;
            // rVector2[0] *= rVectorLen;
            // rVector2[1] *= rVectorLen;
            // rVector2[2] *= rVectorLen;
            
            tVector[0] = tVector[0] + (rVector[0] * cosdeltheta + tVector[0] * sindeltheta);
            tVector[1] = tVector[1] + (rVector[1] * cosdeltheta + tVector[1] * sindeltheta);
            tVector[2] = tVector[2] + (rVector[2] * cosdeltheta + tVector[2] * sindeltheta);
        
            // vectorLen = vectorNorm(tVector, 3);
            // rVectorLen = 1.0/vectorLen;
            // tVector[0] *= rVectorLen;
            // tVector[1] *= rVectorLen;
            // tVector[2] *= rVectorLen;
            
            // copy initial basin (again) to give original particle positions
            // and apply new unit vector to particle positions for new position 1
        
            bas1 = bas;
            //makeDisplacement(mpi, bas1, rVector2, dimerParameters.dimerDisplacement, outStream);
        
            // calculate translational forces for non-frozen particles in position 1
        
            calculateForces(mpi, oneForce, bas1, nbrList, fld, bas1Energy, shortRangeCut, verletShell, doEwald);

            // calculate translational forces for non-frozen particles in position 2
            // (using approximation from initial and position 1 forces)
        
            twoForce[0] = 2.0 * zeroForce[0] - oneForce[0];
            twoForce[1] = 2.0 * zeroForce[1] - oneForce[1];
            twoForce[2] = 2.0 * zeroForce[2] - oneForce[2];
        
            // calculate rotational forces at each end of dimer
        
            oneFdotN2 = dotProduct(oneForce, rVector2, 3);
            twoFdotN2 = dotProduct(twoForce, rVector2, 3);
        
            rotoneForce[0] = oneForce[0] - oneFdotN2 * rVector2[0];
            rotoneForce[1] = oneForce[1] - oneFdotN2 * rVector2[1];
            rotoneForce[2] = oneForce[2] - oneFdotN2 * rVector2[2];
            rottwoForce[0] = twoForce[0] - twoFdotN2 * rVector2[0];
            rottwoForce[1] = twoForce[1] - twoFdotN2 * rVector2[1];
            rottwoForce[2] = twoForce[2] - twoFdotN2 * rVector2[2];

            // find overall scaled rotational force
        
            rotForce[0] = (rotoneForce[0] - rottwoForce[0]) * rDimerDisplace;
            rotForce[1] = (rotoneForce[1] - rottwoForce[1]) * rDimerDisplace;
            rotForce[2] = (rotoneForce[2] - rottwoForce[2]) * rDimerDisplace;

            //rotForce2 = dotProduct(rotForce, tVector, 3);

            // find magnitude of mean rotational force (average between two unit vectors)
            // and derivative w.r.t. angle
        
            //forceScale = 0.5 * (rotForce1 + rotForce2);
            //kurvature = rdtheta * (rotForce2 - rotForce1);
        
            // calculate angle to rotate second dimer to obtain minimum dimer energy
            // (add half pi if derivative is negative)
        
            dtheta = -0.5 * (atan(2.0 * forceScale / kurvature) + dimerParameters.deltaTheta);
            if (kurvature < 0.0) dtheta += halfpi;

            // calculate eigenvalue of rotational curvature (valid if minimum mode found)
        
            eigenValue = 0.5 * rDimerDisplace * (twoFdotN1 - oneFdotN1);
        
            // if not yet close to minimum mode (i.e. if rotational force still
            // too high), rotate second dimer to start over

            foundMode = !(fabs(forceScale) > dimerParameters.forceTol);
            if(!foundMode) {
                sindtheta = sin(dtheta);
                cosdtheta = cos(dtheta);
                rVector[0] = rVector2[0] + (rVector2[0] * cosdtheta + tVector[0] * sindtheta);
                rVector[1] = rVector2[1] + (rVector2[1] * cosdtheta + tVector[1] * sindtheta);
                rVector[2] = rVector2[2] + (rVector2[2] * cosdtheta + tVector[2] * sindtheta);
                // vectorLen = vectorNorm(rVector, 3);
                // rVectorLen = 1.0/vectorLen;
                // rVector[0] *= rVectorLen;
                // rVector[1] *= rVectorLen;
                // rVector[2] *= rVectorLen;
            }
        
            miniter++;
    
        } while (!foundMode && miniter < dimerParameters.N_min);
    
        // adjust eigenvalue if ran out of rotations for minimum mode
    
        if(!foundMode) {
            kurvature = 0.5 * rDimerDisplace * (twoFdotN2 - oneFdotN2);
            rotforce2mag = vectorNorm(rotForce, 3);
            eigenValue = kurvature - 0.5 * rotforce2mag * tan(dtheta - 0.5 * dimerParameters.deltaTheta);
        }

        // now translate dimer along lowest curvature mode to get to saddle point:
        // start by calculating modified force acting in direction of step vector
    
        oneFdotN1 = dotProduct(zeroForce, rVector, 3);
        if(eigenValue>0.0) {
            oneForce[0] = -oneFdotN1 * rVector[0];
            oneForce[1] = -oneFdotN1 * rVector[1];
            oneForce[2] = -oneFdotN1 * rVector[2];
        }
        else {
            oneForce[0] = zeroForce[0] - 2.0 * zeroFdotN * rVector[0];
            oneForce[1] = zeroForce[1] - 2.0 * zeroFdotN * rVector[1];
            oneForce[2] = zeroForce[2] - 2.0 * zeroFdotN * rVector[2];
        }
    
    // find direction of step vector
    
        vectorLen = vectorNorm(oneForce, 3);
        rVectorLen = 1.0/vectorLen;
    
        rVector2[0] = oneForce[0] * rVectorLen;
        rVector2[1] = oneForce[1] * rVectorLen;
        rVector2[2] = oneForce[2] * rVectorLen;
    
    // shift particles along small distance in direction of step vector
    // and calculate forces for that point
    
        bas1 = bas;
       // makeDisplacement(mpi, bas1, rVector2, dimerParameters.lineDisplacement, outStream);
       // calculateForces(mpi, newForce, bas1, nbrList, fld, bas1Energy, shortRangeCut, verletShell, doEwald);

        // calculate modified force for second point
    
        newFdotN = dotProduct(newForce, rVector, 3);

        if(eigenValue>0.0) {
            twoForce[0] = -newFdotN * rVector[0];
            twoForce[1] = -newFdotN * rVector[1];
            twoForce[2] = -newFdotN * rVector[2];
        }
        else {
            twoForce[0] = newForce[0] - 2.0 * newFdotN * rVector[0];
            twoForce[1] = newForce[1] - 2.0 * newFdotN * rVector[1];
            twoForce[2] = newForce[2] - 2.0 * newFdotN * rVector[2];
        }
    
        // find force magnitude and curvature along step vector
    
        oneFdotN1 = dotProduct(oneForce, rVector2, 3);
        twoFdotN1 = dotProduct(twoForce, rVector2, 3);
    
        forceScale = 0.5 * (oneFdotN1 + twoFdotN1);
        kurvature = (twoFdotN1 - oneFdotN1) * rLineDisplace;
    
        // find distance to saddle point: if greater than maximum
        // value (i.e. still some way away), set to this value instead
    
        deltaX = (-forceScale / kurvature + 0.5 * dimerParameters.lineDisplacement);
        if(deltaX > dimerParameters.maxDeltaX) deltaX = dimerParameters.maxDeltaX;
    
        // shift particles to new position ready for next iteration,
        // then calculate new forces and energy
 
        //makeDisplacement(mpi, bas, rVector2, deltaX, outStream);
        calculateForces(mpi, newForce, bas, nbrList, fld, totalEnergyNew, shortRangeCut, verletShell, doEwald);
        
        // check for convergence with energy - if achieved, report and finish
        
        if(dimerParameters.debug && iter % dimerParameters.print == 0 && grpRank == 0)
        {
            outStream << "\n Dimer: iteration " << iter << " " << totalEnergyNew << " " << deltaX << std::endl;
            outStream << " Dimer: force " << forceScale << " curvature " << kurvature << " energy difference "
                << fabs(totalEnergyNew - totalEnergyOld) << std::endl;
        }
        
        deltaE = fabs(totalEnergyOld - totalEnergyNew);
        
        if (deltaE < dimerParameters.energyTol)
        {
            if (dimerParameters.debug && grpRank == 0)
            {
                outStream << "\n Dimer converged : delta E (eV) = " << deltaE << std::endl;
                outStream << "\n Dimer converged : delta E (internal) " << deltaE << " tolerance " 
                          << dimerParameters.energyTol << std::endl;
            }
            convergedDimer = true;
            status.setStatusSuccess();
            finalenergy = totalEnergyNew;
            return;
        }
        
        // if not yet converged, update energies and start new iteration
        
        totalEnergyOld = totalEnergyNew;
        
    }
    */