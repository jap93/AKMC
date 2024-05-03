#include "Relax.h"

//////////////////////////////////////////////////////////////////////
// construction/destruction
//////////////////////////////////////////////////////////////////////

Relax::Relax()
{


}

Relax::~Relax()
{

}


void Relax::RelaxStructure(const MPIComms& mpi, Basis& bas, Field* fld, 
                           int natoms, const RelaxCntrl& rCntrl, Energy& eng, Status& status, std::ofstream& outStream)
{
    if (rCntrl.verletIntegration)
        euler = false;

    if (rCntrl.method == "fire")
    {
        FIRERelax(mpi, bas, fld, natoms, rCntrl, status, eng, outStream);
    }
    else if (rCntrl.method == "fire2")
    {
        FIRE2Relax(mpi, bas, fld, natoms, rCntrl, status, eng, outStream);
    }
    else
    {
        fld->calculateConfig(mpi, bas.posX, bas.posY, bas.posZ, bas.latVector, bas.rcpVector, bas.charge, bas.atmLabel, 
                             bas.frozen, natoms, eng);

        //by default assume that it is successful but need to get some info back eventually
        status.setStatusSuccess(); 
    }

}

void Relax::FIRERelax(const MPIComms& mpi, Basis& bas, Field* fld, int natoms, const RelaxCntrl& rCntrl, Status& status, 
                      Energy& eng, std::ofstream& outStream)
{
    int
        it,
        nStep = 1;

    double
        totalEnergyOld,
        totalEnergyNew,
        timeStep,
        timeStepMax,
        halfStep,
        halfStepSquared,
        alpha = rCntrl.alphaStart,
        OneLessAlpha,
        alphaForce,
        deltaE,
        fNorm,
        vNorm,
        fMax,
        vMax;

    double
        P;

    double
        stress[6] = {0.0},
        *velocityX, // velocity for deltaT/2
        *velocityY,
        *velocityZ,
        *forceNewX,
        *forceNewY,
        *forceNewZ;


    //Set FIRE parameters
    timeStep =  rCntrl.initTimeStep;
    timeStepMax = rCntrl.maxTimeStep;
    

    velocityX = alloc_dvector(natoms, "FIRE", 0.0);
    velocityY = alloc_dvector(natoms, "FIRE", 0.0);
    velocityZ = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewX = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewY = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewZ = alloc_dvector(natoms, "FIRE", 0.0);



    //calulate new forces
    fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceNewX, forceNewY, forceNewZ, stress, bas.latVector, bas.rcpVector, bas.charge, bas.atmLabel, 
                         bas.frozen, natoms, eng);
    //bas.freezeForces(forceNewX, forceNewY, forceNewZ);

    totalEnergyOld = eng.getTotalEnergy();

    for (int i = 0; i < natoms; i++)
    {
        if(bas.frozen[i] == 1)
        {
            forceNewX[i] = 0.0;
            forceNewY[i] = 0.0;
            forceNewZ[i] = 0.0;
        }
    }

    vectorNormMax(forceNewX, forceNewY, forceNewZ, fNorm, fMax, natoms);

    if (rCntrl.debug)
    {
        outStream << "\n FIRE: initial total energy       " << std::endl;
        eng.printEnergy(0, outStream);
        outStream << " FIRE: gnorm " << fNorm << " force max " << fMax << std::endl;
        outStream.flush();
    }

    
    halfStep = 0.5 * timeStep * eng.getTimeStepUnit();
    halfStepSquared = 0.5 * timeStep * timeStep;

    for (it = 1; it <= rCntrl.maxIter; it++)
    {

         // F1 of program
        P = dotProduct(forceNewX, forceNewY, forceNewZ, velocityX, velocityY, velocityZ, natoms);

        OneLessAlpha = 1.0 - alpha;
            //!!FIRE Update
        vectorNormMax(forceNewX, forceNewY, forceNewZ, fNorm, fMax, natoms);
        vectorNormMax(velocityX, velocityY, velocityZ, vNorm, vMax, natoms);

        alphaForce = alpha * vNorm / fNorm;

        if(P > 0)  // downhill move
        {
            // this is part F2 of the paper
            for (int i = 0; i < natoms; i++)
            {
                if(bas.frozen[i] == 0)
                {
                    velocityX[i] = OneLessAlpha * velocityX[i] + (alphaForce * forceNewX[i]);
                    velocityY[i] = OneLessAlpha * velocityY[i] + (alphaForce * forceNewY[i]);
                    velocityZ[i] = OneLessAlpha * velocityZ[i] + (alphaForce * forceNewZ[i]);
                }
            }

            if (nStep > rCntrl.N_min)
            {   //F3

                timeStep = fmin(timeStep * rCntrl.stepIncrease, timeStepMax); // increase the time step
                alpha *= rCntrl.alphaDecrease;                  // decrease alpha
            }

            nStep += 1;

        }
        else // uphill move (P <= 0)
        {   //F4
            nStep = 0;                        // reset counter for the number of steps between N_min
            timeStep *= rCntrl.stepDecrease;         // decrease timestep
            alpha = rCntrl.alphaStart;               // reset alpha to start

            for (int i = 0; i < natoms; i++)
            {
                velocityX[i] = 0.0;              // freeze system
                velocityY[i] = 0.0; 
                velocityZ[i] = 0.0; 
            }

        }

        if (euler == true)
        {
            double fact = timeStep * eng.getTimeStepUnit();
            for (int i = 0; i < natoms; i++)
            {
            // advance coordinates
                if(bas.frozen[i] == 0)
                {
                    bas.posX[i] += timeStep * velocityX[i];
                    bas.posY[i] += timeStep * velocityY[i];
                    bas.posZ[i] += timeStep * velocityZ[i];
                    velocityX[i] += fact * forceNewX[i];
                    velocityY[i] += fact * forceNewY[i];
                    velocityZ[i] += fact * forceNewZ[i];
                }
            } 

        }
        else
        {
            halfStep = 0.5 * timeStep * eng.getTimeStepUnit();  // NB include force conversion factor
            //1: v (t + 1/2 t ) <- v (t ) + 1/2 t· F (x (t ))/ m
            for (int i = 0; i < natoms; i++)
            {
                if(bas.frozen[i] == 0)
                {
                    velocityX[i] += halfStep * forceNewX[i];
                    velocityY[i] += halfStep * forceNewY[i];
                    velocityZ[i] += halfStep * forceNewZ[i];
                }
                else
                {
                    velocityX[i] = 0.0;
                    velocityY[i] = 0.0;
                    velocityZ[i] = 0.0;
                }
            }
        
            //3: x (t + t ) <- x (t ) + t ·v (t + 1/2 t )
            for (int i = 0; i < natoms; i++)
            {
                // advance coordinates
                if(bas.frozen[i] == 0)
                {
                    bas.posX[i] += timeStep * velocityX[i];
                    bas.posY[i] += timeStep * velocityY[i];
                    bas.posZ[i] += timeStep * velocityZ[i];
                }
            } 

        }
           
        //calulate new forces
        fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceNewX, forceNewY, forceNewZ, stress, bas.latVector, bas.rcpVector, bas.charge, bas.atmLabel, 
                         bas.frozen, natoms, eng);
        bas.freezeForces(forceNewX, forceNewY, forceNewZ);

        if (!euler)
        {
             //full time step 6: v (t + t ) <- v (t + 1/2 t ) + 1/2 t ·F (x (t + t ))/m
            for (int i = 0; i < natoms; i++)
            {
                if(bas.frozen[i] == 0)
                {
                    velocityX[i] += halfStep * forceNewX[i];
                    velocityY[i] += halfStep * forceNewY[i];
                    velocityZ[i] += halfStep * forceNewZ[i];
                }
                else
                {
                    velocityX[i] = 0.0;
                    velocityY[i] = 0.0;
                    velocityZ[i] = 0.0;
                    forceNewX[i] = 0.0;
                    forceNewY[i] = 0.0;
                    forceNewZ[i] = 0.0;
                }
            }
        }

        centreOfMass(velocityX, velocityY, velocityZ, natoms);


        totalEnergyNew = eng.getTotalEnergy();

       
        vectorNormMax(forceNewX, forceNewY, forceNewZ, fNorm, fMax, natoms);

        // convergence is tested here
        deltaE = fabs(totalEnergyOld - totalEnergyNew);

        if (rCntrl.debug && it % rCntrl.print == 0)
        {
            outStream << std::dec << std::scientific << std::setprecision(10)  << "\n FIRE: iteration " << it << " " << totalEnergyNew 
                      << " " << eng.getEnergyUnit() << " timestep " 
                      << timeStep << std::endl;
            outStream << std::dec << std::scientific << std::setprecision(10)  << " FIRE: gnorm " << fNorm << " force max " << fMax << " energy difference " <<  
                      fabs(totalEnergyNew - totalEnergyOld) << std::endl;
        }

        if (fNorm < rCntrl.normTol || fMax < rCntrl.forceTol) // || deltaE < rCntrl.energyTol)
        {

            if (rCntrl.debug)
            {
                outStream << std::dec << std::scientific << std::setprecision(10)  << "\n FIRE converged : delta E " << deltaE << std::endl;
                outStream << std::dec << std::scientific << std::setprecision(10)  << " FIRE converged : max force " << fMax << " tolerance " << rCntrl.forceTol << std::endl;
                outStream << std::dec << std::scientific << std::setprecision(10)  << " FIRE converged : gnorm " << fNorm << " tolerance " << rCntrl.normTol << std::endl;
                outStream.flush();
            }

            status.setStatusSuccess();
        
            free_dvector(velocityX);
            free_dvector(velocityY);
            free_dvector(velocityZ);
            free_dvector(forceNewX);
            free_dvector(forceNewY);
            free_dvector(forceNewZ);

            return;

        }

        //can now update energies
        totalEnergyOld = totalEnergyNew;

    }

    //should only get here if too many iterations are required
    status.setStatusFailed();

    free_dvector(velocityX);
    free_dvector(velocityY);
    free_dvector(velocityZ);
    free_dvector(forceNewX);
    free_dvector(forceNewY);
    free_dvector(forceNewZ);
    
}

void Relax::FIRE2Relax(const MPIComms& mpi, Basis& bas, Field* fld, int natoms, const RelaxCntrl& rCntrl, Status& status, 
                      Energy& eng, std::ofstream& outStream)
{
    int
        it,
        nUpHill = 0,
        nStep = 1;

    double
        totalEnergyOld,
        totalEnergyNew,
        timeStep,
        timeStepMax,
        timeStepMin,
        halfStep,
        halfStepSquared,
        alpha = rCntrl.alphaStart,
        OneLessAlpha,
        alphaForce,
        deltaE,
        fNorm,
        vNorm,
        fMax,
        vMax;

    double
        P;

    double
        stress[6] = {0.0},
        *velocityX, // velocity for deltaT/2
        *velocityY,
        *velocityZ,
        *forceNewX,
        *forceNewY,
        *forceNewZ;


    //Set FIRE parameters
    timeStep =  rCntrl.initTimeStep;
    timeStepMax = rCntrl.maxTimeStep;
    timeStepMin = rCntrl.minTimeStep;
    

    velocityX = alloc_dvector(natoms, "FIRE", 0.0);
    velocityY = alloc_dvector(natoms, "FIRE", 0.0);
    velocityZ = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewX = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewY = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewZ = alloc_dvector(natoms, "FIRE", 0.0);



    //calulate new forces
    fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceNewX, forceNewY, forceNewZ, stress, bas.latVector, bas.rcpVector, bas.charge, bas.atmLabel, 
                         bas.frozen, natoms, eng);

    //bas.freezeForces(forceNewX, forceNewY, forceNewZ);

    totalEnergyOld = eng.getTotalEnergy();

    for (int i = 0; i < natoms; i++)
    {
        if(bas.frozen[i] == 1)
        {
            forceNewX[i] = 0.0;
            forceNewY[i] = 0.0;
            forceNewZ[i] = 0.0;
        }
    }

    vectorNormMax(forceNewX, forceNewY, forceNewZ, fNorm, fMax, natoms);

    if (rCntrl.debug)
    {
        outStream << "\n FIRE: initial total energy       " << std::endl;
        eng.printEnergy(0, outStream);
        outStream << " FIRE: gnorm " << fNorm << " force max " << fMax << std::endl;
        outStream.flush();
    }

    
    halfStep = 0.5 * timeStep * eng.getTimeStepUnit();
    halfStepSquared = 0.5 * timeStep * timeStep;

    for (it = 1; it <= rCntrl.maxIter; it++)
    {

         // F1 of program
        P = dotProduct(forceNewX, forceNewY, forceNewZ, velocityX, velocityY, velocityZ, natoms);

        OneLessAlpha = 1.0 - alpha;
            //!!FIRE Update
        vectorNormMax(forceNewX, forceNewY, forceNewZ, fNorm, fMax, natoms);
        vectorNormMax(velocityX, velocityY, velocityZ, vNorm, vMax, natoms);

        alphaForce = alpha * vNorm / fNorm;

        if(P > 0)  // downhill move
        {
            // this is part F2 of the paper. I am a bit confused here as apendix A indicates it should go in MD step?
            /*for (i = 0; i < natoms; i++)
            {
                if(bas.frozen[i] == 0)
                {
                    velocityX[i] = OneLessAlpha * velocityX[i] + (alphaForce * forceNewX[i]);
                    velocityY[i] = OneLessAlpha * velocityY[i] + (alphaForce * forceNewY[i]);
                    velocityZ[i] = OneLessAlpha * velocityZ[i] + (alphaForce * forceNewZ[i]);
                }
            } */

            if (nStep > rCntrl.N_min)
            {   //F3

                timeStep = fmin(timeStep * rCntrl.stepIncrease, timeStepMax); // increase the time step
                alpha *= rCntrl.alphaDecrease;                  // decrease alpha
            }

            nStep += 1;
            nUpHill = 0;

        }
        else // uphill move (P <= 0)
        {   //F4
            nStep = 0;                        // reset counter for the number of steps between N_min
            nUpHill += 1;                     // increment the maximum number of uphill moves

            if (nUpHill > rCntrl.N_maxUphill)
            {
                status.setStatusFailed();
                free_dvector(velocityX);
                free_dvector(velocityY);
                free_dvector(velocityZ);
                free_dvector(forceNewX);
                free_dvector(forceNewY);
                free_dvector(forceNewZ);

                return;
            }

            if (rCntrl.initialDelay)
            {

                if (it >= rCntrl.N_delaySteps)
                {
                    timeStep *= rCntrl.stepDecrease;         // decrease timestep
                    if (timeStep < timeStepMin)
                        timeStep = timeStepMin;

                    alpha = rCntrl.alphaStart;               // reset alpha to start
                }
            }
            else
            {
                timeStep *= rCntrl.stepDecrease;         // decrease timestep
                alpha = rCntrl.alphaStart;               // reset alpha to start
            }
            
            // put positions back ie correct for up hill motion
            for (int i = 0; i < natoms; i++)
            {
            
                if(bas.frozen[i] == 0)
                {
                    bas.posX[i] -= 0.5 * timeStep * velocityX[i];
                    bas.posY[i] -= 0.5 * timeStep * velocityY[i];
                    bas.posZ[i] -= 0.5 * timeStep * velocityZ[i];
                }
            }

            for (int i = 0; i < natoms; i++)
            {
                velocityX[i] = 0.0;              // freeze system
                velocityY[i] = 0.0; 
                velocityZ[i] = 0.0; 
            }

        }

        if (euler == true)  // this is semi-implicit as indicated in appendix A
                            // explicit is not used as indicated
        {
            double fact = timeStep * eng.getTimeStepUnit();
            for (int i = 0; i < natoms; i++)
            {
            // advance coordinates
                if(bas.frozen[i] == 0)
                {
                    if (P > 0.0)
                    {
                        velocityX[i] = OneLessAlpha * velocityX[i] + (alphaForce * forceNewX[i]);
                        velocityY[i] = OneLessAlpha * velocityY[i] + (alphaForce * forceNewY[i]);
                        velocityZ[i] = OneLessAlpha * velocityZ[i] + (alphaForce * forceNewZ[i]);
                    }
                    
                    bas.posX[i] += timeStep * velocityX[i];
                    bas.posY[i] += timeStep * velocityY[i];
                    bas.posZ[i] += timeStep * velocityZ[i];
                
                    velocityX[i] += fact * forceNewX[i];
                    velocityY[i] += fact * forceNewY[i];
                    velocityZ[i] += fact * forceNewZ[i];
                }
            } 

        }
        else
        {
            halfStep = 0.5 * timeStep * eng.getTimeStepUnit();  // NB include force conversion factor
            //1: v (t + 1/2 t ) <- v (t ) + 1/2 t· F (x (t ))/ m
            for (int i = 0; i < natoms; i++)
            {
                if(bas.frozen[i] == 0)
                {
                    if (P > 0.0)
                    {
                        velocityX[i] = OneLessAlpha * velocityX[i] + (alphaForce * forceNewX[i]);
                        velocityY[i] = OneLessAlpha * velocityY[i] + (alphaForce * forceNewY[i]);
                        velocityZ[i] = OneLessAlpha * velocityZ[i] + (alphaForce * forceNewZ[i]);
                    }
                    
                    velocityX[i] += halfStep * forceNewX[i];
                    velocityY[i] += halfStep * forceNewY[i];
                    velocityZ[i] += halfStep * forceNewZ[i];
                }
                else
                {
                    velocityX[i] = 0.0;
                    velocityY[i] = 0.0;
                    velocityZ[i] = 0.0;
                }
            }
        
            //3: x (t + t ) <- x (t ) + t ·v (t + 1/2 t )
            for (int i = 0; i < natoms; i++)
            {
                // advance coordinates
                if(bas.frozen[i] == 0)
                {
                    bas.posX[i] += timeStep * velocityX[i];
                    bas.posY[i] += timeStep * velocityY[i];
                    bas.posZ[i] += timeStep * velocityZ[i];
                }
            } 

        }
           
        //calulate new forces
        fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceNewX, forceNewY, forceNewZ, stress, bas.latVector, bas.rcpVector, bas.charge, bas.atmLabel, 
                         bas.frozen, natoms, eng);
        bas.freezeForces(forceNewX, forceNewY, forceNewZ);

        if (!euler)
        {
             //full time step 6: v (t + t ) <- v (t + 1/2 t ) + 1/2 t ·F (x (t + t ))/m
            for (int i = 0; i < natoms; i++)
            {
                if(bas.frozen[i] == 0)
                {
                    velocityX[i] += halfStep * forceNewX[i];
                    velocityY[i] += halfStep * forceNewY[i];
                    velocityZ[i] += halfStep * forceNewZ[i];
                }
                else
                {
                    velocityX[i] = 0.0;
                    velocityY[i] = 0.0;
                    velocityZ[i] = 0.0;
                }
            }
        }

        centreOfMass(velocityX, velocityY, velocityZ, natoms);


        totalEnergyNew = eng.getTotalEnergy();

       
        vectorNormMax(forceNewX, forceNewY, forceNewZ, fNorm, fMax, natoms);

        // convergence is tested here
        deltaE = fabs(totalEnergyOld - totalEnergyNew);

        if (rCntrl.debug && it % rCntrl.print == 0)
        {
            outStream << std::dec << std::scientific << std::setprecision(10)  << "\n FIRE2: iteration " << it << " " << totalEnergyNew 
                      << " " << eng.getEnergyUnit() << " timestep " 
                      << timeStep << std::endl;
            outStream << std::dec << std::scientific << std::setprecision(10)  << " FIRE2: gnorm " << fNorm << " force max " << fMax << " energy difference " <<  
                      fabs(totalEnergyNew - totalEnergyOld) << std::endl;
        }

        if (fNorm < rCntrl.normTol || fMax < rCntrl.forceTol) // || deltaE < rCntrl.energyTol)
        {

            if (rCntrl.debug)
            {
                outStream << std::dec << std::scientific << std::setprecision(10)  << "\n FIRE2 converged : delta E " << deltaE << std::endl;
                outStream << std::dec << std::scientific << std::setprecision(10)  << " FIRE2 converged : max force " << fMax << " tolerance " << rCntrl.forceTol << std::endl;
                outStream << std::dec << std::scientific << std::setprecision(10)  << " FIRE2 converged : gnorm " << fNorm << " tolerance " << rCntrl.normTol << std::endl;
                outStream.flush();
            }

            status.setStatusSuccess();
        
            free_dvector(velocityX);
            free_dvector(velocityY);
            free_dvector(velocityZ);
            free_dvector(forceNewX);
            free_dvector(forceNewY);
            free_dvector(forceNewZ);

            return;

        }

        //can now update energies
        totalEnergyOld = totalEnergyNew;

    }

    //should only get here if too many iterations are required
    status.setStatusFailed();

    free_dvector(velocityX);
    free_dvector(velocityY);
    free_dvector(velocityZ);
    free_dvector(forceNewX);
    free_dvector(forceNewY);
    free_dvector(forceNewZ);
    
}
/*
void Relax::FIRERelax(const MPIComms& mpi, Basis& bas, Field* fld, int natoms, const RelaxCntrl& rCntrl, Status& status, 
                      Energy& eng, std::ofstream& outStream)
{
    int
        i,
        it,
        nStep = 1;

    double
        totalEnergyOld,
        totalEnergyNew,
        timeStep,
        timeStepMax,
        halfStep,
        halfStepSquared,
        alpha = rCntrl.alphaStart,
        OneLessAlpha,
        deltaE,
        fNorm,
        vNorm,
        fMax,
        vMax;

    double
        P;

    double
        *velocityX, // velocity for deltaT/2
        *velocityY,
        *velocityZ,
        *forceOldX,
        *forceOldY,
        *forceOldZ,
        *forceNewX,
        *forceNewY,
        *forceNewZ;


    //Set FIRE parameters
    timeStep =  rCntrl.initTimeStep;
    timeStepMax = rCntrl.maxTimeStep;
    

    velocityX = alloc_dvector(natoms, "FIRE", 0.0);
    velocityY = alloc_dvector(natoms, "FIRE", 0.0);
    velocityZ = alloc_dvector(natoms, "FIRE", 0.0);
    forceOldX = alloc_dvector(natoms, "FIRE", 0.0);
    forceOldY = alloc_dvector(natoms, "FIRE", 0.0);
    forceOldZ = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewX = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewY = alloc_dvector(natoms, "FIRE", 0.0);
    forceNewZ = alloc_dvector(natoms, "FIRE", 0.0);



    //calulate new forces
    fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceNewX, forceNewY, forceNewZ, bas.latVector, bas.rcpVector, bas.charge, bas.atmLabel, 
                         bas.frozen, natoms, eng);

    //bas.freezeForces(forceNewX, forceNewY, forceNewZ);

    totalEnergyOld = eng.getTotalEnergy();

    for (int i = 0; i < natoms; i++)
    {
        if(bas.frozen[i] == 1)
        {
            forceNewX[i] = 0.0;
            forceNewY[i] = 0.0;
            forceNewZ[i] = 0.0;
        }
    }

    vectorNormMax(forceNewX, forceNewY, forceNewZ, fNorm, fMax, natoms);

    if (rCntrl.debug)
    {
        outStream << "\n FIRE: initial total energy       " << std::endl;
        eng.printEnergy(outStream);
        outStream << " FIRE: gnorm " << fNorm << " force max " << fMax << std::endl;
        outStream.flush();
    }

    
    halfStep = 0.5 * timeStep * eng.getTimeStepUnit();
    halfStepSquared = 0.5 * timeStep * timeStep;

    for (it = 1; it <= rCntrl.maxIter; it++)
    {

        for (i = 0; i < natoms; i++)
        {
            forceOldX[i] = forceNewX[i];
            forceOldY[i] = forceNewY[i];
            forceOldZ[i] = forceNewZ[i];
        }

        //1: v (t + 1/2 t ) <- v (t ) + 1/2 t· F (x (t ))/ m
        for (int i = 0; i < natoms; i++)
        {
            if(bas.frozen[i] == 0)
            {
                velocityX[i] += halfStep * forceOldX[i];
                velocityY[i] += halfStep * forceOldY[i];
                velocityZ[i] += halfStep * forceOldZ[i];
            }
            else
            {
                velocityX[i] = 0.0;
                velocityY[i] = 0.0;
                velocityZ[i] = 0.0;
            }
        }
        //i = 0;
        //outStream << "before " << bas.rpos[i].x << " " << bas.rpos[i].y  << " " << bas.rpos[i].z << std::endl;
        centreOfMass(velocityX, velocityY, velocityZ, natoms);

    

        //!!FIRE Update
        vectorNormMax(forceOldX, forceOldY, forceOldZ, fNorm, fMax, natoms);
        vectorNormMax(velocityX, velocityY, velocityZ, vNorm, vMax, natoms);



        // F1 of program
        P = dotProduct(forceOldX, forceOldY, forceOldZ, velocityX, velocityY, velocityZ, natoms);

        OneLessAlpha = 1.0 - alpha;
        //Norm = alpha / (fNorm * vNorm);

        if(P > 0)  // downhill move
        {
            // this is part F2 of the paper
            for (i = 0; i < natoms; i++)
            {
                if(bas.frozen[i] == 0)
                {
                    velocityX[i] = OneLessAlpha * velocityX[i] + (alpha * forceOldX[i] * vNorm / fNorm);
                    velocityY[i] = OneLessAlpha * velocityY[i] + (alpha * forceOldY[i] * vNorm / fNorm);
                    velocityZ[i] = OneLessAlpha * velocityZ[i] + (alpha * forceOldZ[i] * vNorm / fNorm);
                }
            }

            if (nStep > rCntrl.N_min)
            {   //F3

                timeStep = fmin(timeStep * rCntrl.stepIncrease, timeStepMax); // increase the time step
                alpha *= rCntrl.alphaDecrease;                  // decrease alpha
            }

            nStep += 1;

        }
        else // uphill move (P <= 0)
        {   //F4
            nStep = 0;                        // reset counter for the number of steps between N_min
            timeStep *= rCntrl.stepDecrease;         // decrease timestep
            alpha = rCntrl.alphaStart;               // reset alpha to start

            for (i = 0; i < natoms; i++)
            {
                velocityX[i] = 0.0;              // freeze system
                velocityY[i] = 0.0; 
                velocityZ[i] = 0.0; 
            }

        }


        //3: x (t + t ) <- x (t ) + t ·v (t + 1/2 t )
        for (i = 0; i < natoms; i++)
        {
            // advance coordinates
            if(bas.frozen[i] == 0)
            {
                bas.posX[i] += timeStep * velocityX[i];
                bas.posY[i] += timeStep * velocityY[i];
                bas.posZ[i] += timeStep * velocityZ[i];
            }
        } 
           
        //calulate new forces
        fld->calculateForces(mpi, bas.posX, bas.posY, bas.posZ, forceNewX, forceNewY, forceNewZ, bas.latVector, bas.rcpVector, bas.charge, 
                            bas.atmLabel, bas.frozen, natoms, eng);
        //bas.freezeForces(forceNewX, forceNewY, forceNewZ);

        totalEnergyNew = eng.getTotalEnergy();

        //full time step 6: v (t + t ) <- v (t + 1/2 t ) + 1/2 t ·F (x (t + t ))/m
        for (int i = 0; i < natoms; i++)
        {
            if(bas.frozen[i] == 0)
            {
                velocityX[i] += halfStep * forceNewX[i];
                velocityY[i] += halfStep * forceNewY[i];
                velocityZ[i] += halfStep * forceNewZ[i];
            }
            else
            {
                velocityX[i] = 0.0;
                velocityY[i] = 0.0;
                velocityZ[i] = 0.0;
                forceNewX[i] = 0.0;
                forceNewY[i] = 0.0;
                forceNewZ[i] = 0.0;
            }
        }
    

        // convergence is tested here
        deltaE = fabs(totalEnergyOld - totalEnergyNew);

        if (rCntrl.debug && it % rCntrl.print == 0)
        {
            outStream << "\n FIRE: iteration " << it << " " << totalEnergyNew << " " << eng.getEnergyUnit() << " timestep " 
                      << timeStep << std::endl;
            outStream << " FIRE: gnorm " << fNorm << " force max " << fMax << " energy difference " 
                      << fabs(totalEnergyNew - totalEnergyOld) << std::endl;
        }

        if (fNorm < rCntrl.normTol || fMax < rCntrl.forceTol) // || deltaE < rCntrl.energyTol)
        {

            if (rCntrl.debug)
            {
                outStream << "\n FIRE converged : delta E " << deltaE << std::endl;
                outStream << " FIRE converged : max force " << fMax << " tolerance " << rCntrl.forceTol << std::endl;
                outStream << " FIRE converged : gnorm " << fNorm << " tolerance " << rCntrl.normTol << std::endl;
                outStream.flush();
            }

            status.setStatusSuccess();
        
            free_dvector(velocityX);
            free_dvector(velocityY);
            free_dvector(velocityZ);
            free_dvector(forceOldX);
            free_dvector(forceOldY);
            free_dvector(forceOldZ);
            free_dvector(forceNewX);
            free_dvector(forceNewY);
            free_dvector(forceNewZ);

            return;

        }

        //the MD part of FIRE - note mass is set to 1.0 for FIRE
        halfStep = 0.5 * timeStep * eng.getTimeStepUnit();
        halfStepSquared = 0.5 * timeStep * timeStep;

        //can now update energies
        totalEnergyOld = totalEnergyNew;

    }

    //should only get here if too many iterations are required
    status.setStatusFailed();

    free_dvector(velocityX);
    free_dvector(velocityY);
    free_dvector(velocityZ);
    free_dvector(forceOldX);
    free_dvector(forceOldY);
    free_dvector(forceOldZ);
    free_dvector(forceNewX);
    free_dvector(forceNewY);
    free_dvector(forceNewZ);
    
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// utility functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void Relax::vectorNormMax(double* vecX, double* vecY, double* vecZ, double& vecNorm, double& vecMax, int n)
{
    int
        i;
   
    double
        vecSq;

    vecNorm = 0.0;
    vecMax = 0.0;

    for (i = 0; i < n; i++)
    {
        vecSq = vecX[i] * vecX[i] + vecY[i] * vecY[i] + vecZ[i] * vecZ[i];
   
        if (vecSq > vecMax) vecMax = vecSq;
      
        vecNorm += vecSq;
    }
   
    vecNorm = sqrt(vecNorm);
    vecMax = sqrt(vecMax);
}

double Relax::dotProduct(double* x1, double* y1, double* z1, double* x2, double* y2, double* z2, int n)
{
    int 
        i;

    double
        dot = 0.0;

    for (i = 0; i < n; i++)
        dot += x1[i] * x2[i] + y1[i] * y2[i] + z1[i] * z2[i];

    return dot;
}

void Relax::centreOfMass(double* velocityX, double* velocityY, double* velocityZ, int natoms)
{
    int
        i;

    double
        vx = 0.0,
        vy = 0.0,
        vz = 0.0;

    for (i = 0; i < natoms; i++)
    {

        vx += velocityX[i];
        vy += velocityY[i];
        vz += velocityZ[i];

    }

    vx = vx / double(natoms);
    vy = vy / double(natoms);
    vz = vz / double(natoms);

    for (i = 0; i < natoms; i++)
    {

        velocityX[i] -= vx;
        velocityY[i] -= vy;
        velocityZ[i] -= vz;

    }

}
