/** 
*/
#include "KineticMonteCarlo.h"

double getRandomNumber(void);

KineticMonteCarlo::KineticMonteCarlo()
{
    
}

KineticMonteCarlo::~KineticMonteCarlo()
{
    
}

/** 
The main MD driver function
*/
void KineticMonteCarlo::runKMC(const MPIComms& mpi, int& restartIteration, double& restartTime, const Species& spec, Field* fld,  
                               const JobControl& job, Basis& basin1, std::ofstream& outStream)
{
    int
        ievent,
        iteration = 0,
        numProcs = mpi.getNumProcs(),
        worldRank = mpi.getRank();

    double
        actEnergy = 0.0,
        timeIncrease = 0.0,
        kmcTime = 0.0,
        ranNumber = 0.0,
        totalRate = 0.0,
        tempRate = 0.0;

    bool
        converged = false,
        prime = false;

    double
        basin1Energy;

    Energy
        relaxedEnergy;

    Event
        ev;

    Relax
        rel;

    Status
        relStatus;

    std::vector<int>
        excludedList;

    std::vector<Event>
        eventList;

    std::ofstream
        archiveStream;

    if (worldRank == 0)
        prime = true;

    if (job.kmcParameters.superBasin)
    {
        //excludedList.reserve(basin1.numberOfAtoms);
        for (int i = 0; i < basin1.numberOfAtoms; i++)
            excludedList.push_back(0);
    }

    if (worldRank == 0)
    {
        if (job.restart)
        {
            archiveStream.open("archive", std::ios::out | std::ofstream::app);
        }
        else
        {
            archiveStream.open("archive", std::ios::out);
        }  
    }

    if (worldRank == 0 && job.restart)
    {
        iteration = restartIteration;
        kmcTime = restartTime / job.kmcParameters.kmcPreFactor;

        outStream << "\n\n *****************************************************************************************************" << std::endl;
        outStream << "\n Restart undertaken at iteration : " << iteration << " for a kmc time of : " << kmcTime << std::endl;
        outStream << " *****************************************************************************************************" << std::endl;
        outStream.flush();
    }

    //start the number of iterations-NB the initial position on basin1 has been read
    //in by all cores
    while (!converged) //for (iteration = 1; ; iteration++)
    {
        iteration++;
        
        outStream << "\n\n *****************************************************************************************************" << std::endl;
        outStream << " Kinetic Monte Carlo iteration : " << iteration << std::endl;
        outStream << " *****************************************************************************************************" << std::endl;
        outStream.flush();
            
        //relax the atoms to get 1st basin energy - even prime has a copy of basin1Energy
        basin1.unFreezeAtoms();
        basin1.freezeAtomTypes(spec, job.frozenTypes);  //if atom types have been frozen in main job description
        
        rel.RelaxStructure(mpi, basin1, fld, basin1.numberOfAtoms, 
                               job.minParameters, relaxedEnergy, relStatus, outStream); 

        basin1Energy = relaxedEnergy.getTotalEnergy(); 

        //check that it did actually relax and interrupts
        if (relStatus.getStatus() == FAILED)
        {
            outStream << "\n failed to relax atoms at the start of the KMC iteration" << std::endl;
            outStream.flush();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        } 
        
        outStream  << std::scientific << std::setprecision(10) 
                   << "\n total energy of initial relaxation " << basin1Energy << std::endl;
        
        //auto startTime = system_clock::now(); 
        auto startTime = MPI_Wtime();

        // create the list of events
        if (numProcs == 1)
        {
            calculateRateConstantSingleCore(mpi, spec, fld, job, basin1, relaxedEnergy, eventList, iteration, outStream);
        }
        else
        {
            calculateRateConstantExperimental(mpi, spec, fld, job, basin1, basin1Energy, eventList, excludedList, iteration, outStream);
        }
                        
        //at this point eventList will be different for each processor so we need to amalgamate the rates
        if (prime)
        {
            
            totalRate = 0.0;
            //for (int it = 0 ; it < eventList.size(); it++)
            for (Event v : eventList)
            {
                //tempRate = eventList[it].getRate();
                //actEnergy = eventList[it].getActivationEnergy();
                tempRate = v.getRate();
                actEnergy = v.getActivationEnergy();
                outStream << "\n activation energy, rate " << actEnergy << "    " << tempRate << std::endl;
                totalRate += tempRate;
            }

            ranNumber = getRandomNumber() * totalRate;
            
            tempRate = 0.0;
            for (int it = 0 ; it < eventList.size(); it++)
            {
                tempRate += eventList[it].getRate();
                if (tempRate > ranNumber) // this is to selct event
                {
                    ievent = it;
                    break;
                }
            }
            outStream << "\n\n kmc iteration : " << iteration << " chosen event : " << ievent << std::endl;

            timeIncrease = -1.0 * calculateTime(totalRate);
            
            kmcTime += timeIncrease;
            // write out data
            outStream << std::scientific << std::setprecision(5) << " total rate : " << totalRate << std::endl;
            outStream << std::scientific << std::setprecision(5) << " kinetic monte time step : " << timeIncrease << " now :" << kmcTime << std::endl;

            //auto finishTime = system_clock::now();
            //auto diff = finishTime - startTime;
            //outStream << " elapsed time for this iteration : " << duration<double, std::milli>(diff).count() / 1000 << " seconds" << std::endl;
            auto finishTime = MPI_Wtime();
            auto diff = finishTime - startTime;
            outStream << std::scientific << std::setprecision(5) << " elapsed time for this iteration : " << diff << " seconds" << std::endl;
            startTime = finishTime;

            if (worldRank == 0 && job.kmcParameters.writeEvents)
                writeEventList(mpi, spec, basin1, eventList, kmcTime * job.kmcParameters.kmcPreFactor, basin1Energy, ievent, iteration); 

            //update the excluded list before creating new basin 1 position
            if (job.kmcParameters.superBasin)
                updateExcludedList(basin1, eventList[ievent], excludedList, job.kmcParameters.kmcBasinRadius, job.kmcParameters.historySize, outStream);
	    
            
            //set the new basin to the positions of the selected event
            eventList[ievent].setPositions(basin1, outStream);

            basin1Energy = eventList[ievent].getBasin2Energy();

            //either delete events or recycle events
            clearEventList(eventList, basin1.latVector, basin1.rcpVector, ievent, job.kmcParameters.recycleRadius, job.kmcParameters.recycle, 
                            job.kmcParameters.useCluster, job.kmcParameters.recycleType, job.kmcParameters.recyclePercentage, basin1, 
			               excludedList, job.kmcParameters.kmcBasinRadius, job.kmcParameters.superBasin, outStream);
     
        }

        //update all cores with new positions
        broadCastPositions(mpi, basin1, outStream);
       
        if (iteration >= job.kmcParameters.kmcSteps)
            converged = true;
        

        if(worldRank == 0)
        {
            std::ofstream restartStream;
            restartStream.open("restart", std::ofstream::out);
            restartStream << "restart " << basin1.title << std::endl;
            restartStream << " 0   3 " << basin1.numberOfAtoms << " " << basin1.numberOfAtoms << std::endl;

	        basin1.dumpBasis(mpi, spec, kmcTime * job.kmcParameters.kmcPreFactor, basin1Energy, iteration, restartStream);
	        restartStream.flush();
	        restartStream.close();

            basin1.dumpBasis(mpi, spec, kmcTime * job.kmcParameters.kmcPreFactor, basin1Energy, iteration, archiveStream);
            archiveStream.flush();

	        if (job.kmcParameters.superBasin)
	        {
	            std::ofstream exListStream;
                exListStream.open("excludedlist", std::ofstream::out);

                writeExcludedList(excludedList, exListStream);
		        exListStream.flush();
		        exListStream.close();
            }

        }

            //make all processors wait here
        int ierr;
        #ifdef DEBUG
        outStream << std::fixed << std::setprecision(10) << "\n waiting at mpi_barrier on iteration " << iteration << std::endl;
        outStream.flush();
        #endif
        ierr = MPI_Barrier(MPI_COMM_WORLD);
        
    } // done enough iterations

    //write the restart file
    if(worldRank == 0)
    {
        std::ofstream restartStream;
        restartStream.open("restart", std::ofstream::out);
        restartStream << "restart " << basin1.title << std::endl;
        restartStream << " 0   3 " << basin1.numberOfAtoms << " " << basin1.numberOfAtoms << std::endl;
	    basin1.dumpBasis(mpi, spec, kmcTime * job.kmcParameters.kmcPreFactor, basin1Energy, iteration - 1, restartStream);
	    restartStream.flush();
	    restartStream.close();

        basin1.dumpBasis(mpi, spec, kmcTime * job.kmcParameters.kmcPreFactor, basin1Energy, iteration, archiveStream);
        archiveStream.flush();
        archiveStream.close();
    }

    //make all processors wait here
    int ierr;
    #ifdef DEBUG
    outStream << std::fixed << std::setprecision(10) << "\n waiting at mpi_barrier" << std::endl;
    outStream.flush();
    #endif
    ierr = MPI_Barrier(MPI_COMM_WORLD);

}

void KineticMonteCarlo::calculateRateConstant(const MPIComms& mpi, const Species& spec, Field* fld, const JobControl& job, Basis& basin1, 
                                              double& basin1Energy, std::vector<Event>& eventList, const std::vector<int>& excludedList, 
                                              int iteration, std::ofstream& outStream)
{
    int
        source = -1,
        command,
        numBusy = 0,
        numEvents = 0,
        grpRank = mpi.getWorkGroupRank(),
        grpID = mpi.getGroupID(),
        numGrps = mpi.getNumOfWorkGroups();

    double
        actEnergy = 0.0;

    Energy
        startingEnergy,
        basin2Energy,
        saddleEnergy;
    
    bool
        continueSearch = false,
        interrupted = false,
        foundNewBasin = false;

    Basis
        saddle,
        basin2;

    double
        centre[3];

    std::vector<int>
        atmList;

    std::vector<int>
        workerList;

    Event
        ev;

    Relax
        rel;

    MPI_Status
        statMPI;

    SaddleSearch  // TODO convert to a unique pointer
        *method = nullptr;

    if (job.kmcParameters.kmcMethod == "art")
    {
        method = new ART(); //new ART();
    }
    else if (job.kmcParameters.kmcMethod == "dimer")
    {
        method = new Dimer(); //new Dimer();
    }
    else
    {
        outStream << "\n The kmc method is not recognised" << std::endl;
        outStream.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    if (numGrps > job.kmcParameters.kmcImages)
    {
        outStream << "\n The number of events required is less than the number of work groups" << std::endl;
        outStream.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    
    if (grpID == 0)  // the parent process - choreographs the calculation of rates
    {

        bool finish = false;
        int attemptedSearches = 0;
        int failedSearches = 0;
        int msgTag = 0;

        //busyNode = alloc_ivector(numGrps, "allocate the number of groups ", 0);

        for (int ig = 1; ig <= numGrps; ig++)
        {
            if (grpRank == 0)
            {
                mpi.sendCommandToGroup(STARTSEARCH, ig-1, outStream);  // initiate all groups
                attemptedSearches++;
                numBusy++;
                //busyNode[ig] = 1;
            }
            
        }
              
        while (!finish) //keep on waiting for data to be sent back
        {

            #ifdef DEBUG
            outStream << "\n\n waiting to get message back from worker";
            outStream.flush();
            #endif
            mpi.checkForMessage(command, source, msgTag, outStream);

            #ifdef DEBUG
            outStream << "\n parent has received command " << command << " proc " << source << " grp " << mpi.getGroupID(source);
            outStream.flush();
            #endif

            //look for incoming message - first should contain the size of data, data (or error code if negative) and tag should be the iteration
            int grp = mpi.getGroupID(source);
            //busyNode[grp] = 0; // set the group to not busy
            if (command > SUCCESS) // SUCCESS == 0. command is the size of the message expected
            {
                #ifdef DEBUG
                outStream << "\n about to get event";
                outStream.flush();
                #endif

                //tag is the rank of the processor sending the message
                if (grpRank == 0)
                {
                    //receivePositions(mpi, ev, basin1.numberOfAtoms, job.kmcParameters.kmcPreFactor, job.kmcParameters.kmcTemperature,
                    //                          job.kmcBoltzmannConst, source, msgTag, outStream);
                    receiveEvent(mpi, ev, job.kmcParameters.kmcPreFactor, job.kmcParameters.kmcTemperature, command, 
                                source, msgTag, outStream);
                    eventList.push_back(ev);

                    outStream << "\n activation energy " << ev.getActivationEnergy()
                              << " from proc and work group " << source << " " << grp
                              << " the number of events " << eventList.size() << std::endl;
                    numBusy--;
                    
                }
                
            }
            else if (command == FAILED)
            {
                failedSearches++;
                if (grpRank == 0)
                    numBusy--;
            }

            mpi.barrierGroup();
            numEvents = eventList.size();

            if ((numBusy + numEvents) < job.kmcParameters.kmcImages)
            {
                #ifdef DEBUG
                outStream << "\n parent sending new command to group " << grp;
                outStream.flush();
                #endif
                if (grpRank == 0)
                {
                    mpi.sendCommandToGroup(STARTSEARCH, grp, outStream);  // initiate all groups
                    attemptedSearches++;
                    numBusy++;
                    //busyNode[grp] = 1;
                }
            }

            if (numEvents == job.kmcParameters.kmcImages)
                finish = true;

            mpi.barrierGroup();
        }   //finished collecting data
        
        //enough events have been found so interrupt all the processes
        
        if (grpRank == 0)
        {
            for (int ig = 1; ig <= numGrps; ig++)
            {
                mpi.sendCommandToGroup(TERMINATED, ig-1, outStream);  // terminate all groups
            }
        } 

        if (grpRank == 0)
        {
            outStream << "\n the number of searches started " << attemptedSearches << std::endl;
            outStream << "\n the number of failed searches  " << failedSearches << std::endl;
        }

        //free_ivector(busyNode);

    }
    else
    {
        

        #ifdef DEBUG
        outStream << "\n in worker" << std::endl;
        outStream.flush();
        #endif


        do
        {

            Status
                sadStatus;

            int
                saddleSearchStatus = 0;

            //start time
            auto startTime = system_clock::now();

            //else do another event search
            foundNewBasin = false;
            int command = mpi.waitForCommand();
            
            // tag = statMPI.MPI_TAG;
            source = statMPI.MPI_SOURCE;
            
            
            if (command == TERMINATED)
            {
                //if (method != nullptr)
                //{
                //    delete method;
                //}
                interrupted = true;
                //mpi.sendMessageToParent(TERMINATED);
                break;
            }
            else if (command == STARTSEARCH || command == CONTINUESEARCH)
            {	
                //if (mpi.checkForInterrupt(outStream))
                //{
                //    outStream << "\n interrupted in KineticMonteCarlo.cpp!" << std::endl;
                //    outStream.flush();
                //    goto exitFunction; 
                //}   
                         
                try
                {
                    outStream << "\n starting search " << std::endl;
                    outStream.flush();
                    sadStatus.setStatusFailed();
                    //find a transition state/saddle
	                method->findSaddlePoint(mpi, spec, job, fld, basin1, saddle, sadStatus, startingEnergy, saddleEnergy, 
                                            centre, continueSearch, outStream);
                }

                catch(const std::exception& e)
                {
                    outStream << e.what() << std::endl;
                }

                if (sadStatus.getStatus() == FAILED)
                {
                    mpi.sendMessageToParent(FAILED);
                    continue;
                }

                
                sadStatus.setStatusFailed();
                foundNewBasin = checkForNewBasin(mpi, spec, job, fld, basin1, basin2, saddle, sadStatus, basin2Energy, outStream);
                
                if (sadStatus.getStatus() == FAILED)
                {
                    mpi.sendMessageToParent(FAILED);
                    continue;
                }

                //calculate activation energy 
                actEnergy = (saddleEnergy.getTotalEnergy() - basin1Energy);

                if (foundNewBasin == false)
                {
                    mpi.sendMessageToParent(FAILED);
                }
                else if (actEnergy >= job.kmcParameters.kmcMinCap && actEnergy <= job.kmcParameters.kmcMaxCap && foundNewBasin == true)
                {
                    #if DEBUG
                    outStream << "\n about to send positions. saddlesearch status " << saddleSearchStatus << " interrupt " << interrupted;
                    outStream.flush();
                    #endif

                    int iter = 0;
                    //send the new positions back to rank 0 and ask for more work to do
                    interrupted = sendEvent(mpi, basin1, basin2, saddle, centre, basin1Energy, basin2Energy.getTotalEnergy(), 
                                            saddleEnergy.getTotalEnergy(), actEnergy, iter, outStream);

                    if (interrupted)
                    {
                        #ifdef DEBUG
                        outStream << "\n interrupted after send position!" << std::endl;
                        outStream.flush();
                        #endif

                        goto exitFunction;
                    }
                       
                    outStream << "\n\n*** valid activation energy " 
                              << "\n basin 1 energy      : " << basin1Energy
                              << "\n basin 2 energy      : " << basin2Energy.getTotalEnergy()
                              << "\n saddle point energy : " << saddleEnergy.getTotalEnergy()
                              << "\n activation energy   : " << actEnergy << std::endl;

                    auto finishTime = system_clock::now();
                    auto diff = finishTime - startTime;
                    outStream << "\n time to calculate event : " << duration<double, std::milli>(diff).count() / 1000 << " seconds" << std::endl;

                }
                else
                {
                    outStream << std::scientific << std::setprecision(5) << "\n invalid activation energy " << actEnergy << std::endl;
                    outStream << std::scientific << std::setprecision(5) << "\n limits are " << job.kmcParameters.kmcMinCap << " to " << job.kmcParameters.kmcMaxCap;
                    outStream.flush();
                    
                    mpi.sendMessageToParent(FAILED);
                    //send an error code back to the parent

                }
            }

        } while (!interrupted);

    }


 
exitFunction:
    #ifdef DEBUG
    outStream << "\n leaving calculate rate constant for proc " << worldRank;
    outStream.flush();
    #endif

    if (method != nullptr) delete method; //get rid of the memory.
    method = nullptr;
}
     
void KineticMonteCarlo::calculateRateConstantExperimental(const MPIComms& mpi, const Species& spec, Field* fld, const JobControl& job, Basis& basin1, 
                                              double& basin1Energy, std::vector<Event>& eventList, const std::vector<int>& excludedList, 
                                              int iteration, std::ofstream& outStream)
{
    int
        tag = 0,
        source = -1,
        command,
        numBusy = 0,
        numEvents = 0,
        grpRank = mpi.getWorkGroupRank(),
        grpID = mpi.getGroupID(),
        numGrps = mpi.getNumOfWorkGroups();

    double
        actEnergy = 0.0;

    Energy
        startingEnergy,
        basin2Energy,
        saddleEnergy;
    
    bool
        continueSearch = false,
        interrupted = false,
        foundNewBasin = false;

    Basis
        saddle,
        basin2;

    double
        centre[3];

    std::vector<int>
        atmList;

    std::vector<int>
        workerList;

    Event
        ev;

    Relax
        rel;

    MPI_Status
        statMPI;

    SaddleSearch  // TODO convert to a unique pointer
        *method = nullptr;

    if (job.kmcParameters.kmcMethod == "art")
    {
        method = new ART(); //new ART();
    }
    else if (job.kmcParameters.kmcMethod == "dimer")
    {
        method = new Dimer(); //new Dimer();
    }
    else
    {
        outStream << "\n The kmc method is not recognised" << std::endl;
        outStream.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    if (numGrps > job.kmcParameters.kmcImages)
    {
        outStream << "\n The number of events required is less than the number of work groups" << std::endl;
        outStream.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    
    if (grpID == 0)  // the parent process - choreographs the calculation of rates
    {

        bool finish = false;
        int attemptedSearches = 0;
        int failedSearches = 0;
        int terminatedSearches = 0;
        int msgTag = 0;

        for (int ig = 1; ig <= numGrps; ig++)
        {
            if (grpRank == 0)
            {
                outStream << "\n sending command to start " << ig-1 << "  " << numGrps;
                mpi.sendCommandToGroup(STARTSEARCH, ig-1, outStream);  // initiate all groups
                attemptedSearches++;               
            }
            numBusy++;
            int active = 1;
            workerList.push_back(active);
            
        }
              
        while (!finish) //keep on waiting for data to be sent back
        {

            #ifdef DEBUG
            outStream << "\n\n waiting to get message back from worker";
            outStream.flush();
            #endif
            mpi.checkForMessage(command, source, msgTag, outStream);

            #ifdef DEBUG
            outStream << "\n parent has received command " << command << " proc " << source << " grp " << mpi.getGroupID(source);
            outStream.flush();
            #endif

            //look for incoming message - first should contain the size of data, data (or error code if negative) and tag should be the iteration
            int grp = mpi.getGroupID(source);
            workerList[grp] = 0; // set the group to not busy
            numBusy--;

            if (command > SUCCESS) // SUCCESS == 0. command is the size of the message expected
            {
                #ifdef DEBUG
                outStream << "\n about to get event";
                outStream.flush();
                #endif

                //tag is the rank of the processor sending the message
                if (grpRank == 0)
                {
                    //receivePositions(mpi, ev, basin1.numberOfAtoms, job.kmcParameters.kmcPreFactor, job.kmcParameters.kmcTemperature,
                    //                          job.kmcBoltzmannConst, source, msgTag, outStream);
                    receiveEvent(mpi, ev, job.kmcParameters.kmcPreFactor, job.kmcParameters.kmcTemperature, command, 
                                source, msgTag, outStream);
                    eventList.push_back(ev);

                    outStream << "\n activation energy " << ev.getActivationEnergy()
                              << " from proc and work group " << source << " " << grp
                              << " the number of events " << eventList.size() << std::endl;
                   
                    
                }
                
            }
            else if (command == FAILED)
            {
                failedSearches++;
            }

            mpi.barrierGroup();
            numEvents = eventList.size();

            if (numEvents == job.kmcParameters.kmcImages)
                finish = true;

            if (finish == false)
            {
                //if ((numBusy + numEvents) < (job.kmcParameters.kmcImages + job.kmcParameters.kmcExcessSearches))
                if (numEvents < job.kmcParameters.kmcImages)
                {
                    #ifdef DEBUG
                    outStream << "\n parent sending new command to group " << grp;
                    outStream.flush();
                    #endif
                    if (grpRank == 0)
                    {
                        mpi.sendCommandToGroup(STARTSEARCH, grp, outStream);  // initiate all groups
                        attemptedSearches++; 
                    }
                
                    numBusy++;
                    workerList[grp] = 1;
                }
            }

           
            mpi.barrierGroup();
        }   //finished collecting data
        
        //enough events have been found so interrupt all the processes
        
        if (grpRank == 0)
        {
            for (int ig = 1; ig <= numGrps; ig++)
            {
                if (workerList[ig-1] == 1)
                {
                    #ifdef DEBUG
                    outStream << "\n sending final interrupted to active work group " << ig << std::endl;
                    #endif
                    terminatedSearches++;
                    
                }
                #ifdef DEBUG
                if (workerList[ig-1] == 0)
                {
                    outStream << "\n sending final interrupted to work group " << ig << std::endl;
                }
                #endif

                mpi.sendCommandToGroup(TERMINATED, ig-1, outStream);  // initiate all groups
                
            }
        } 

        if (grpRank == 0)
        {
            outStream << "\n the number of searches started " << attemptedSearches << std::endl;
            outStream << "\n the number of failed searches  " << failedSearches << std::endl;
            outStream << "\n the number of stopped searches " << terminatedSearches << std::endl;
        }

    }
    else
    {
        
        do
        {
            Status
                sadStatus;

            int
                saddleSearchStatus = 0;

            //start time
            auto startTime = system_clock::now();

            //else do another event search
            foundNewBasin = false;
            
            int command = mpi.waitForCommand();
            #ifdef DEBUG
            outStream << "\n command received " << command << std::endl;
            outStream.flush();
            #endif
            
            
            if (command == TERMINATED)
            {
                interrupted = true;
                //mpi.sendMessageToParent(TERMINATED);
                break;
            }
            else if (command == STARTSEARCH || command == CONTINUESEARCH)
            {	
                
                         
                try
                {
                    #ifdef DEBUG
                    outStream << "\n starting search " << std::endl;
                    outStream.flush();
                    #endif 

                    sadStatus.setStatusFailed();
                    //find a transition state/saddle
	                method->findSaddlePoint(mpi, spec, job, fld, basin1, saddle, sadStatus, startingEnergy, saddleEnergy, 
                                            centre, continueSearch, outStream);
                }

                catch(const std::exception& e)
                {
                    outStream << e.what() << std::endl;
                }

                if (sadStatus.getStatus() == FAILED)
                {
                    mpi.sendMessageToParent(FAILED);
                    continue;
                }

                //mpi.checkForInterrupt(sadStatus, outStream);

                if (sadStatus.getStatus() == TERMINATED)
                {
                    outStream << "\n interrupted during ART!" << std::endl;
                    outStream.flush();
                    interrupted = true;
                    break; 
                }

                
                sadStatus.setStatusFailed();
                foundNewBasin = checkForNewBasin(mpi, spec, job, fld, basin1, basin2, saddle, sadStatus, basin2Energy, outStream);
                
                if (sadStatus.getStatus() == FAILED)
                {
                    mpi.sendMessageToParent(FAILED);
                    continue;
                }

                //calculate activation energy 
                actEnergy = (saddleEnergy.getTotalEnergy() - basin1Energy);

                mpi.checkForInterrupt(sadStatus, outStream);

                if (sadStatus.getStatus() == TERMINATED)
                {
                    outStream << "\n interrupted in KineticMonteCarlo.cpp!" << std::endl;
                    outStream.flush();
                    interrupted = true;
                    break; 
                }

                if (foundNewBasin == false)
                {
                    mpi.sendMessageToParent(FAILED);
                }
                else if (actEnergy >= job.kmcParameters.kmcMinCap && actEnergy <= job.kmcParameters.kmcMaxCap && foundNewBasin == true)
                {
                    #if DEBUG
                    outStream << "\n about to send positions. saddlesearch status " << saddleSearchStatus << " interrupt " << interrupted;
                    outStream.flush();
                    #endif

                    int iter = 0;
                    //send the new positions back to rank 0 and ask for more work to do
                    interrupted = sendEvent(mpi, basin1, basin2, saddle, centre, basin1Energy, basin2Energy.getTotalEnergy(), 
                                            saddleEnergy.getTotalEnergy(), actEnergy, iter, outStream);

                    if (interrupted)
                    {
                        #ifdef DEBUG
                        outStream << "\n interrupted after before position!" << std::endl;
                        outStream.flush();
                        #endif

                        interrupted = true;
                        break;
                    }
                       
                    outStream << "\n\n*** valid activation energy " 
                              << "\n basin 1 energy      : " << basin1Energy
                              << "\n basin 2 energy      : " << basin2Energy.getTotalEnergy()
                              << "\n saddle point energy : " << saddleEnergy.getTotalEnergy()
                              << "\n activation energy   : " << actEnergy << std::endl;

                    auto finishTime = system_clock::now();
                    auto diff = finishTime - startTime;
                    outStream << "\n time to calculate event : " << duration<double, std::milli>(diff).count() / 1000 << " seconds" << std::endl;

                }
                else
                {
                    outStream << std::scientific << std::setprecision(5) << "\n invalid activation energy " << actEnergy << std::endl;
                    outStream << std::scientific << std::setprecision(5) << "\n limits are " << job.kmcParameters.kmcMinCap << " to " << job.kmcParameters.kmcMaxCap;
                    outStream.flush();
                    
                    mpi.sendMessageToParent(FAILED);
                    //send an error code back to the parent

                }
            }

        } while (!interrupted);

    }

workerList.clear();
 
exitFunction:
    #ifdef DEBUG
    outStream << "\n leaving calculate rate constant for proc " << grpID;
    outStream.flush();
    #endif

    if (method != nullptr) delete method; //get rid of the memory.
    method = nullptr;
}

void KineticMonteCarlo::calculateRateConstantSingleCore(const MPIComms& mpi, const Species& spec, Field* fld, const JobControl& job, Basis& basin1, 
                                              Energy& basin1Energy, std::vector<Event>& eventList, int iteration, std::ofstream& outStream)
{
    int
        tag,
        dummy,
        ierr,
        numBusy = 0,
        numEvents = 0,
        numSteps = 1,
        rank = mpi.getRank(),
        numProcs = mpi.getNumProcs();

    double
        actEnergy = 0.0,
        centre[3];

    Energy
        basin2Energy,
        saddleEnergy;

    bool
        continueSearch = false,
        interrupted = false,
        prime = false,
        newBasin = false,
        foundNewBasin = false;

    Basis
        saddle,
        basin2;

    SaddleSearch  // TODO convert to a unique pointer
        *method;

    Status
        status;

    Event
        ev;

    Relax
        rel;

    if (rank == 0)
        prime = true;

    eventList.clear();

    if (job.kmcParameters.kmcMethod == "art")
    {
        method = new ART();
    }
    else if (job.kmcParameters.kmcMethod == "dimer")
    {
        method = new Dimer();
    }
    else
    {
        outStream << "\n The kmc method is not recognised" << std::endl;
        outStream.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    while (numEvents < job.kmcParameters.kmcImages) //keep on waiting for data to be sent back
    {

        //start time
        auto startTime = system_clock::now();
 
        //else do another event search
        foundNewBasin = false;

        do
        {

            try
            {
                //find a transition state/saddle
                method->findSaddlePoint(mpi, spec, job, fld, basin1, saddle, status, basin1Energy, saddleEnergy, centre, 
                                        continueSearch, outStream);
                outStream << "\n just finished saddle search " << status.getStatus() << std::endl;
            }
            catch(const std::exception& e)
            {
                outStream << e.what() << std::endl;
            }

            if (status.getStatus() == FAILED)
            {
                outStream << "\n search failed " << std::endl;
                continue;
            }

            status.setStatusFailed();
            foundNewBasin = checkForNewBasin(mpi, spec, job, fld, basin1, basin2, saddle, status, basin2Energy, outStream);
                
            if (status.getStatus() == FAILED)
            {
                continue;
            }

            //calculate activation energy
            actEnergy = (saddleEnergy.getTotalEnergy() - basin1Energy.getTotalEnergy());

            if (actEnergy >= job.kmcParameters.kmcMinCap && actEnergy <= job.kmcParameters.kmcMaxCap)
            {

                foundNewBasin = true;
                outStream << "\n\n*** valid activation energy "
                              << "\n basin 1 energy      : " << basin1Energy.getTotalEnergy()
                              << "\n basin 2 energy      : " << basin2Energy.getTotalEnergy()
                              << "\n saddle point energy : " << saddleEnergy.getTotalEnergy()
                              << "\n activation energy   : " << actEnergy << std::endl;

                auto finishTime = system_clock::now();
                auto diff = finishTime - startTime;
                outStream << " time to calculate event : " << duration<double, std::milli>(diff).count() / 1000 << " seconds" << std::endl;

                float* vecf = alloc_fvector(6*saddle.numberOfAtoms+7, " singlecore: vecf");
                int* veci = alloc_ivector(saddle.numberOfAtoms, " singlecore veci");
                //copy the positions to the buffer
                for (int i = 0; i < saddle.numberOfAtoms; i++)
                {
                    vecf[6*i] = basin2.posX[i];
                    vecf[6*i+1] = basin2.posY[i];
                    vecf[6*i+2] = basin2.posZ[i];
                    vecf[6*i+3] = saddle.posY[i];
                    vecf[6*i+4] = saddle.posY[i];
                    vecf[6*i+5] = saddle.posZ[i];
                    veci[i] = saddle.globalNo[i];
                }
                int dim = 6 * saddle.numberOfAtoms;
                vecf[dim] = float(basin1Energy.getTotalEnergy());
                vecf[dim+1] = float(basin2Energy.getTotalEnergy());
                vecf[dim+2] = float(saddleEnergy.getTotalEnergy());
                vecf[dim+3] = float(actEnergy);  
    
                vecf[dim+4] = float(centre[0]);
                vecf[dim+5] = float(centre[1]);
                vecf[dim+6] = float(centre[2]);

                ev = Event(veci, vecf, saddle.numberOfAtoms, job.kmcParameters.kmcPreFactor, job.kmcParameters.kmcTemperature);

                eventList.push_back(ev);
                numEvents = eventList.size();

                free_ivector(veci);
                free_fvector(vecf);

            }
            else
            {
                outStream << "\n invalid activation energy " << actEnergy << std::endl;
            }



        } while (!foundNewBasin);

    }

    delete method; //get rid of the memory.
        
}

/*
 * advance monte carlo time
 */
double KineticMonteCarlo::calculateTime(double totalRate)
{
    double
        deltaT = 0.0,
        ranNumber = getRandomNumber();

    if (abs(totalRate) > 1.0e-18)
        deltaT = (1.0 / totalRate) * log(ranNumber);

    return deltaT;
}

void KineticMonteCarlo::clearEventList(std::vector<Event>& eventList, double* latVector, double* rcpVector, int chosenEvent, double radius, bool recycle, 
                                       bool global, int mode, double percent, Basis& basin, std::vector<int>& excludedList, double basinRadius,
				                       bool superBasin, std::ofstream& outStream)
{
    int
        numEvents = eventList.size();

    double
        *ac,
        *bc,
        choice = 0.0,
        ratio = 1.2,  // set to a value that is always accepted to initialise
        xx, yy, zz,
        rx, ry, rz, r;

    bool
        found = false;

    std::vector<Event> tmpList;

    if (mode == 2)
        ratio = percent / 100;

    if (recycle == true)  // if global events are done then the centre is always set to zero
    {
        // recycling of saddle points is done based on centre to centre distance
        ac = eventList[chosenEvent].getCentre();

	outStream << "\n centre point of chosen event " << ac[0] << " " << ac[1] << " " << ac[2] << std::endl;

        for (int i = 0; i < numEvents; i++)
        {
            if (i == chosenEvent)
                continue;

            bc = eventList[i].getCentre();

            //distance between the centres of events
            rx = ac[0] - bc[0];
            ry = ac[1] - bc[1];
            rz = ac[2] - bc[2];
            
            xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
            yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
            zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            r = sqrt(rx * rx + ry * ry + rz * rz);

            if (r >= radius)
            {

		if (superBasin)
		{
		    found = false;
		    //add the atoms that have just moved
                    //first work out atoms that have moved greater than the basin distance
                    std::vector<int> atmList = eventList[i].getDisplacedAtoms(basin, basinRadius);
                    for (int ex : excludedList)
                    {
			for (int atm : atmList)
			{
			    if (ex == atm)
				found = true;
			}
                    }
                }
                if (mode == 2)
                    choice = getRandomNumber();
                
                if (choice < ratio && found == false)
                    tmpList.push_back(eventList[i]);  //add the event to the temporary list
            }

        }

        //clear the event list and copy the temporary list back into it
        eventList.clear();
        for (int it = 0; it < tmpList.size(); it++)
        {
            eventList.push_back(tmpList[it]);
            ac = tmpList[it].getCentre();
            outStream << std::scientific << std::setprecision(7) 
                      << "\n saved event with centre " << ac[0] << "  " << ac[1] << "  " << ac[2] << std::endl;
        }

        outStream << "\n the number of events that have been recycled = " << tmpList.size();

        tmpList.clear();

    }
    else 
    {
        eventList.clear();
    }
}

bool KineticMonteCarlo::checkEventValidity(Basis& basin, Event& ev, const std::vector<int>& excludedList, int iteration, double basinRadius, std::ofstream& outStream)
{
    bool
        valid = true;

    //atoms that have displaced
    std::vector<int> atoms = ev.getDisplacedAtoms(basin, basinRadius);

    //check if the atom is in the excluded list
    for (int i = 0; i < atoms.size(); i++)
    {
        int atm = atoms[i];
        outStream << "\n checking atom " << atm << "in the excluded list";
        if (excludedList[atm] > 0)
        {
            outStream << "\n atom " << atm << "found in the excluded list";
            return false;  // the atom is in both the event and excluded list so can not move
        }

    }

    return valid;
}

void KineticMonteCarlo::updateExcludedList(Basis& basin, Event& ev, std::vector<int>& excludedList, double basinRadius, int history, std::ofstream& outStream)
{
    //decrease the counter on all atoms (provided they are frozen)
    for (int i = 0; i < excludedList.size(); i++)
    {
        if (excludedList[i] > 0)
            excludedList[i] -= 1;
    }

    //add the atoms that have just moved
    //first work out atoms that have moved greater than the basin distance
    std::vector<int> atmList = ev.getDisplacedAtoms(basin, basinRadius);

    //second add the history length to the atoms in the excluded list
    for (int i = 0; i < atmList.size(); i++)
    {
        int j = atmList[i];
        excludedList[j] += history;
        outStream << "\n atom " << j << " has been added to the excluded list " << excludedList[j];
    }

    //now remove from any saved events

}
void KineticMonteCarlo::broadCastPositions(const MPIComms& mpi, Basis& basin1, std::ofstream& outStream)
{
    int
        pid = mpi.getRank(),
        ierr;

    double
        *buffer;

    buffer = alloc_dvector(3 * basin1.numberOfAtoms, " broadcast basin", 0.0);

    if (pid == 0)
    {
        //copy the positions to the buffer
        for (int i = 0; i < basin1.numberOfAtoms; i++)
        {
            buffer[3*i] = basin1.posX[i];
            buffer[3*i+1] = basin1.posY[i];
            buffer[3*i+2] = basin1.posZ[i];
        }
    }

    #ifdef DEBUG
    outStream << "\n broadcasting new set of positions";
    outStream.flush();
    #endif
    //broadcast positions to all cores for the next iteration
    ierr = MPI_Bcast(buffer, 3 * basin1.numberOfAtoms, MPI_DOUBLE, 0, MPI_COMM_WORLD); //mpi.getCommunicator());
    //mpi.sumDoubleWorld(buffer, 3 * basin1.numberOfAtoms);

    if (ierr != 0)
    {
        #ifdef DEBUG
        std::cerr << "\n mpi error in bcast of new positions " << std::endl;
        std::cerr.flush();
        #endif
        
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    //copy back to basin and free memory
    for (int i = 0; i < basin1.numberOfAtoms; i++)
    {
        basin1.posX[i] = buffer[3*i];
        basin1.posY[i] = buffer[3*i+1];
        basin1.posZ[i] = buffer[3*i+2];
    }    

    free_dvector(buffer);
}

void KineticMonteCarlo::sendErrorCode(const MPIComms& mpi, int code, int tag, std::ofstream& outStream)
{

    int errCode = -1 * code;

    //mpi.sendIntVector(&errCode, 1, 0, tag, outStream);
    MPI_Send(&errCode, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);

}
bool KineticMonteCarlo::sendEvent(const MPIComms& mpi, Basis& basin1, Basis& basin2, Basis& saddle, double* centre, double basin1Energy, double basin2Energy, 
                                double saddleEnergy, double actEnergy, int tag, std::ofstream& outStream)
{
    int
        //tag = mpi.getRank(),  //tag is used to indicate where the positions have come from
        size = saddle.numberOfAtoms,
        destination = 0,
        ierr;

    Status
        stat;

    int
        *bufferi = nullptr;

    float
        *bufferf = nullptr;

    bool interrupt = false;
    
    mpi.checkForInterrupt(stat, outStream);
    if (stat.getStatus() == TERMINATED)
    {
        interrupt = true;
        return interrupt;
    } 
    
    //#ifdef DEBUG
    outStream << "\n sending event, size, destination, tag) : " << size << " " << destination << " " << tag << std::endl;
    outStream.flush();
    //#endif
    MPI_Send(&size, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);

    bufferi = alloc_ivector(size, " broadcast basin", 0);
    for (int i = 0; i < size; i++)
    {
        bufferi[i] = saddle.globalNo[i];
    }

    #ifdef DEBUG
    outStream << "\n sendEvent : sending global number " << tag;
    outStream.flush();
    #endif
    MPI_Send(bufferi, size, MPI_INT, destination, tag, MPI_COMM_WORLD);
    
    bufferf = alloc_fvector(6*size+7, " broadcast basin", 0.0);

    //copy the positions to the buffer
    int j = 0;
    for (int i = 0; i < size; i++)
    {
        bufferf[j] = basin2.posX[i];
        bufferf[j+1] = basin2.posY[i];
        bufferf[j+2] = basin2.posZ[i];
        bufferf[j+3] = saddle.posX[i];
        bufferf[j+4] = saddle.posY[i];
        bufferf[j+5] = saddle.posZ[i]; 

        j += 6;
        
    }

    int dim = size * 6;
    bufferf[dim] = basin1Energy;
    bufferf[dim+1] = basin2Energy;
    bufferf[dim+2] = saddleEnergy;
    bufferf[dim+3] = actEnergy;
    bufferf[dim+4] = centre[0];
    bufferf[dim+5] = centre[1];
    bufferf[dim+6] = centre[2];


    //broadcast
    #ifdef DEBUG
    outStream << "\n sendEvent : sending positions " << tag;
    outStream.flush();
    #endif

    MPI_Send(bufferf, 6*size+7, MPI_FLOAT, destination, tag, MPI_COMM_WORLD);

    
    free_fvector(bufferf);
    free_ivector(bufferi);

    return interrupt;
}

void KineticMonteCarlo::receivePositions(const MPIComms& mpi, Event& ev, int numberOfAtoms, double preFactor, double temperature, double boltz,
                                        int source, int msgTag, std::ofstream& outStream)
{
    int
        tag,
        err;

    double
        basin1Energy,
        basin2Energy,
        saddleEnergy,
        actEnergy,
        *buffer;

    MPI_Status
        status;

    buffer = alloc_dvector(9 * numberOfAtoms + 7, " receive position ", 0.0);

    //broadcast
    err = MPI_Recv(buffer, 9 * numberOfAtoms + 7, MPI_DOUBLE, MPI_ANY_SOURCE, msgTag, MPI_COMM_WORLD, &status);  // should be done in MPIComms

    if (err != MPI_SUCCESS)
    {
        int eclass;
        MPI_Error_class(err, &eclass);

        std::cerr << "\n KineticMonteCarlo::receivePositions " << err << " " << eclass << " rank " << mpi.getRank() << std::endl;
        std::cerr.flush();
        mpi.commsAbortWorld();
        exit(EXIT_FAILURE);
    }

    basin1Energy = buffer[9*numberOfAtoms];
    basin2Energy = buffer[9*numberOfAtoms+1];
    saddleEnergy = buffer[9*numberOfAtoms+2];
    actEnergy = buffer[9*numberOfAtoms+3];

    //tag = status.MPI_TAG;
    //Event(int* veci, float* vecf, int natoms, double prefactor, double temperature)
    //ev = Event(buffer, numberOfAtoms, actEnergy, preFactor, temperature, boltz, basin1Energy, basin2Energy, saddleEnergy);

    free_dvector(buffer);

    //return tag;
}

void KineticMonteCarlo::receiveEvent(const MPIComms& mpi, Event& ev, double preFactor, double temperature, int& size, 
                                    int& source, int& tag, std::ofstream& outStream)
{
    int
        *bufferi = nullptr;

    double
        basin1Energy,
        basin2Energy,
        saddleEnergy,
        actEnergy;

    float
        *bufferf = nullptr;

    MPI_Status
        status;

    bufferi = alloc_ivector(size, " receive position ", 0.0);
    bufferf = alloc_fvector(6*size+7, " receive position ", 0.0);


    #ifdef DEBUG
    outStream << "\n receiveEvent: about to receive a integer vector " << size << "  " << source << "  " << tag;
    outStream.flush();
    #endif
    MPI_Recv(bufferi, size, MPI_INT, source, tag, MPI_COMM_WORLD, &status);

    #ifdef DEBUG
    outStream << "\n receiveEvent: about to receive a double vector " << size << "  " << source << "  " << tag;
    outStream.flush();
    #endif
    MPI_Recv(bufferf, 6*size+7, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    
    #ifdef DEBUG
    outStream << "\n receiveEvent: calling event constructor ";
    outStream.flush();
    #endif
    ev = Event(bufferi, bufferf, size, preFactor, temperature);
    #ifdef DEBUG
    outStream << "\n receiveEvent: created event ";
    outStream.flush();
    #endif

    free_ivector(bufferi);
    free_fvector(bufferf);

}



void KineticMonteCarlo::writeEventList(const MPIComms& mpi, const Species& spec, Basis& basin1, std::vector<Event>& eventList, double simulationTime, 
                                       double basin1Energy, int chosenEvent, int iteration)
{
    int
        lb;

    Element
        ele;

    std::ofstream
        outStream;

    std::string
        name,
        fileName = "events." + std::to_string(iteration);

    outStream.open(fileName, std::ofstream::out);

    //write out basin 1 energy and positions 
    outStream << "events   " << eventList.size() << " chosen event " << chosenEvent << " time " << simulationTime << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
              << std::setprecision(10) << std::setw(15) << "basin1 " << basin1Energy << std::endl;

    basin1.printVectors(outStream);

    for (int i = 0; i < basin1.numberOfAtoms; i++)
    {            
        lb = basin1.atmLabel[i];
        ele = spec.getSpecies(lb);
        name = ele.name;

        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
            << std::setprecision(10)
            << std::setw(5) << name << "    " << i << "\n"
            << std::setw(15) << basin1.posX[i]
            << std::setw(15) << basin1.posY[i]
            << std::setw(15) << basin1.posZ[i] << std::endl;
    }
    
    //write out the activation energy, saddle positions and energy and basin2 positions and energy
    for (int i = 0; i < eventList.size(); i++)
    {
        outStream << " event " << i << std::endl;
        outStream << " activation energy " << eventList[i].getActivationEnergy() << std::endl;

        eventList[i].writeCentre(outStream);

        eventList[i].writeSaddle(spec, basin1.atmLabel, outStream);

        eventList[i].writeBasin2(spec, basin1.atmLabel, outStream);
    }

    outStream.flush();
    outStream.close();
}

void KineticMonteCarlo::writeExcludedList(std::vector<int> excludedList, std::ofstream& outStream)
{
    for (int i = 0; i < excludedList.size(); i++)
    {
        outStream << i << "     " << excludedList[i] << std::endl;
    }
}

bool KineticMonteCarlo::checkForNewBasin(const MPIComms& mpi, const Species& spec, const JobControl& job, Field* fld, Basis& basin1, Basis& basin2, 
                                         Basis& saddle, Status& relStatus, Energy& basin2Energy, std::ofstream& outStream)
{
    int
        rank = mpi.getWorkGroupRank();
    double
        rx, ry, rz,
        xx, yy, zz,
        deltaPosX[basin1.numberOfAtoms],
        deltaPosY[basin1.numberOfAtoms],
        deltaPosZ[basin1.numberOfAtoms];

    bool 
        newBasinFound = false;

    Relax
        rel;

    relStatus.setStatusFailed();
    //copy the saddle point and give the positions a kick before relaxation to help prevent it falling back into basisn1
    //actually this should be basin1 otherwise just the local region is relaxed - but make sure to use the saddle positions further down
    //allocate memory for positions here - ie after the return
    basin2 = saddle;
    for (int i = 0; i < basin1.numberOfAtoms; i++)
    {
        rx = saddle.posX[i] - basin1.posX[i];  //take into account periodic boundary
        ry = saddle.posY[i] - basin1.posY[i];
        rz = saddle.posZ[i] - basin1.posZ[i];

        xx = rx * basin1.rcpVector[0] + ry * basin1.rcpVector[1] + rz * basin1.rcpVector[2];
        yy = rx * basin1.rcpVector[3] + ry * basin1.rcpVector[4] + rz * basin1.rcpVector[7];
        zz = rx * basin1.rcpVector[6] + ry * basin1.rcpVector[5] + rz * basin1.rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        deltaPosX[i] = xx * basin1.latVector[0] + yy * basin1.latVector[1] + zz * basin1.latVector[2];
        deltaPosY[i] = xx * basin1.latVector[3] + yy * basin1.latVector[4] + zz * basin1.latVector[7];
        deltaPosZ[i] = xx * basin1.latVector[6] + yy * basin1.latVector[5] + zz * basin1.latVector[8]; 
    }

    double norm = 0.0;
    for (int i = 0; i < basin1.numberOfAtoms; i++)
    {
        norm += deltaPosX[i] * deltaPosX[i] + deltaPosY[i] * deltaPosY[i] + deltaPosZ[i] * deltaPosZ[i];
    }
    norm = sqrt(norm) * job.kmcParameters.basin2Delta;
    for (int i = 0; i < basin1.numberOfAtoms; i++)
    {
        if (saddle.frozen[i] == 0)
        {
            basin2.posX[i] = saddle.posX[i] + norm * deltaPosX[i];  //NB use of saddle positions
            basin2.posY[i] = saddle.posY[i] + norm * deltaPosY[i];
            basin2.posZ[i] = saddle.posZ[i] + norm * deltaPosZ[i];
        }
    }

    //relax to find basin2
    relStatus.setStatusFailed();
    
    rel.RelaxStructure(mpi, basin2, fld, basin2.numberOfAtoms, job.minParameters, basin2Energy, relStatus, outStream);

    
    //basin1.reportBasisDifference(mpi, basin2, spec, job.kmcBasinRadius, outStream);
    //check that it hasnt gone back into basin1
    newBasinFound = basin1.checkBasin(basin2, job.kmcParameters.kmcBasinRadius, job.kmcParameters.kmcMaxRadius); 
        
    if (newBasinFound)
    {
        relStatus.setStatusSuccess();
        if (rank == 0)
        {
            outStream << "\n relaxed into second basin " << std::endl;
            outStream << "\n basin difference for basin 2 : " << std::endl;
            
            basin1.reportBasisDifference(mpi, basin2, spec, job.kmcParameters.kmcBasinRadius, outStream);
        }
        //foundSaddle = true;
    }
    else
    {
        relStatus.setStatusFailed();
        newBasinFound = false;
    }

    return newBasinFound;
}
