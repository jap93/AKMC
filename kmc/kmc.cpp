/*
The  KMC module controls the calling of parasite programs that actually determine the activateion energy/rate.
I dont think it needs to calculate energies - although the field is read to get species.
*/

#include <random>
#include <chrono>
#include <time.h>
#include <iostream>
#include <fstream>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "Basis.h"
#include "JobControl.h"
#include "Species.h"
#include "Field.h"
#include "MPICommunicator.h"


#include "KineticMonteCarlo.h"

//random number generator

std::mt19937_64 generator;
std::uniform_real_distribution<double> dist{ 0.0, 1.0 };

int main(int argc, char* argv[])
{
    int
        rank = 0,
        numProcs = 1,
        restartIteration = 0;

    double
        restartTime = 0.0;

    JobControl
        job;

    Species
        spec;

    Basis
        bas;

    Field
        *fld = nullptr;

    NbrListPBC
        nbrList;

    KineticMonteCarlo
        kMC;

    MPIComms
        mpi;

    std::string
        basisFileName = "basis",
        fieldFileName = "potentials",
        logFileName = "kmc.log",
        jobFileName = "control";

    std::ofstream
        outStream;

    std::ifstream
        inStream;


    //
    //  Initialize MPI.
    //
    #ifdef _OPENMP
        int mpiInitLevel = 0;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpiInitLevel);
    #else
        MPI_Init(&argc, &argv);
    #endif
 
    //MPI_Comm groupComm = mui::mpi_split_by_app();
    //Get the processor number
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //  Get the number of processors.
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int err = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    mpi.setRank(rank);
    mpi.setNumProcs(numProcs);

    //open rank 0 only
    if (rank == 0)
    {
        //logFileName = logFileName + std::to_string(rank);
        outStream.open("kmc.out", std::ofstream::out | std::ofstream::app);
    }

    #ifdef _OPENMP
        //int MPI_Init_thread(  int * argc, char ** argv[],int thread_level_required, int * thead_level_provided);
        //int MPI_Query_thread( int * thread_level_provided);
        //int MPI_Is_main_thread(int * flag);
        int numThread = 1;
        int tid = 0;
        #pragma omp parallel firstprivate(tid, numThread)
        {
            numThread = omp_get_num_threads();
            tid = omp_get_thread_num();
            outStream << "\n mpi and openmp " << rank  << " " << numProcs << " " << tid << " " << numThread << std::endl;
            outStream.flush();
        }

        if (mpiInitLevel != MPI_THREAD_FUNNELED)
        {
            outStream << "\n MPI_THREAD FUNNELED is required for hybrid parallelisation"  << std::endl;
            outStream.flush();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }
        //make sure control files are correct - efectively a single core per worker
        job.ptile = 1;
        
    #endif

    

    //read the job control parameters and then write out
    inStream.open(jobFileName, std::ios::in);

    if (inStream.fail())              // check to see file is there
    {
        outStream << "\n*** could not find control file" << std::endl;
        outStream.flush();
        outStream.close();
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    job.readJobControl(inStream, outStream);
    inStream.close();
    
    if (rank == 0)
        job.writeSimulationDetails(mpi, outStream);

    //set parallel strategy
    mpi.redistributeCommunicators(job.numGrps, job.parentCores, "kmc", outStream);

    // open all other files (if required)
    if (rank != 0)  //open up files for other workgroups
    {	
        //open log file - use append and have each core produce a file. Also add interface name so we can track the files
        logFileName = logFileName + std::to_string(rank);
        outStream.open(logFileName, std::ofstream::out);
    }

    //set kmc parameters
    job.kmcParameters.kmcMaxCap = job.kmcParameters.kmcMinCap + (job.kmcParameters.kmcWindow * job.kmcParameters.kmcTemperature * BOLTZMANN);
    
    if (rank == 0)
        outStream << "\n maximum activation energy : " <<  job.kmcParameters.kmcMaxCap << std::endl;


    //create the random number generator - the seed has to be identical within all members of the group
    //and you cannot rely on the clock being the same!
    int seed1 = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    if (seed1 < 0)
        seed1 *= -1;
	
    //use the group id to make sure the different work groups have different trajectories
    seed1 += rank;
    generator.seed(seed1);

    //only the species present are required - just for data processing
    inStream.open(fieldFileName, std::ios::in);

    if (inStream.fail())              // check to see file is there
    {
            outStream << "\n*** could not find potentials file" << std::endl;
            outStream.flush();
            outStream.close();
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }
    
    fld = new Field;
    fld->readPotential(inStream, outStream, spec);
    inStream.close();
    

    //read data from a file
    
        inStream.open(basisFileName, std::ios::in);

        if (inStream.fail())              // check to see file is there
        {
            outStream << "\n*** could not find basis file" << std::endl;
            outStream.flush();
            outStream.close();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }

        //what about lattice parameters
        bas.inputBasis(restartIteration, restartTime, spec, job, inStream, outStream);
        bas.resetSimulationCell();
        bas.dumpBasis(mpi, spec, 0.0, 0.0, 0, outStream);
    
        inStream.close();
    if (rank == 0)
    {
        bas.printVectors(outStream);
        bas.printRcpVec(outStream);
        bas.printBasis(spec, outStream);
        
    }

    //setup field
    fld->setup(mpi, bas.latVector, bas.rcpVector, bas.minDimension(), bas.cellSize(), bas.numberOfAtoms, spec, outStream);
    
    auto startTime = MPI_Wtime();

    //call the KMC class that handles all the calculation
    kMC.runKMC(mpi, restartIteration, restartTime, spec, fld, job, bas, outStream);

    auto finishTime = MPI_Wtime();
    auto diff = finishTime - startTime;
    if (rank == 0)
        outStream << std::scientific << std::setprecision(7) << "\n total time for kmc : " << diff << " seconds" << std::endl;

    //shut down and destroy all memory

    if (rank == 0)
    {
        outStream.flush();
        outStream.close();
    }

    delete fld;

    MPI_Finalize();
    return EXIT_SUCCESS;
}

/*
Rerturn a random number
*/
double getRandomNumber(void)
{
    return dist(generator);
}

