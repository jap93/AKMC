/**
 * @file MPICommunicator.h
 * @author John Purton (john.purton@stfc.ac.uk)
 * @brief handles MPI communication within the ACDC program. Most functions are of "const" type to allow pass by reference
 * @version 0.1
 * @date 2021-06-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once 

#include <fstream>
#include <iomanip>

#include "Memory.h"
#include "Status.h"
#include "WorkGroup.h"

#ifdef NOMPI
    struct MPI_Comm
    {
        int comm;
    };
#else
    #include "mpi.h"
#endif
class MPIComms
{

private:

    #ifdef NOMPI
        MPI_Comm MPI_COMM_NULL, MPI_COMM_WORLD;
    #endif

    MPI_Comm 
        /** group communicator when MPI_COMM_WORLD is split */      
        groupComm = MPI_COMM_NULL;

    int
        /** the total number of processors/cores being used */
        numProcs = 1,

        /** the processor/core rank */
        rank = 0,

        /** list of the work group belonging to each core */
        *groupIdentity = nullptr,

        /** the group number for this core */
        idGrp = 0,

        /** the rank of a node within a group */
        grpRank = 0,

        /** the number of nodes in each workgroup */
        numWorkGrpProcs = 1,

        /** ther number of work groups */
        numWorkGroups = 1;

    
    std::vector<WorkGroup> 
        /** list of thread numbers in all workgroups */
        workGroups;

public:


    /**
     * @brief Construct a new MPIComms object
     * 
     */
    MPIComms();

    /**
     * @brief Destroy the MPIComms object
     * 
     */
    virtual ~MPIComms();
    /** 
     * @brief this should initialise the MPI communication but not functional yet
     */
    void initialise(void);

    /**
     * @brief call MPI_Abort with the world communicator
     * 
     */
    void commsAbortWorld(void) const;

    /**
     * @brief releases/frees the group communicator
     * 
     */
    void freeGroupComm(void);

    /**
     * @brief splits the  comm_world communicator to allow work groups
     * 
     * @param ptile (int) : 
     * @param offset (int&) : placeholder for the number of cores for the parent or zero'th workgroup that controls workgroups
     * @param mode (std::string) : placeholder for mode when different programs are used 
     * @param outStream (ofstream&) : stream for output file
     */
    void redistributeCommunicators(int ptile, int offset, std::string mode, std::ofstream& outStream);

    /**
     * @brief Get the Communicator object
     * 
     * @return MPI_Comm 
     */
    MPI_Comm getCommunicator(void) const    
        {return groupComm;}

    /**
     * @brief returns the processor rank on MPI_COMM_WORLD
     * 
     * @return rank (int) 
     */
    int getRank(void) const 
        {return rank;} 
    
    /**
     * @brief Get the number of processors on MPI_COMM_WORLD
     * 
     * @return numProcs (int) 
     */
    int getNumProcs(void) const
        {return numProcs;}

    /**
     * @brief Set the processor rank for MPI_COMM_WORLD
     * 
     * @param pid (int) : sets the rank of a processor
     */
    void setRank(int pid){rank = pid;}

    /**
     * @brief returns the rank of the thread within the group mpi communicator
     * 
     * @return groupRank (int) 
     */
    int getWorkGroupRank(void) const 
        {return grpRank;}

    /**
     * @brief returns the number of threads within the group mpi communicator
     * 
     * @return numWorkGrpProcs (int) 
     */
    int getNumWorkGroupProcs(void) const
        {return numWorkGrpProcs;} 

    /**
     * @brief Set the number of processors for MPI_COMM_WORLD
     * 
     * @param nump (int) : number of processors
     */
    void setNumProcs(int nump){numProcs = nump;}

    /**
     * @brief returns the group number for a given thread rank
     * 
     * @param p (int) : rank
     * @return workGroup (int)
     */
    int getGroupID(int p) const
        {return groupIdentity[p];}
    
    /**
     * @brief returns the group number
     * 
     * @return workGroup (int)
     */
    int getGroupID(void) const
        {return idGrp;}

    /**
     * @brief returns the number of work groups
     * 
     * @return numWorkGroups (int) 
     */
    int getNumOfWorkGroups(void) const
        {return numWorkGroups;}

# ifndef NOMPI

    /**
     * @brief sends a message to all threads that are not within the callers work group.
     * 
     * @param outStream (ofstream&) : stream for writing to a file
     * @return true 
     * @return false 
     */
    bool sendSignalToAll(std::ofstream& outStream) const;

    /**
     * @brief sends a message to a thread "i" that are not within the callers work group.
     * 
     * @param i (int) : the number of the thread
     * @param outStream (ofstream&) : stream for writing to a file
     * @return true 
     * @return false 
     */
    bool sendInterruptToSingleCore(int i, std::ofstream& outStream) const;

    void sendMessageToParent(int command) const;

    /**
     * @brief uses IProbe to check whether a message has been sent and then reads the message. Used to interrupt saddle search/and or relaxation
     * 
     * @param outStream  (ofstream&) : stream for writing to a file
     * @return true 
     * @return false 
     */
    void checkForInterrupt(Status& stat,std::ofstream& outStream) const;

    /**
     * @brief mpi barrier using the group communicator i.e. syncs a group
     * 
     */
    void barrierGroup(void) const;

    /**
     * @brief mpi barrier using the world communicator 
     * 
     */
    void barrierWorld(void) const;

    /**
     * @brief checks to see whether as been posted using MPI_IProbe 
     * 
     * @param command (int&) : the command sent for action
     * @param source (int&) : the source of the message
     * @param msgTag (int&) : the tag from the MPI message
     */
    void checkForMessage(int& command, int& source, int& msgTag, std::ofstream&) const;

    /**
     * @brief receives a command in the form of an integer from the parent
     * 
     * @return int : ineger value of the command
     */
    int waitForCommand(void) const;

    /**
     * @brief mpi barrier using the group communicator i.e. syncs a group
     * 
     * @param command (int) : the integer command to be sent
     * @param grp (int) : the group number
     * @param outStream (ofstream&) : stream for writing to a file
     */
    void sendCommandToGroup(int command, int grp, std::ofstream& outStream) const;

    /**
     * @brief not used
     * 
     * @return int 
     */
    int cancelRequest(void) const;

    /**
     * @brief broadcasts an integer to the work group
     * 
     * @param val (int) : the value being sent
     * @param source (int) : the source of the mpi call
     */
    void broadCastIntGroup(int& val, int source) const;

    /**
     * @brief broadcasts an integer to the world communicator
     * 
     * @param val (int) : the value being sent
     * @param source (int) : the source of the mpi call
     */
    void broadCastIntWorld(int& val, int source) const;

    /**
     * @brief broadcasts a long integer to the work group
     * 
     * @param val (long) : the value being sent
     * @param source (int) : the source of the mpi call
     */
    void broadCastLongGroup(long& val, int source) const;

    /**
     * @brief broadcasts a long integer to the world communicator
     * 
     * @param val (long) : the value being sent
     * @param source (int) : the source of the mpi call
     */
    void broadCastLongWorld(long& val, int source) const;

    /**
     * @brief sums a double value within a work group
     * 
     * @param dIn (double) : the input value to be summed
     * @return sum (double) 
     */
    double sumDoubleGroup(double dIn) const;

    /**
     * @brief sums a double value within world communicator
     * 
     * @param dIn (double) : the input value to be summed
     * @return sum (double) 
     */
    double sumDoubleWorld(double dIn) const;

    /**
     * @brief sums a double vector within a work group
     * 
     * @param dIn (double) : the input vector to be summed
     */
    void sumDoubleVectorGroup(double* dIn, int n) const;

    /**
     * @brief sums a double vector within world communicator
     * 
     * @param dIn (double) : the input vector to be summed
     * @param n (int) : the length of the vector
     */
    void sumDoubleVectorWorld(double* dIn, int n) const;

    /**
     * @brief sums an integer value within a work group
     * 
     * @param iIn (int) : the input value to be summed
     * @param n (int) : the length of the vector
     */
    int sumIntGroup(int iIn) const;

    /**
     * @brief sums an integer value within world communicator
     * 
     * @param iIn (int) : the input value to be summed
     * @return sum (int) 
     */    
    int sumIntWorld(int iIn) const;

    /**
     * @brief sums an integer vector within work group
     * 
     * @param iIn (int) : the input vector to be summed
     * @param n (int) : the length of the vector 
     */  
    void sumIntGroup(int* iIn, int n) const;

    /**
     * @brief sums an integer vector within world communicator
     * 
     * @param iIn (int) : the input vector to be summed
     * @param n (int) : the length of the vector 
     */
    void sumIntWorld(int* iIn, int n) const;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// dummy routines if no MPI
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# else
    void barrierGroup(void) const
        {return;}

    void barrierWorld(void) const
        {return;}

    double sumDoubleGroup(double dIn) const
        { return dIn;}

    double sumDoubleWorld(double dIn) const
        { return dIn;}

    void sumDoubleVectorGroup(double* dIn, int n) const
        {return;}

    void sumDoubleVectorWorld(double* dIn, int n) const
        {return;}
    
    int sumIntGroup(int iIn) const
        {return iIn;}
  
    int sumIntWorld(int iIn) const
        {return iIn;}
  
    void sumIntGroup(int* iIn, int n) const
        {return;}

    void sumIntWorld(int* iIn, int n) const
        {return;}

#endif

};



