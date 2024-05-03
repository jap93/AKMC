#pragma once
#include<unistd.h>
#include <chrono>
#include <thread>

#include "Basis.h"
#include "JobControl.h"
#include "Species.h"
#include "Field.h"

#include "MPICommunicator.h"

#include "ART.h"
#include "Dimer.h"

#include "Event.h"

using namespace std::chrono;
using namespace std::this_thread;     // sleep_for, sleep_until
using namespace std::chrono_literals;

class KineticMonteCarlo
{
private:

    void printData(double kineticEnergy, double rcpEnergy, double realEnergy, double vdwEnergy, int natoms, std::ofstream& outStream);

    /**
     * @brief function that selects region to search over and coordinates activity of workers
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param restartIteration (int&) : iteration number on restart
     * @param restartTime (double&) : kmc time on restart
     * @param spec (const Species&) : species information
     * @param fld (Field*) : force-field parameters etc for force calculation
     * @param job (const JobControl&) : run parameters
     * @param bas (Basis&) : atomic positions
     * @param outStream (std::ofstream&) : output streams
     */
    void driver(const MPIComms& mpi, int& restartIteration, double& restartTime, const Species& spec, Field* fld, const JobControl& job, Basis& bas, 
                std::ofstream& outStream);

    /**
     * @brief function that undertakes saddle searches
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param spec (const Species&) : species information
     * @param fld (Field*) : force-field parameters etc for force calculation
     * @param job (const JobControl&) : run parameters
     * @param bas (Basis&) : atomic positions
     * @param outStream (std::ofstream&) : output streams
     */
    void worker(const MPIComms& mpi, const Species& spec, Field* fld, const JobControl& job, std::ofstream& outStream);

    /**
     * @brief calculates the events of a sinle KMC iteration - multithread version
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param spec (const Species&) : information on atomic species
     * @param fld (Field*) : class for energy/force evaluation
     * @param job (onst JobControl&) : run parameters
     * @param basin1 (Basis&) : the initial configuration
     * @param basin1Energy (double&) : energy of the initial basis
     * @param eventList (std::vector<Event>&) : list of all the events
     * @param excludeList (std::vector<int>&) : allows an event to be excluded is it has moved previously
     * @param iteration (int) : the kmc iteration number (used for the mpi tag)
     * @param outStream (std::ofstream&) : the output stream
     */
    
    void calculateRateConstant(const MPIComms& mpi, const Species& spec, Field* fld, const JobControl& job, Basis& basin1, 
                               double& basin1Energy, std::vector<Event>& eventList, const std::vector<int>& excludeList, int iteration, 
                               std::ofstream& outStream);

    /**
     * @brief calculates the events of a sinle KMC iteration - multithread version but experimental so do not use it unless you want to code
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param spec (const Species&) : information on atomic species
     * @param fld (Field*) : class for energy/force evaluation
     * @param job (onst JobControl&) : run parameters
     * @param basin1 (Basis&) : the initial configuration
     * @param basin1Energy (double&) : energy of the initial basis
     * @param eventList (std::vector<Event>&) : list of all the events
     * @param excludeList (std::vector<int>&) : allows an event to be excluded is it has moved previously
     * @param iteration (int) : the kmc iteration number (used for the mpi tag)
     * @param outStream (std::ofstream&) : the output stream
     */
    
    void calculateRateConstantExperimental(const MPIComms& mpi, const Species& spec, Field* fld, const JobControl& job, Basis& basin1, 
                               double& basin1Energy, std::vector<Event>& eventList, const std::vector<int>& excludeList, int iteration, 
                               std::ofstream& outStream);
    /**
     * @brief calculates the events of a sinle KMC iteration - single thread version (mostly for testing functionality in other parts of the code)
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param spec (const Species&) : information on atomic species
     * @param fld (Field*) : class for energy/force evaluation
     * @param job (onst JobControl&) : run parameters
     * @param basin1 (Basis&) : the initial configuration
     * @param basin1Energy (double&) : energy of the initial basis
     * @param eventList (std::vector<Event>&) : list of all the events
     * @param iteration (int) : the kmc iteration number (used for the mpi tag)
     * @param outStream (std::ofstream&) : the output stream
     */
    void calculateRateConstantSingleCore(const MPIComms& mpi, const Species& spec, Field* fld, const JobControl& job, Basis& basin1, 
                               Energy& basin1Energy, std::vector<Event>& eventList, int iteration, std::ofstream& outStream);

    /**
     * @brief calculates the KMC time increment
     * 
     * @param totalRate (double) : the sum of all the event rates
     * @return kmcTime (double) 
     */
    double calculateTime(double totalRate);

    /**
     * @brief resets the list of saddle points and conducts recucling of saddle points if required
     * 
     * @param eventList (std::vector<Event>&) : vectors containing the events
     * @param chosenEvent (int) : the chosen event whose centre will be used to check others
     * @param radius (double) : the minimum distance allowed between centres
     * @param recycle (bool) : flag to indicate whether all or recycling is employed
     * @param global (bool) : flag for global relaxation (if true then list is cleared)
     * @param mode (int) : the style of recycling
     * @param percent (double) : if a random cull is also done then this is the percentage to remove
     * @param basin (Basis&) : the original basin 
     * @param excludedList (std::vector<int>&) : the excluded list of atoms that can move
     * @param basinRadius (double) : the minimum distance that must be moved
     * @param superBasin (bool) : flag to indicate super basin method being used
     * @param outStream (std::ofstream&) : output stream
     */
    void clearEventList(std::vector<Event>& eventList, double* latVector, double* rcpVector, int chosenEvent, double radius, bool recycle, 
                        bool global, int mode, double percent, Basis& basin, std::vector<int>& excludedList, double basinRadius, bool superBasin, std::ofstream& outStream);

    /**
     * @brief determines whether the atom displaced is allowed under superbasin rules
     * 
     * @param basin (Basis&) : the initial configuration for camparison with the event
     * @param ev (Event&) : the Event for validation
     * @param excludedList (const std::vector<int>&) : the list of atoms not allowed to move
     * @param iteration (int) : the kmc iteration number
     * @param basinRadius (double) : the minimum distance that must be moved
     * @param outStream  (std::ofstream&) : output stream
     * @return true 
     * @return false 
     */
    bool checkEventValidity(Basis& basin, Event& ev, const std::vector<int>& excludedList, int iteration, double basinRadius, std::ofstream& outStream);

    /**
     * @brief determines whether the atom displaced is allowed under superbasin rules
     * 
     * @param basin (Basis&) : the initial configuration for camparison with the event
     * @param ev (Event&) : the chosen event
     * @param excludedList (const std::vector<int>&) : the list of atoms not allowed to move
     * @param basinRadius (double) : the distance required for an atom to move to be considered in a new basin
     * @param history (int) : the number of iterations that the atom is not allowed to move
     * @param basinRadius (double) : the minimum distance that must be moved
     * @param outStream  (std::ofstream&) : output stream
     */
    void updateExcludedList(Basis& basin, Event& ev, std::vector<int>& excludedList, double basinRadius, int history, std::ofstream& outStream);

    /**
     * @brief the zero'th workgroup sends basin 1 positions to other workgroups
     * 
     * @param mpi  (const MPIComms&) : MPI workgroup and parallel information
     * @param basin1 (Basis&) : the initial configuration
     * @param outStream (std::ofstream&) : the output stream
     */
    void broadCastPositions(const MPIComms& mpi, Basis& basin1, std::ofstream& outStream);

    void sendErrorCode(const MPIComms& mpi, int code, int tag, std::ofstream& outStream);

    /**
     * @brief sends the positions of basins and saddlepoint to rank 0
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param basin1 (Basis&) : the initial configuration
     * @param basin2 (Basis&) : the basin 2 configuration
     * @param saddle (Basis&) : the saddle point configuration
     * @param centre (double*) : search centre point
     * @param basin1Energy (double&) : energy of the initial basin
     * @param basin2Energy (double&) : energy of basin 2
     * @param saddleEnergy (double&) : energy of saddle configuration
     * @param actEnergy (double&) : activation energy
     * @param destination (int) :destination rank
     * @param outStream (std::ofstream&) : the output stream
     * 
     * @return true/false whether an interrupt signal was detected.
     */
    bool sendEvent(const MPIComms& mpi, Basis& basin1, Basis& basin2, Basis& saddle, double* centre, double basin1Energy, double basin2Energy, 
                                      double saddleEnergy, double actEnergy, int destination, std::ofstream& outStream);

    /**
     * @brief function to receive an event (on rank 0)
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param ev (Event&) : the event information
     * @param preFactor (double) : kmc prefactor
     * @param temperature (double) : simulation temperature
     * @param size (int) : length of vector to receive
     * @param source (int) : the sorce of the mpi message
     * @param tag (int) : the MPI tag
     * @param outStream  (std::ofstream&) : the output stream
     */
    void receiveEvent(const MPIComms& mpi, Event& ev, double preFactor, double temperature, int& size, 
                                    int& source, int& tag, std::ofstream& outStream);

    /**
     * @brief function to receive an event (on rank 0)
     *
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param ev (Event&) : the event information
     * @param numberOfAtoms (int) : number of atoms
     * @param preFactor (double) : kmc prefactor
     * @param temperature (double) : simulation temperature
     * @param boltz (double) : boltzmann constant
     * @param source (int) : the source of the expected positions
     * @param msgTag (int) : the tag of the expected message
     * @param outStream  (std::ofstream&) : the output stream
     */
    void receivePositions(const MPIComms& mpi, Event& ev, int numberOfAtoms, double preFactor, double temperature, double boltz, int source, int msgTag,
                         std::ofstream& outStream);
                         
    /**
     * @brief writes out the complete eventlist to a file events.#
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param spec (const Species&) : species information
     * @param basin1 (Basis&) : the initial configuration
     * @param eventList (std::vector<Event>&) : list of all the events
     * @param simulationTime (double) : clock time
     * @param basin1Energy (double) : basin 1 energy
     * @param chosenEvent (int) : the event number that was selected by kmc routine
     * @param iteration (int) : kmc iteration
     */
    void writeEventList(const MPIComms& mpi, const Species& spec, Basis& basin1, std::vector<Event>& eventList, double simulationTime, 
                        double basin1Energy, int chosenEvent, int iteration);
    
    /**
     * @brief writes out the excluded list for super basins
     *
     * @param excludedList (std::Vector<int>) : the excluded list
     * @param outStream (std::ofstrem&) : the output stream fro the excluded list
     */
    void writeExcludedList(std::vector<int> excludedList, std::ofstream& outStream);

    bool checkForNewBasin(const MPIComms& mpi, const Species& spec, const JobControl& job, Field* fld, Basis& basin1, Basis& basin2, 
                                         Basis& saddle, Status& relStatus, Energy& basin2Energy, std::ofstream& outStream);

public:

    /**
     * @brief public function to undertake the simulation
     * 
     * @param mpi (const MPIComms&) : MPI workgroup and parallel information
     * @param restartIteration (int&) : iteration number on restart
     * @param restartTime (double&) : kmc time on restart
     * @param spec (const Species&) : species information
     * @param fld (Field*) : force-field parameters etc for force calculation
     * @param job (const JobControl&) : run parameters
     * @param bas (Basis&) : atomic positions
     * @param outStream (std::ofstream&) : output streams
     */
    void runKMC(const MPIComms& mpi, int& restartIteration, double& restartTime, const Species& spec, Field* fld, const JobControl& job, Basis& bas, 
                std::ofstream& outStream);

    KineticMonteCarlo();
    virtual ~KineticMonteCarlo();
};
