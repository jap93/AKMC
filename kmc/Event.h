#pragma once 


#include <fstream>
#include <vector>
#include <iomanip>
#include <iterator>
#include "mpi.h"


#include "../Common/MPICommunicator.h"
#include "../Common/Memory.h"
#include "../Common/Basis.h"

class Event
{

private:

    float
        /** vector containing positions of saddle*/
        *sadPosX = nullptr,
        *sadPosY = nullptr,
        *sadPosZ = nullptr,
        /** vector containing positions of basis*/
        *bas2PosX = nullptr,
        *bas2PosY = nullptr,
        *bas2PosZ = nullptr;
        
    int
        /** the total number of atoms sent back - may not be the same as the no in the basis class*/
        numberOfAtomsSent;

    int
        /** the global number of the atom */
        *globNumber = nullptr;

    double
        /** basin 1 energy */
        basin1Energy = 0.0,

        /** basin 2 energy */
        basin2Energy = 0.0,

        /** saddle energy */
        saddleEnergy = 0.0,
         
        /** activation energy for the event */
        actEnergy = 0.0,

        /** the rate of the processes */
        rate = 0.0;
        
    double
        /** the centre for "cluster" searches */
        *centre = nullptr;


public:

    /** 
     * @brief returns the calculated rate 
     * @return rate (double)
     */
    double getRate(void)
    { return rate;}

    /** 
     * @brief returns the activation energy 
     *
     * @return actEnergy (double)
     */
    double getActivationEnergy(void)
    { return actEnergy;}

    /** 
     * @brief returns the basin 1 energy 
     * 
     * @return basin1Energy (double)
     */
    double getBasin1Energy(void)
    { return basin1Energy;}

    /** 
     * @brief returns the basin 2 energy 
     * 
     * @return basin2Energy (double)
     */
    double getBasin2Energy(void)
    { return basin2Energy;}

    /** 
     * @brief returns the activation energy 
     * 
     * @return saddleEnergy (double)
     */
    double getSaddleEnergy(void)
    { return saddleEnergy;}

    double* getCentre(void)
    { return centre;}

    /**
     * @brief write to file the structure at the saddlepoint
     * 
     * @param spec (const Species&) : element information
     * @param atmLabel (int*) : integer label of atom type
     * @param outStream (std::ofstream&) : the output stream
     */
    void writeSaddle(const Species& spec, int* atmLabel, std::ofstream& outStream);


    /**
     * @brief write to file the structure of the second basin
     * 
     * @param spec (const Species&) : element information
     * @param atmLabel (int*) : integer label of atom type
     * @param outStream (std::ofstream&) : the output stream
     */
    void writeBasin2(const Species& spec, int* atmLabel, std::ofstream& outStream);

    /**
     * @brief write to file the centre point of the search
     * 
     * @param outStream (std::ofstream&) : the output stream
     */
    void writeCentre(std::ofstream& outStream);

    /** 
     * @brief copies the positions of an event to those of a basin and the saddle point
     * 
     * @param bas (Basis&) : the basin 
     * @param outStream (std::ofstream&) : the output stream
     */
    void setPositions(Basis& bas, std::ofstream& outStream);

    /**
     * @brief Get the Displaced Atoms object
     * 
     * @param basin (const Basis&) : the atom positions in the original basin for comparison
     * @param basinRadius (double) : the minimum distance required between atoms to consider new basin
     * @return std::vector<int> : list of atoms that have been displaced more than basinRadius
     */
    std::vector<int> getDisplacedAtoms(const Basis& basin, double basinRadius);

    /**
    * constructor - makes sure all pointers are set to nullptr and energies are zero
    */
    Event();

    /**
     * @brief Construct a new Event object
     * 
     * @param bas1 (Basis&) : the basis 1 positions etc
     * @param bas2 (Basis&) : the basis 2 positions etc 
     * @param sad  (Basis&) : the saddle point positions etc
     * @param eact (double) : activation energy
     * @param prefactor (double) : pre-factor of the adaptive kinetic MC
     * @param temperature (double) : temperature
     * @param bas1Energy (double) : basin 1 energy
     * @param bas2Energy (double) : basin 2 energy
     * @param sadEnergy (double) : energy of the saddle point energy config
     */
    //Event(Basis& bas1, Basis& bas2, Basis& sad, double eact, double prefactor, double temperature, double bas1Energy, double bas2Energy, double sadEnergy);

    /**
     * @brief Construct a new Event object
     * 
     * @param veci (int*) : contains the global no passed over mpi
     * @param vecf (float*) : contains the configuration information (often passed over mpi)
     * @param natoms (int) : the number of atoms
     * @param prefactor (double) : pre-factor of the adaptive kinetic MC
     * @param temperature (double) : temperature
     * 
     */
    Event(int* veci, float* vecf, int natoms, double prefactor, double temperature);

    /**
    * @brief Copy constructor. Allocates memory to arrays and performs a deep copy.
    * 
    * @param oldEvent (const Event&) : is the positions etc of the old/original object/configuration.
    */
    Event(const Event& oldEvent);

    /**
    * @brief Destructor. Deallocates dynamically allocated arrays and sets pointers back to nullptr
    */
    virtual ~Event();

    /**
    * @brief overloaded assignment operator. Performs a deep copy
    * 
    * @param src (const Event&) : is the atomic positions and rate being copied
    * @return returns a copy of the Event positions, atom identities, energies
    */
    Event& operator=(const Event& src);

};


