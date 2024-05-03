# pragma once

/**********************************************************
* Interface to the Dimer method class
* The Original method was developed by Henkelman and Jónsson J. Chem.
* Phys. vol 111 pages 7010-7022 (1999). However, the implementation
* here follows more closely that of Olsen et al J. Chem. Phys. vol 121
* pages 9776-9792, using fewer force calculations and extrapolating
* forces for second dimer point.
*
* originally coded in by Michael Seatom and tidied up by John Purton
*/

#include "SaddleSearch.h"

#include "Memory.h"
#include "Basis.h"
#include "JobControl.h"

enum dimerStatus
{
    DIMERFAILED = 5, DIMERINTERRUPTED = 6
};

class Dimer : public SaddleSearch
{
private:

    
    int
        /* the status of the saddle search */
        dimerStatus = DIMERFAILED;


    /** calculates the forces on atoms and returns a single vector for total force acting on all atoms in dimer
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param forceVector (double*) : the single vector (total dimer force) that is returned
     * @param bas (Basis&) : atomic information
     * @param fld (Field*) : interaction parameter
     * @param totalEnergy (double&) : the calculated total energy
     */
    void calculateForces(const MPIComms& mpi, double* forceVector, Basis& bas, Field* fld, double& totalEnergy);



    
    void FIRESaddlePoint(const MPIComms& mpi, Basis& bas, double* rVector, double& finalenergy, Field* fld, const DimerCntrl& dimerParameters, 
                         Status& status, bool goPerpendicular, std::ofstream& outStream);

public:

        /** 
     * @brief findSaddlePoint is the main routine to determine the activation energy and saddle point of a configuration. Some child classes also have publicly exposed classes
     * to allow for additional functionality and/or testing.
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param spec (const Species&) : element information
     * @param job (const JobControl&) : run parameters controlling execution
     * @param fld (Field*) : the method/style of force calculation
     * @param basin1 (Basis&) : atomic information for the first basin
     * @param saddle (Basis&) : atomic information for the saddle point
     * @param status (Status&) : validity of the calculation
     * @param basin1Energy (Energy&) : the calculated total energy of the first basin
     * @param saddleEnergy (Energy&) : the calculated total energy of the saddle point
     * @param centre (double*) : the midpoint of the region used to calculate saddle points
     * @param continueSearch (bool) : indicates that search should be resumed not started
     * @param outStream (ofstream&) : output stream for the workgroup
     */
    virtual void findSaddlePoint(const MPIComms& mpi, const Species& spec, const JobControl& job, Field* fld, Basis& basin1, Basis& saddle,
                                 Status& status, Energy& basin1Energy, Energy& saddleEnergy, double* centre, bool continueSearch,
                                 std::ofstream& outStream);

    Dimer();
    virtual ~Dimer();
};
