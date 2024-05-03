# pragma once 

/**********************************************************
* Interface to the Activation Relaxation Technique (ART) class
* The Original method was developed by Malek and Mousseau Phys
* Rev E vol 62 pages 7723-7728 (not this has been expanded/improved
* on by subsequent papers of this group). However, the implementation
* here follows more closely that of Olsen et al J. Chem. Phys. vol 121
* pages 9776-9792. The variable nomenclature follows their paper as far as
* possible, but I have attempted to give more descriptive names in addition!
*
* updated June 2020 in Lanczos to get rid of Eigen and the way initial iteration in Lanczos is set up (ie realised I could do it in the loop)
*/

#include "SaddleSearch.h"

#include "Memory.h"
#include "Basis.h"
#include "JobControl.h"



enum artStatus 
{
    LANCZOSCONVPOSEIGENVALUE = 3, LANCZOSCONVNEGEIGENVALUE = 4, LANCZOSFAILED = 5, LANZOSINTERRUPTED = 6
};

class ART : public SaddleSearch
{
private:

    
    int
        /* the iteration number for the art relaxation */
        itART = 0,
        /* number of the Lanczos vectors */
        numLanczosVectors = 0;

    double
	ARTCentre[3];

    bool
        /* flag to indicate restart */
        ARTrestart = false;

    /**
     * @brief checks that the number of vectors in the Lanczos minimisation is sufficient
     * 
     */
    void checkVectorSize(void);

    /** 
     * @brief calculates the forces on atoms and returns a single vector for the lin algebra rather than the three used in MD etc.
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param forceVector (double*) : the single vector that is returned
     * @param bas (Basis&) : atomic information
     * @param fld (Field*) : interaction parameter
     * @param outStream (ofstream&) : output stream for the workgroup
     * @param totalEnergy (Energy&) : the calculated total energy
     */
    void calculateForces(const MPIComms& mpi, double* forceVector, Basis& bas, Field* fld, Energy& totalEnergy, std::ofstream& outStream);




    /**
     * @brief minimises to the saddle point using the FIRE method
     * 
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param bas (Basis&) : atomic information
     * @param spec (const Species&) : used to dump trajectory
     * @param centre (double*) : center point required for restart
     * @param finalenergy (Energy&) : the energy of the saddlepoint
     * @param fld (Field*) : the method/style of force calculation
     * @param artParameters (const ARTCntrl&) : run parameters
     * @param status (Status&) : validity of the calculation
     * @param goPerp (bool) : use perpendicular forces only
     * @param outStream (ofstream&) : output stream for the workgroup
     */
    void FIRESaddlePoint(const MPIComms& mpi, Basis& bas, const Species& spec, double* centre, Energy& finalenergy,  Field* fld, const ARTCntrl& artParameters, 
                         Status& status, bool goPerp, std::ofstream& outStream);

    /**
     * @brief calculates the forces using Lanczos method
     * 
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param artParameters (const ARTCntrl&) : run parameters
     * @param bas (Basis&) : atomic information
     * @param eigenValue (double&) : the lowest eigenvalue
     * @param forcePar (double*) : forces parallel to the basin
     * @param forcePerp (double*) : forces perpendicular to the basin
     * @param fld (Field*) : the method/style of force calculation
     * @param totalEnergy (Energy&) : the energy of the cell on this iteration
     * @param outStream (ofstream&) : output stream for the workgroup
     * @return converged (bool) : whether the Lanczos has converged
     */
    bool calculateLanczosForce(const MPIComms& mpi, const ARTCntrl& artParameters, Basis& bas, double& eigenValue, double* forcePar, 
                               double* forcePerp, Field* fld, Energy& totalEnergy, std::ofstream& outStream);

    /**
     * @brief the interface to lapack (dstev) for the diagonalisation of the tridiagonal matrix
     * 
     * @param diagonalT (double*) : diagonal element of the matrix
     * @param subDiagT (double*) : subdiagonal elements of the matrix
     * @param eigVal (double&) : the lowest eigenvalue
     * @param eigVec (double*) : the eigenvector of the lowest eigenvalue
     * @param n (double) : the size of the matrix
     */
    void diagonalise(double* diagonalT, double* subDiagT, double& eigVal, double* eigVec, int n);


    void saveARTDetails(const MPIComms& mpi, int* frozen, double timeStep, int natom, std::ofstream& outStream);

    
    void loadARTDetails(int* frozen, double& timeStep, int natom, std::ifstream& inStream);


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

    ART();
    virtual ~ART();
};
