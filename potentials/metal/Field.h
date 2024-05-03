#pragma once

#include <string>
#include <stdio.h>
#include <algorithm>

#include "Species.h"
#include "Memory.h"
#include "Constants.h"
#include "MPICommunicator.h"
#include "Energy.h"

//#include "Ewald.cuh"
#include "NbrListPBC.h"
#include "ManyBody.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

class Field
{
private:

    int
        
        // the number of cells to be searched over when not using the image convention
        cellX = 1,
        cellY = 1,
        cellZ = 1,

        /** the size of the density work arrays */
        wrkSize = 0,

        /** the maximum number of neighbours */
        maximumNbrs = 0;


    bool
        /** the neighbour list is used if true - faster but can have problems in ART*/
        useNbrList = false,

        /** flag to indicate that the nbr list requires updating */
        recalcNbrList = false,

        /** turns of image convention if true */
        minimumImage = true;

    double

		/** limit of extent of two-body interaction */
	    shortRangeCut = 15.0,

        
		/** the extent of the shell around the cutoff */
		verletShell = -0.5,
        
        /** conversion unit for all energies */
        energyUnit = 1.0;

    double
        /** temp storage of atom density */
        *tmpRho = nullptr,

        /** work storage of rho */
        *wrkRho = nullptr;


    
    NbrListPBC
        nbrList; 
        
    

    ManyBody
        mbdy;
  

    void sumForces(const MPIComms& mpi, double*, double*, double*, int);

    //void resetSimulationCell(double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int numatoms);

    void resetAtomSimulationCell(int i, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int numAtoms);

    void realSpaceForceImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                             double* forceX, double* forceY, double* forceZ, double* stress, double* latVector,
                             double* rcpVector, double* atomcharge, int* intlabel, int numatoms);

    void realSpaceForceImageNbrList(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                             double* forceX, double* forceY, double* forceZ, double* stress, double* latVector,
                             double* rcpVector, double* atomcharge, int* intlabel, int numatoms);

    void realSpaceForceNoImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                             double* forceX, double* forceY, double* forceZ, double* stress, double* latVector,
                             double* rcpVector, double* atomcharge, int* intlabel, int numatoms);

    void manyBodyForceImage(const MPIComms& mpi, double& mbdyEnergy, double* posX, double* posY, double* posZ, 
                                double* forceX, double* forceY, double* forceZ, double* latVector,
                                double* rcpVector, double* atomcharge, int* intlabel, int numatoms);



    void realSpaceEnergyImage(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                              double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int numatoms);

    void realSpaceEnergyImageNbrList(const MPIComms& mpi, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                     double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int numatoms);

    //void manyBodyEnergyImage(double& mbdyEnergy, double* posX, double* posY, double* posZ, 
    //                        double* latVector, double* rcpVector, double* atomcharge, int* intlabel, NbrListPBC& nbrList, int numatoms);

    void realSpaceAtomEnergyImage(const MPIComms& mpi, int aatom, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                             double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int numatoms);

    void realSpaceAtomEnergyImageNbrList(const MPIComms& mpi, int aatom, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                             double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int numatoms);

    void realSpaceAtomEnergyRemoveImage(const MPIComms& mpi, int aatom, double& manyEnergy, double& twoenergy, double* posX, double* posY, double* posZ, 
                                    double* latVector, double* rcpVector, double* atomCharge, int* atmLabel, int numAtoms);

    double calculateAtomRho(int aatom, double* posX, double* posY, double* posZ, 
                                 double* latVector, double* rcpVector, int* atmLabel, int numAtoms);

    double calculateAtomRhoNoImage(int aatom, double* posX, double* posY, double* posZ, 
                                 double* latVector, double* rcpVector, int* atmLabel, int numAtoms);
public:

    // constructor
    Field();

    // default destructor
    virtual ~Field();

    /**
     * return the type of periodic boundary as bool: True for minimum image convention
    */
    bool getImageConvention()
        {return minimumImage;}
        
    void resetSimulationCellNoImage(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int natoms);

    void resetSimulationCellMinImage(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int natoms);


    /**
     * @brief Set the External Communicator within a library
     * 
     * @param mpi (const MPIComms&) : the mpi communicatior
     * @param outStream (std::ofstream&) : the output stream
     */
    void setExternalCommunicator(const MPIComms& mpi, std::ofstream& outStream);

    /**
     * sets up the environment to underatke a calculation using the desired calculator
     * @param mpi (const MPIComms&) : the mpi communicatior
     * @param outStream (ofstream&) : the output stream
     */
    void createFileSystem(const MPIComms& mpi, std::string rootName, std::string externalFile, bool restart, std::ofstream& outStream);

    /**
     * sets up the environment to underatke a calculation using the desired calculator
     * @param mpi (const MPIComms&) : the mpi communicatior
     * @param latVector (double*) : lattice vectors
     * @param rcpVector (double*) : reciprocal lattice vectors
     * @param minDimension (double) : the smallest dimension of the cell
     * @param volume (double) | the cell volume
     * @param spec (Species&) : the elemental data
     * @param outStream (ofstream&) : the output stream
     */
    void setup(const MPIComms& mpi, double* latVector, double* rcpVector, double minDimension, double volume, int numAtoms, const Species& spec, std::ofstream& outstream);

    /**
     * tidies up after a calculation for example may move restart files etc
     * @param mpi (const MPIComms&) : the mpi communicatior
     */
    void tidy(const MPIComms& mpi);

    /**
     * copies any work files required
     * @param mpi (const MPIComms&) : the mpi communicatior
     */
    void copyWorkFiles(const MPIComms& mpi);

    /**
     * reads any information reuired to calculate the energy. If an incorrect keyword is used an invalid argument is thrown.
     * @param inStream (ifstream&) : the input stream
     * @param outStream (ofstream&) : the output stream
     * @param spec (Species&) : the elemental data
     */
    void readPotential(std::ifstream& inStream, std::ofstream& outStream, Species& spec);

    void calculateForces(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ, 
                         double* stress, double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int* frozen, int numatoms, Energy& eng);

    void calculateConfig(const MPIComms& mpi, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, double* atomcharge, 
                         int* intlabel, int* frozen, int numatoms, Energy& eng);

    void calculateEnergy(const MPIComms& mpi, double* posX, double* posY, double* posZ, 
                        double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int* frozen, int numatoms, 
                        Energy& eng);

    void calculateAtomEnergy(const MPIComms& mpi, int aatom, double* posX, double* posY, double* posZ, 
                        double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int* frozen, int numatoms, 
                        Energy& eng);

    void calculateAtomEnergyRemove(const MPIComms& mpi, int aatom, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, 
                                double* atomCharge, int* atmLabel, int* frozen, int numAtoms, Energy& eng);

    void calculateAtomEnergyDiff(const MPIComms& mpi, int aatom, double* oldPos, double* newPos, double* posX, double* posY, double* posZ, 
                            double* latVector, double* rcpVector, double* atomcharge, int* intlabel, int* frozen, int numatoms, 
                            Energy& eng);

    void print(const MPIComms& mpi, std::ofstream& outstream);

    void updatePhaseSum(void);

    void resetAtomMove(int atm);

};

