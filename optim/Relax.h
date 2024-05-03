// minim.h: interface for the relax class
// this controls input of potential parameters
// and the calculation of minim
//
//////////////////////////////////////////////////////////////////////
#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>


#include "JobControl.h"
#include "Memory.h"
#include "Status.h"
#include "MPICommunicator.h"
#include "Basis.h"
#include "RelaxCntrl.h"
#include "Field.h"

class Relax
{

private:

    bool
        euler = true;

    /** 
     * @brief Relax a configuration using the FIRE method.
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param bas (Grid&) : atomic information for the config
     * @param fld (Field*) : the method/style of force calculation
     * @param natoms (int) : the number of atoms
     * @param rCntrl (const RelaxCntrl&) : parameters to control the relaxation
     * @param status (Status&) : validity of the calculation
     * @param eng (Energy*) : pointer to storage of energy
     * @param outStream (ofstream&) : output stream for the workgroup
     */
    void FIRERelax(const MPIComms& mpi, Basis& bas, Field* fld, int natoms, 
                   const RelaxCntrl& rCntrl, Status& status, Energy& eng, std::ofstream& outStream);

       /** 
     * @brief Relax a configuration using the FIRE method.
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param bas (Grid&) : atomic information for the config
     * @param fld (Field*) : the method/style of force calculation
     * @param natoms (int) : the number of atoms
     * @param rCntrl (const RelaxCntrl&) : parameters to control the relaxation
     * @param status (Status&) : validity of the calculation
     * @param eng (Energy*) : pointer to storage of energy
     * @param outStream (ofstream&) : output stream for the workgroup
     */
    void FIRE2Relax(const MPIComms& mpi, Basis& bas, Field* fld, int natoms, 
                   const RelaxCntrl& rCntrl, Status& status, Energy& eng, std::ofstream& outStream);

    /**
     * @brief calculate the vector norm and maximum
     * 
     * @param vecX (double*) : x component of vector
     * @param vecY (double*) : y component of vector
     * @param vecZ (double*) : z component of vector
     * @param vecNorm (double&) : norm of the vector
     * @param vecMax (double&) : maximum of the vector
     * @param n (int) : length of the vectors
     */
    void vectorNormMax(double* vecX, double* vecY, double* vecZ, double& vecNorm, double& vecMax, int n);

    /**
     * @brief dot product of two vectors
     * 
     * @param x1 (double*) : x component of vector 1
     * @param y1 (double*) : y component of vector 1
     * @param z1 (double*) : z component of vector 1
     * @param x2 (double*) : x component of vector 2
     * @param y2 (double*) : y component of vector 2
     * @param z2 (double*) : z component of vector 2
     * @param n (int) : length of the vectors
     * @return dotProduct (double) 
     */
    double dotProduct(double* x1, double* y1, double* z1, double* x2, double* y2, double* z2, int n);    

    /**
     * @brief 
     * 
     * @param velocityX (double*) : x component of vector
     * @param velocityY (double*) : y component of vector
     * @param velocityZ (double*) : z component of vector
     * @param natoms (int) : number of atoms/length of the vectors
     */
    void centreOfMass(double* velocityX, double* velocityY, double* velocityZ, int natoms);

    

public:
	
    // constructor
    Relax();

    // default destructor
    virtual ~Relax();

    /** 
     * @brief driver routine for relaxation/energy minimisation of structure.
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param bas (Grid&) : atomic information for the config
     * @param fld (Field*) : the method/style of force calculation
     * @param natoms (int) : the number of atoms
     * @param rCntrl (const RelaxCntrl&) : parameters to control the relaxation
     * @param eng (Energy*) : pointer to object containing the energy of the system
     * @param status (Status&) : validity of the calculation
     * @param outStream (ofstream&) : output stream for the workgroup
     */
    void RelaxStructure(const MPIComms& mpi, Basis& bas, Field* fld,  
                        int natoms, const RelaxCntrl& rCntrl, Energy& eng, Status& status, std::ofstream& outStream);

};
