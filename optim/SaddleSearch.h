/**********************************************************
* Virtual class for the calculation of SaddlePoints
*/

#pragma once

#include "Memory.h"
#include "Basis.h"
#include "JobControl.h"
#include "Relax.h"
#include "Status.h"

#include "Field.h"
#include "Energy.h"

class SaddleSearch
{

protected:

    int
        /* the dimension of the vector arrays used within ART once the atoms have been mapped to vectors*/
        vectorSize = 0;

    double
        *rVector = nullptr;
        
    /**
     * @brief displaces the atoms by a random amount
     * 
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param bas (Basis&) : atomic information
     * @param initialDisplacement (double) : size of displacement
     * @param outStream (ofstream&) : output stream for the workgroup
     * @return r (double*) : vector containing the random displacements 
     */
    void makeRandomDisplacement(const MPIComms& mpi, Basis& bas, double initialDisplacement, std::ofstream& outStream);

    /**
     * @brief displaces the atoms by a random amount that is weighted by a Gaussian at centre
     * 
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param bas (Basis&) : atomic information
     * @param centre (double*) : epicentre of the displacement
     * @param initialDisplacement (double) : size of displacement
     * @param gaussWidth (double) : the width of the Gaussian
     * @param outStream (ofstream&) : output stream for the workgroup
     * @return r (double*) : vector containing the random displacements 
     */
    void makeGaussDisplacement(const MPIComms& mpi, Basis& bas, double* centre, double initialDisplacement, double gaussWidth, std::ofstream& outStream);

    /**
     * @brief determines which atoms are activated and returns the centre for a cluster (zero if all atoms are activated)
     * 
     * @param mpi (const MPIComms&) : mpi communication controller
     * @param bas (Basis&) : atomic information
     * @param spec (const Species&) : element information
     * @param job (const JobControl&) : run parameters controlling execution
     * @param centre (double*) : the midpoint of the region used to calculate saddle points
     * @param outStream (ofstream&) : output stream for the workgroup
     */
    void constructActiveRegion(const MPIComms& mpi, Basis& bas, const Species& spec, const JobControl& job, double* centre, std::ofstream& outStream);

/**
* @brief calculates dot product of two vectors
* @param v1 (double*)  : the vector 1
* @param v2 (double*)  : the vector 2
* @param n (int) : the number of elements in the vector
* @return d (double) : the vector dot product
*/    
double dotProduct(double* __restrict__ v1, double* __restrict__ v2, int n)
{
    double
       d = 0.0;

    for (int i = 0; i < n; i++)
        d = d + v1[i] * v2[i];

    return d;
}    

/**
* @brief returns the norm of a vector
* @param v (double*)  : the actual vector
* @param n (int) : the number of elements in the vector
* @return norm (double) : the vector norm
*/
double vectorNorm(double* __restrict__ v, int n)
{
    double
        norm = 0.0;

    for (int i = 0; i < n; i++)
        norm += v[i] * v[i];

    return sqrt(norm);
}

/**
* @brief returns the both the norm and the maximum of a vector
* @param v (double*)  : the actual vector
* @param vecNorm (double) : the vector norm
* @param vecMax (double) : the maximum of the vector
* @param n (int) : the number of elements in the vector
* @return norm (double) : the vector norm
**/
void vectorNormMax(double* __restrict__ vec, double& vecNorm, double& vecMax, int n)
{
    int
        i;
   
    double
        v = 0.0,
        vecSq = 0.0;

    vecNorm = 0.0;
    vecMax = 0.0;

    for (i = 0; i < n; i++)
    {
        v = vec[i];
        vecSq = v * v;
   
        if (vecSq > vecMax) vecMax = vecSq;
      
        vecNorm += vecSq;
    }
   
    vecNorm = sqrt(vecNorm);
    vecMax = sqrt(vecMax);
}

/**
 * @brief Calculates the centre of mass of the vector (the vector is updated)
 * 
 * @param vec (double*) : vector that is used to calculate the centre of mass
 * @param n (int) : length of the vector
 */
void centreOfMass(double* vec, int n)
{
    int
        i, j,
        natoms = n / 3;

    double
        sumX = 0.0,
        sumY = 0.0,
        sumZ = 0.0;

    j = 0;  //take account that atom positions are in separate vectors

    //determine the centre of mass
    for (i = 0; i < natoms; i++)
    {
        sumX += vec[j];
        sumY += vec[j+1];
        sumZ += vec[j+2];

        j += 3;
    }
    
    sumX /= natoms;
    sumY /= natoms;
    sumZ /= natoms;

    //subtract the total
    j = 0;
    for (i = 0; i < natoms; i++)
    {
        vec[j] -= sumX;
        vec[j+1] -= sumY;
        vec[j+2] -= sumZ;

        j += 3;
    }
    
}


public:

    /** 
     * @brief findSaddlePoint is the main routine to determine the activation energy and saddle point of a configuration and has
     * access points in the child classes for MD and ART (and other methods). Some child classes also have publicly exposed classes
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
                                 std::ofstream& outStream) = 0;

                          
    SaddleSearch();
    virtual ~SaddleSearch();
};
