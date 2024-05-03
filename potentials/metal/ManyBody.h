#pragma once


#include <fstream>
#include <vector>
#include <iomanip>
#include <string>

#include "Constants.h"
#include "Species.h"
#include "MetalPotential.h"
#include "Memory.h"

enum keymetal
{
    sutchen = 1,
    fnsc,
    xfnsc,
    gupt,
    eam,
    meam,
    metalend
};

class ManyBody
{
private:

    int 
        maxParameter = 9; // the max number of energy values per potential unless EAM where this is changed

    int
        numberPotentials = 0;

    int
        maxInteraction = 0,
        numSpec = 0,
        numPot = 0;

    std::vector<MetalPotential>
        metPots;  // stores the parameters

    int
        maxRepParameter = 0,
        maxDensParameter = 0,
        maxEmbedParameter = 0;
    int
        *meshType = nullptr,
        *meshTable = nullptr;                 // links atom type to look up mesh

    double
        cutoff = 0.0,
        *mesh = nullptr;
        /**distRepEng = nullptr,
        *distDensEng = nullptr,
        *distEmbedEng = nullptr,
        *meshRepEng = nullptr,
        *meshDensEng = nullptr,
        *meshEmbedEng = nullptr,
        *meshRepFrc = nullptr,
        *meshDensFrc = nullptr,
        *meshEmbedFrc = nullptr*/            /**< the actual parameters used in the calculation */

    bool
        potSutChen = false,
        potFinSin = false,
        potGupta = false,
        potEAM = false;
       

    int potenKey(std::string keyWord);

public:

    void calculateEnergyMesh(const Species& ele, double shortrangecut, std::ofstream& outStream);

    void calculateEAMEnergyMesh(const Species& ele, double shortrangecut, std::ofstream& outStream);

    double calculateManyBodyDensEnergy(int ltypea, int ltypeb, double r, double rsq);

    double embed(int ltypea, int ltypeb, double rho);

    void calculateManyBodyPairEnergy(int ltypea, int ltypeb, double r, double rsq, double& vPair);

    void calculateManyBodyPairForce(int ltypea, int ltypeb, double r, double rsq, double& vPair, double& forc);

    void calculateManyBodyForce(int ltypea, int ltypeb, double r, double rsq, double rhoi, double rhoj, double& fDens);

    
    /********************************************************************
    returns pair potential energy and force
    *******************************************************************/
double eamMeshEnergy(double* energyMesh, double* distMesh, double r, double dr, int index);

double eamPairForce(double r, double dr, int potnum, int index)
{
//force_pair: this function returns the pair potential part of the force
    double forc = 0.0;
    //double forc = -(meshRepFrc[potnum + index] + (meshRepFrc[potnum + index + 1] - meshRepFrc[potnum + index])*(r - distRepEng[potnum + index]) / dr);
    return forc;
}

    
    void checkMesh(void);

    ManyBody& operator=(const ManyBody& src);

    ManyBody(const ManyBody& src);

    ManyBody();
    virtual ~ManyBody();

    void loadPotential(int, std::ifstream&, std::ofstream&);          // dimension and load potential parameters

    void calculateDerivativeMesh(double* pot, double* deriv, int atomType);

    void printPotential(std::ofstream&); 	                    // print potentials  

    void convertPotential(void);

};

