#pragma once


#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <algorithm>


#include "Species.h"
#include "Memory.h"
#include "Constants.h"

#include "Potential.h"

enum potenkeyword
{
    buckingham = 1,
    morse,
    ljones,
    simple,
    stweber,
    sw2,
    wca,
    ew,
    bond,
    bhm,
    sevnsix,
    potends
};

class TwoBody
{
private:

    const int MAXPAR = 6; // the max number of parameters per potential

    int
        /** the number of potential parameter sets */
        numberPotentials = 0,
        
        /** the maximum number of interactions */
        maxInteraction = 0;

    int
        /** local copy of the number of species */
        numSpec;

    

    std::vector<Potential>
        vdwPots;  // stores the parameters

    double
        delMesh,                   // spacin of points on look up table
        cutoff;                    // real space cutoff

    int
        *meshType,          /**< the type of potential used in the calculation */
        *meshTable;                 // links atom type to look up mesh

    double
        *mesh;             /**< the actual parameters used in the calculation */

public:

    void shortRange(int ltypea, int ltypeb, double r, double rsq, double& energy)
    {
        int
            pot,
            potnum;

        double
            temp1,
            eng1,
            eng2;

        energy = 0.0;
        
        if (r < 1.0e-8)
            return;

        
        pot = ltypea * numSpec + ltypeb;

        potnum = meshTable[pot];

        if (potnum < 0)
            return;

        switch (meshType[potnum])
        {

        case buckingham:   // buckingham potential interaction
        {
            potnum = potnum * MAXPAR;
            eng1 = mesh[potnum] * exp(-r / mesh[potnum+1]);
            eng2 = -mesh[potnum+2] / (rsq * rsq * rsq);
            energy = eng1 + eng2;

            break;

        }

        case morse:   // morse potential interaction
        {
            potnum = potnum * MAXPAR;
            eng1 = exp(-mesh[potnum+1] * (r - mesh[potnum+2]));
            eng2 = mesh[potnum+3] / pow(rsq, 6);
            energy = eng2 + mesh[potnum+0] * (1.0 - eng1) * (1.0 - eng1) - mesh[potnum+0];

            break;
        
        }

        case ljones:   // len-jones potential interaction
        {
            potnum = potnum * MAXPAR;
            temp1 = mesh[potnum+1] / r;
            eng1 = 4.0 * mesh[potnum+0] * pow(temp1, 12.0);
            eng2 = -4.0 * mesh[potnum+0] * pow(temp1, 6.0);

            energy = eng1 + eng2;

            break;
        }
        case bhm: // born-huggins-meyer exp - 6 - 8 potential
        {
            potnum = potnum * MAXPAR;
            double r6 = rsq * rsq * rsq;
            double r8 = r6 * rsq;
            energy = mesh[potnum] * exp(mesh[potnum+1] * (mesh[potnum+2] - r)) - mesh[potnum+3] / r6 - mesh[potnum+4] / r8;

            break;
        }

        case sevnsix: //  severn six
        {
            potnum = potnum * MAXPAR;
            temp1 = mesh[potnum+1] / r;
            eng1 = mesh[potnum] * pow(temp1, 7.0);
            eng2 = -2.0 * mesh[potnum] * pow(temp1, 6.0);

            energy = eng1 + eng2;

            break;

        }

        default:

            std::cout << "\n\n*** error in short range energy";
            std::cout << "\n pot " << ltypea << " " << ltypeb << " " << pot << " " << potnum << " " << meshType[pot];
            std::cout.flush();
            exit(EXIT_FAILURE);

        }

    }

    /**
     * @brief calculates the two-body energy and force between two atoms
     * 
     * @param ltypea (int) : The type for atom A
     * @param ltypeb  (int) : The type for atom B
     * @param r (double) : distance
     * @param rsq (double) : distance squared
     * @param energy (double) : energy
     * @param force (double) : force
     */
    void shortrangeone(int ltypea, int ltypeb, double r, double rsq, double& energy, double& force)
    {
        int
            pot,
            potnum;

        double
            temp1,
            eng1,
            eng2,
            f1,
            f2;

        energy = 0;
        force = 0;

        if (r < 1.0e-8)
            return;

        
        pot = ltypea * numSpec + ltypeb;

        potnum = meshTable[pot];

        if (potnum < 0)
            return;

        switch (meshType[potnum])
        {

        case buckingham:   // buckingham potential interaction
        {
            potnum = potnum * MAXPAR;
            eng1 = mesh[potnum] * exp(-r / mesh[potnum+1]);
            f1 = eng1 / (r * mesh[potnum+1]);
            eng2 = -mesh[potnum+2] / (rsq * rsq * rsq);
            f2 = 6 * eng2 / rsq;
            energy = eng1 + eng2;
            force = f1 + f2;

            break;

        }

        case morse:   // morse potential interaction
        {
            potnum = potnum * MAXPAR;
            eng1 = exp(-mesh[potnum+1] * (r - mesh[potnum+2]));
            energy = mesh[potnum+0] * (1.0 - eng1) * (1.0 - eng1) - mesh[potnum+0];
            force = -(2.0 * mesh[potnum+0] * mesh[potnum+1] * (1.0 - eng1) * eng1);
            force = force / r;
            
            break;
        
        }

        case ljones:   // len-jones potential interaction
        {
            potnum = potnum * MAXPAR;
            temp1 = mesh[potnum+1] / r;
            eng1 = 4.0 * mesh[potnum+0] * pow(temp1, 12.0);
            f1 = 12.0 * eng1 / rsq;
            eng2 = -4.0 * mesh[potnum+0] * pow(temp1, 6.0);
            f2 = 6.0 * eng2 / rsq;

            energy = eng1 + eng2;
            force = f1 + f2;
            force = -1.0 * force;

            break;
        }
        case bhm: // born-huggins-meyer exp - 6 - 8 potential
        {
            potnum = potnum * MAXPAR;
            double r6 = rsq * rsq * rsq;
            double r8 = r6 * rsq;
            energy = mesh[potnum] * exp(mesh[potnum+1] * (mesh[potnum+2] - r)) - mesh[potnum+3] / r6 - mesh[potnum+4] / r8;
            force = r * mesh[potnum] * mesh[potnum+1] * exp(mesh[potnum+1] * (mesh[potnum+2] - r)) - 6.0 * mesh[potnum+3] / r6 - 8.0 * mesh[potnum+4] / r8;

            force = force / rsq;

            break;
        }

        case sevnsix: //  severn six
        {
            potnum = potnum * MAXPAR;
            temp1 = mesh[potnum+1] / r;
            eng1 = mesh[potnum] * pow(temp1, 7.0);
            eng2 = -2.0 * mesh[potnum] * pow(temp1, 6.0);
            double r7 = rsq * rsq * rsq * r;
            double r8 = r7 * r;
            double s7 = pow(mesh[potnum+1], 7.0);
            double s6 = pow(mesh[potnum+1], 6.0);
            f1 = -mesh[potnum] * 7.0 * s7 / r8;
            f2 = -mesh[potnum] * 12.0 * s6 / r7;
            force = f1 + f2;
            force = -1.0 * force / rsq;

            energy = eng1 + eng2;

            break;

        }

        default:

            std::cout << "\n\n*** error in short range energy";
            std::cout << "\n pot " << ltypea << " " << ltypeb << " " << pot << " " << potnum << " " << meshType[pot];
            std::cout.flush();
            exit(EXIT_FAILURE);

        }

    }

    /**
     * @brief makes a copy of two body potentials
     * 
     * @param src (const TwoBody&) : the original potential set
     * @return TwoBody& 
     */
    TwoBody& operator=(const TwoBody& src);

    /**
     * @brief reads in potential parameters
     * 
     * @param instream (std::ifstream&) : input stream
     * @param outstream (std::ofstream&) : output stream
     * @param num (int) : the number of potentials
     */
    void loadPotential(std::ifstream& instream, std::ofstream& outstream, int num);

    /**
     * @brief writes out potential parameters
     * 
     * @param outstream (std::ofstream&) : output stream
     */
    void printPotential(std::ofstream& outstream);

    /**
     * @brief converts the potential name into an integer type
     * 
     * @param keyword (std::string) : the potential name
     * @return key (int) : the type of potential 
     */
    int potenKey(std::string keyword);

    /**
     * @brief calculates the energy mesh from std::vector into a lookup table of parameters
     * 
     * @param spec (const Species&) : atom type data
     * @param shortrangecut (double) : short-range cutoff
     * @param outstream (std::ofstream&) : output stream
     */
    void calculateEnergyMesh(const Species& spec, double shortrangecut, std::ofstream& outstream);

    /**
     * @brief converts energy units of the potential
     * 
     * @param unit (double) : the energy unit for the potential parameter
     */
    void convertPotential(double unit);

    TwoBody();
    virtual ~TwoBody();

    

};

