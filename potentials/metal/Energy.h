/**
 * @brief class to pass energy across the interface. The only place this should be changed is in the field file!!!!!
 * 
 */
#pragma once

#include <fstream>
#include <stdio.h>
#include <iomanip>

#include "Constants.h"

class Energy
{
public:

    double
        /** the total energy of the simulation cell */
        totalEnergy = 0.0,

        /** reciprocal space coulomb energy */
        rcpEnergy = 0.0,

        /** the self energy */
        selfEnergy = 0.0,
        
        /** real space coulomb energy */
        realEnergy = 0.0,
        
        /** real space coulomb energy */
        manyEnergy = 0.0,

        /** two body emprical potentail energy */
        vdwEnergy = 0.0,

        /** the external energy */
        extEnergy = 0.0,
        
        /** converts the time step into that used by MD */
        timeStepUnit = EVTODL,
        
        /** Boltzmann constant for this energy unit*/
        boltzmann = 0.00008617333262145;

    std::string
        /** the label for energy to be printed after energy */
        energyUnit = " eV";

    std::string getEnergyUnit(void)
        {return energyUnit;}

    void setEnergyUnit(std::string unit)
        {energyUnit = unit;}
        
    void setTotalEnergy(double eng)
        { totalEnergy = eng;}

    double getTotalEnergy(void)
        {return rcpEnergy + realEnergy + vdwEnergy + manyEnergy + extEnergy;}

    /**
     * @brief Get the Time Step Unit for Molecular Dynamics
     * 
     * @return double 
     */
    double getTimeStepUnit(void)
        {return timeStepUnit;}

    void setTimeStepUnit(double unit)
        {timeStepUnit = unit;}

    /**
     * @brief Retirns the Boltzmann Constant in the appropriate energy units
     * 
     * @return double 
     */
    double getBoltzmannConstant(void)
        {return boltzmann;}

    void zero(void);

    void printEnergy(int box, std::ofstream& outStream);

    Energy() {};
    virtual ~Energy() {};

    /**
    * @brief Copy constructor. 
    * @param src (const Energy&) : is the old/original object/configuration.
    */
    Energy(const Energy& src);

    /**
    * @brief overloaded assignment operator. Performs a deep copy
    * 
    * @param src (const Energy&) : is the types being copied
    * @return returns a copy of the energies
    */
    Energy& operator=(const Energy& src);


    Energy operator-(const Energy& src);

    Energy operator+(const Energy& src) const;

    Energy& operator+=(const Energy& src);

};
