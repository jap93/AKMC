#include "Energy.h"

Energy::Energy(const Energy& src)
{
    totalEnergy = src.totalEnergy;
    rcpEnergy = src.rcpEnergy;
    realEnergy = src.realEnergy;
    vdwEnergy = src.vdwEnergy;
    manyEnergy = src.manyEnergy;
    vdwEnergy = src.vdwEnergy;
    extEnergy = src.extEnergy;
}

Energy& Energy::operator=(const Energy& src)
{
    
    if (this == &src) //ie self assignment
    {
        return *this;
    }

    totalEnergy = src.totalEnergy;
    rcpEnergy = src.rcpEnergy;
    realEnergy = src.realEnergy;
    vdwEnergy = src.vdwEnergy;
    manyEnergy = src.manyEnergy;
    vdwEnergy = src.vdwEnergy;
    extEnergy = src.extEnergy;
    
    return *this;
}

Energy Energy::operator-(const Energy& src)
{
    Energy result;
    result.totalEnergy = (this->totalEnergy - src.totalEnergy);
    result.rcpEnergy = (this->rcpEnergy - src.rcpEnergy);
    result.realEnergy = (this->realEnergy - src.realEnergy);
    result.selfEnergy = (this->selfEnergy - src.selfEnergy);
    result.manyEnergy = (this->manyEnergy - src.manyEnergy);
    result.vdwEnergy = (this->vdwEnergy - src.vdwEnergy);
    result.extEnergy = (this->extEnergy - src.extEnergy);
    
    return result;
}

Energy Energy::operator+ (const Energy& src) const
{
    Energy result;
    result.totalEnergy = (this->totalEnergy + src.totalEnergy);
    result.rcpEnergy = (this->rcpEnergy + src.rcpEnergy);
    result.realEnergy = (this->realEnergy + src.realEnergy);
    result.selfEnergy = (this->selfEnergy + src.selfEnergy);
    result.manyEnergy = (this->manyEnergy + src.manyEnergy);
    result.vdwEnergy = (this->vdwEnergy + src.vdwEnergy);
    result.extEnergy = (this->extEnergy + src.extEnergy);
    
    return result;
}

Energy& Energy::operator+=(const Energy& src)
{

    this->totalEnergy += src.totalEnergy;
    this->manyEnergy += src.manyEnergy;
    this->realEnergy += src.realEnergy;
    this->vdwEnergy += src.vdwEnergy;
    this->selfEnergy += src.selfEnergy;
    this->extEnergy += src.extEnergy;
    
    return *this;
}

void Energy::zero(void)
{
    totalEnergy = 0.0;
    manyEnergy = 0.0;
    rcpEnergy = 0.0;
    selfEnergy = 0.0;
    realEnergy = 0.0;
    vdwEnergy = 0.0;
    extEnergy = 0.0;
}

void Energy::printEnergy(int box, std::ofstream& outStream)
{
    totalEnergy = rcpEnergy + realEnergy + vdwEnergy + manyEnergy + extEnergy;
    outStream << std::dec << std::scientific << std::setprecision(10)  << "\n\n energies of box " << box 
     << std::setw(13) << "\n total (internal)energy       " << std::setw(25) << totalEnergy
     << std::setw(13) << "\n total madelung energy        " << std::setw(25) << rcpEnergy + realEnergy
     << std::setw(13) << "\n recip space energy           " << std::setw(25) << rcpEnergy 
     << std::setw(13) << "\n real space energy            " << std::setw(25) << realEnergy
     << std::setw(13) << "\n non-bonded vdw energy        " << std::setw(25) << vdwEnergy
     << std::setw(13) << "\n many body energy             " << std::setw(25) << manyEnergy
     << std::setw(13) << "\n external potential energy    " << std::setw(25) << extEnergy << std::endl;
}
