#include "Event.h"

Event::Event()
{

    numberOfAtomsSent = 0;          // initalise to zero
    
    bas2PosX = nullptr;
    bas2PosY = nullptr;
    bas2PosZ = nullptr;
    sadPosX = nullptr;
    sadPosY = nullptr;
    sadPosZ = nullptr;
    globNumber = nullptr;
    centre = nullptr;
    actEnergy = 0.0;
    basin1Energy = 0.0;
    basin2Energy = 0.0;
    saddleEnergy = 0.0;
    rate = 0.0;

}
/*
Event::Event(Basis& bas1, Basis& bas2, Basis& sad, double eact, double prefactor, double temperature, double bas1Energy, double bas2Energy, double sadEnergy)
{
    int
        i;

    numberOfAtoms = bas1.numberOfAtoms;
    actEnergy = eact;
    basin1Energy = bas1Energy;
    basin2Energy = bas2Energy;
    saddleEnergy = sadEnergy;

    //calculate the rate constant
    double beta = 1.0 / (BOLTZMANN * temperature);

    rate = prefactor * exp(-eact * beta);

    globNumber = alloc_uivector(numberOfAtoms, "copy constructor: position");
    bas1PosX = alloc_dvector(numberOfAtoms, "copy constructor: position");
    bas1PosY = alloc_dvector(numberOfAtoms, "copy constructor: position");
    bas1PosZ = alloc_dvector(numberOfAtoms, "copy constructor: position");
    bas2PosX = alloc_dvector(numberOfAtoms, "copy constructor: position");
    bas2PosY = alloc_dvector(numberOfAtoms, "copy constructor: position");
    bas2PosZ = alloc_dvector(numberOfAtoms, "copy constructor: position");
    sadPosX = alloc_dvector(numberOfAtoms, "copy constructor: position");
    sadPosY = alloc_dvector(numberOfAtoms, "copy constructor: position");
    sadPosZ = alloc_dvector(numberOfAtoms, "copy constructor: position");

    for (i = 0; i < numberOfAtoms; i++)
    {
        bas1PosX[i] = bas1.posX[i];
        bas1PosY[i] = bas1.posY[i];
        bas1PosZ[i] = bas1.posZ[i];
        bas2PosX[i] = bas2.posX[i];
        bas2PosY[i] = bas2.posY[i];
        bas2PosZ[i] = bas2.posZ[i];
        sadPosX[i] = sad.posX[i];
        sadPosY[i] = sad.posY[i];
        sadPosZ[i] = sad.posZ[i];
    }
        
}
*/
Event::Event(int* veci, float* vecf, int natoms, double prefactor, double temperature)
{
    numberOfAtomsSent = natoms;

    //calculate the rate constant
    double beta = 1.0 / (BOLTZMANN * temperature);

    bas2PosX = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    bas2PosY = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    bas2PosZ = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosX = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosY = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosZ = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    globNumber = alloc_ivector(numberOfAtomsSent, "copy constructor");
    centre = alloc_dvector(3, "copy contructor", 0.0);

    for (int i = 0; i < numberOfAtomsSent; i++)
    {
        bas2PosX[i] = vecf[6*i];
        bas2PosY[i] = vecf[6*i+1];
        bas2PosZ[i] = vecf[6*i+2];
        sadPosX[i] = vecf[6*i+3];
        sadPosY[i] = vecf[6*i+4];
        sadPosZ[i] = vecf[6*i+5];
        globNumber[i] = int(veci[i]);
    }

    int dim = 6 * numberOfAtomsSent;
    basin1Energy = double(vecf[dim]);
    basin2Energy = double(vecf[dim+1]);
    saddleEnergy = double(vecf[dim+2]);
    actEnergy = double(vecf[dim+3]);  
    
    centre[0] = double(vecf[dim+4]);
    centre[1] = double(vecf[dim+5]);
    centre[2] = double(vecf[dim+6]);

    rate = prefactor * exp(-actEnergy * beta);

}
Event::Event(const Event& src)
{
    int
        i;

    numberOfAtomsSent = src.numberOfAtomsSent;
    actEnergy = src.actEnergy;
    basin1Energy = src.basin1Energy;
    basin2Energy = src.basin2Energy;
    saddleEnergy = src.saddleEnergy;
    rate = src.rate;

    bas2PosX = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    bas2PosY = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    bas2PosZ = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosX = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosY = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosZ = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    globNumber = alloc_ivector(numberOfAtomsSent, "copy constructor");

    for (i = 0; i < src.numberOfAtomsSent; i++)
    {
        bas2PosX[i] = src.bas2PosX[i];
        bas2PosY[i] = src.bas2PosY[i];
        bas2PosZ[i] = src.bas2PosZ[i];
        sadPosX[i] = src.sadPosX[i];
        sadPosY[i] = src.sadPosY[i];
        sadPosZ[i] = src.sadPosZ[i];
        globNumber[i] = src.globNumber[i];
    }

    centre = alloc_dvector(3, "copy contructor", 0.0);
    for (i = 0; i < 3; i++)
        centre[i] = src.centre[i];
        
}


Event::~Event()
{
    free_ivector(globNumber);
    free_fvector(bas2PosX);
    free_fvector(bas2PosY);
    free_fvector(bas2PosZ);
    free_fvector(sadPosX);
    free_fvector(sadPosY);
    free_fvector(sadPosZ);
    free_dvector(centre);
    numberOfAtomsSent = 0;
}

Event& Event::operator=(const Event& src)
{
    int
        i;

    if (this == &src) //ie self assignment
    {
        return *this;
    }

    if (globNumber != nullptr) free_ivector(globNumber);
    if (bas2PosX != nullptr) free_fvector(bas2PosX);
    if (bas2PosY != nullptr) free_fvector(bas2PosY);
    if (bas2PosZ != nullptr) free_fvector(bas2PosZ);
    if (sadPosX != nullptr) free_fvector(sadPosX);
    if (sadPosY != nullptr) free_fvector(sadPosY);
    if (sadPosZ != nullptr) free_fvector(sadPosZ);
    if (centre != nullptr) free_dvector(centre);
    numberOfAtomsSent = 0;


    numberOfAtomsSent = src.numberOfAtomsSent;
    actEnergy = src.actEnergy;
    basin1Energy = src.basin1Energy;
    basin2Energy = src.basin2Energy;
    saddleEnergy = src.saddleEnergy;
    rate = src.rate;

    bas2PosX = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    bas2PosY = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    bas2PosZ = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosX = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosY = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    sadPosZ = alloc_fvector(numberOfAtomsSent, "copy constructor: position");
    globNumber = alloc_ivector(numberOfAtomsSent, "copy constructor");
    centre = alloc_dvector(3, "copy contructor", 0.0);
    
    for (i = 0; i < 3; i++)
        centre[i] = src.centre[i];

    for (i = 0; i < src.numberOfAtomsSent; i++)
    {
        bas2PosX[i] = src.bas2PosX[i];
        bas2PosY[i] = src.bas2PosY[i];
        bas2PosZ[i] = src.bas2PosZ[i];
        sadPosX[i] = src.sadPosX[i];
        sadPosY[i] = src.sadPosY[i];
        sadPosZ[i] = src.sadPosZ[i];
        globNumber[i] = src.globNumber[i];
    }
        
    return *this;

}


void Event::setPositions(Basis& bas, std::ofstream& outStream)
{
    int
        i = 0;


    for (int j = 0; j < numberOfAtomsSent; j++)
    {

        i = globNumber[j];

        if (i >= bas.numberOfAtoms)
        {
            outStream << "\n atom " << i << " global number " << i << " is out of bounds ";
            outStream.flush();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }
        bas.posX[i] = double(bas2PosX[j]);
        bas.posY[i] = double(bas2PosY[j]);
        bas.posZ[i] = double(bas2PosZ[j]);
    }
}

void Event::writeCentre(std::ofstream& outStream)
{
    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
              << std::setprecision(7) << std::setw(15) << "centre " << centre[0] << " " << centre[1] << " " << centre[2] << std::endl;
}

void Event::writeSaddle(const Species& spec, int* atmLabel, std::ofstream& outStream)
{
    int
        lb;

    std::string
        name;

    Element
        ele;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
              << std::setprecision(10) << std::setw(15) << "saddle " << saddleEnergy << "  " << numberOfAtomsSent << std::endl;

    for (int i = 0; i < numberOfAtomsSent; i++)
    {            
        lb = atmLabel[i];
        ele = spec.getSpecies(lb);
        name = ele.name;

        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
            << std::setprecision(10)
            << std::setw(5) << name << "    " << globNumber[i] << "\n"
            << std::setw(15) << sadPosX[i]
            << std::setw(15) << sadPosY[i]
            << std::setw(15) << sadPosZ[i] << std::endl;
    }
}



void Event::writeBasin2(const Species& spec, int* atmLabel, std::ofstream& outStream)
{
    int
        lb;

    std::string
        name;

    Element
        ele;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
              << std::setprecision(10) << std::setw(15) << "basin2 " << basin2Energy << "  " << numberOfAtomsSent << std::endl;

    for (int i = 0; i < numberOfAtomsSent; i++)
    {            
        lb = atmLabel[i];
        ele = spec.getSpecies(lb);
        name = ele.name;

        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
            << std::setprecision(10)
            << std::setw(5) << name << "    " << globNumber[i] << "\n"
            << std::setw(15) << bas2PosX[i]
            << std::setw(15) << bas2PosY[i]
            << std::setw(15) << bas2PosZ[i] << std::endl;
    }
}

std::vector<int> Event::getDisplacedAtoms(const Basis& basin, double basinRadius)
{
    double
        rx, ry, rz, rsq,
        xx, yy, zz,
        radius = basinRadius * basinRadius;

    std::vector<int>
        atmList;

    for (int j = 0; j < numberOfAtomsSent; j++)
    {

        int i = globNumber[j];
    
        rx = bas2PosX[j] - basin.posX[i];
        ry = bas2PosY[j] - basin.posY[i];
        rz = bas2PosZ[j] - basin.posZ[i];

        xx = rx * basin.rcpVector[0] + ry * basin.rcpVector[1] + rz * basin.rcpVector[2];
        yy = rx * basin.rcpVector[3] + ry * basin.rcpVector[4] + rz * basin.rcpVector[5];
        zz = rx * basin.rcpVector[6] + ry * basin.rcpVector[7] + rz * basin.rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * basin.latVector[0] + yy * basin.latVector[3] + zz * basin.latVector[6];
        ry = xx * basin.latVector[1] + yy * basin.latVector[4] + zz * basin.latVector[7];
        rz = xx * basin.latVector[2] + yy * basin.latVector[5] + zz * basin.latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        if (rsq > radius)
            atmList.push_back(i);
      
    }

    return atmList;
}

