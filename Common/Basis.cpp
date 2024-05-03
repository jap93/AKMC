#include "Basis.h"
// basis.cpp: implementation of the basis class.
//
//////////////////////////////////////////////////////////////////////


double get_randonumber(void);
std::vector<std::string> split(std::string s);

Basis::Basis()
{

    numberOfAtoms = 0;          // initalise to zero
    

    atmLabel = nullptr;
    globalNo = nullptr;
    posX = nullptr;
    posY = nullptr;
    posZ = nullptr;
    mass = nullptr;
    charge = nullptr;

    // dimension vectors
    latVector = nullptr;
    rcpVector = nullptr;

    stress = nullptr;
    volume = 0.0;

    numFrozenAtoms = 0;
    frozen = nullptr;

    firstWrite = true;
}


Basis::Basis(const Basis& obas)
{
  
    volume = obas.volume;

    numberOfAtoms = obas.numberOfAtoms;
    numFrozenAtoms = obas.numFrozenAtoms;
    title = obas.title;

    firstWrite = obas.firstWrite;

    atmLabel = alloc_ivector(numberOfAtoms, "copy constructor: atlabel", -1);
    globalNo = alloc_ivector(numberOfAtoms, "copy constructor: atlabel", 0);
    posX = alloc_dvector(numberOfAtoms, "copy constructor: position");
    posY = alloc_dvector(numberOfAtoms, "copy constructor: position");
    posZ = alloc_dvector(numberOfAtoms, "copy constructor: position");
    mass = alloc_dvector(numberOfAtoms, "copy constructor: mass");
    charge = alloc_dvector(numberOfAtoms, "copy constructor: charge");
    frozen = alloc_ivector(numberOfAtoms, "copy constructor frozen");

    latVector = alloc_dvector(9, "lattice vectors");
    rcpVector = alloc_dvector(9, "lattice vectors");

    stress = alloc_dvector(9, "lattice vectors");

    for (int i = 0; i < obas.numberOfAtoms; i++)
    {

        atmLabel[i] = obas.atmLabel[i];
        globalNo[i] = obas.globalNo[i];
        posX[i] = obas.posX[i];
        posY[i] = obas.posY[i];
        posZ[i] = obas.posZ[i];
        mass[i] = obas.mass[i];
        charge[i] = obas.charge[i];
        frozen[i] = obas.frozen[i];
    }

    for (int i = 0; i < 9; i++)
    {

        latVector[i] = obas.latVector[i];
        rcpVector[i] = obas.rcpVector[i];

    }

 
    
}


Basis::~Basis()
{

    free_ivector(atmLabel);
    free_ivector(globalNo);
    free_dvector(posX);
    free_dvector(posY);
    free_dvector(posZ);
    free_dvector(mass);
    free_dvector(charge);
    free_ivector(frozen);

    //deallocate the memory used for lattice vectors
    free_dvector(latVector);
    free_dvector(rcpVector);
    free_dvector(stress);

}

void Basis::allocateArrays(int nAtoms)
{
    numberOfAtoms = nAtoms;
    numFrozenAtoms = 0;
    //dimension atom arrays
    atmLabel = alloc_ivector(numberOfAtoms, "atlabel", -1);
    posX = alloc_dvector(numberOfAtoms, "position x");
    posY = alloc_dvector(numberOfAtoms, "position y");
    posZ = alloc_dvector(numberOfAtoms, "position z");
    mass = alloc_dvector(numberOfAtoms, "mass");
    charge = alloc_dvector(numberOfAtoms, "charge");
    frozen = alloc_ivector(numberOfAtoms, "frozen");
    latVector = alloc_dvector(9, "lattice vectors");
    rcpVector = alloc_dvector(9, "lattice vectors");
    globalNo = alloc_ivector(numberOfAtoms, "copy constructor: atlabel", 0);
}

Basis& Basis::operator=(const Basis& src)
{

    if (this == &src) //ie self assignment
    {
        return *this;
    }

    
    if (atmLabel != nullptr) free_ivector(atmLabel);
    if (globalNo != nullptr) free_ivector(globalNo);
    if (posX != nullptr) free_dvector(posX);
    if (posY != nullptr) free_dvector(posY);
    if (posZ != nullptr) free_dvector(posZ);
    if (mass != nullptr) free_dvector(mass);
    if (charge != nullptr) free_dvector(charge);
    if (frozen != nullptr) free_ivector(frozen);

    //deallocate the memory used for lattice vectors
    if (latVector != nullptr) free_dvector(latVector);
    if (rcpVector != nullptr) free_dvector(rcpVector);
    if (stress != nullptr) free_dvector(stress); 

    numberOfAtoms = src.numberOfAtoms;
    numFrozenAtoms = src.numFrozenAtoms;
    title = src.title;

    firstWrite = src.firstWrite;

    atmLabel = alloc_ivector(numberOfAtoms, "copy constructor: atlabel", -1);
    globalNo = alloc_ivector(numberOfAtoms, "copy constructor: atlabel", 0);
    posX = alloc_dvector(numberOfAtoms, "copy constructor: position");
    posY = alloc_dvector(numberOfAtoms, "copy constructor: position");
    posZ = alloc_dvector(numberOfAtoms, "copy constructor: position");
    mass = alloc_dvector(numberOfAtoms, "copy constructor: mass");
    charge = alloc_dvector(numberOfAtoms, "copy constructor: charge");
    frozen = alloc_ivector(numberOfAtoms, "copy constructor frozen");

    latVector = alloc_dvector(9, "lattice vectors");
    rcpVector = alloc_dvector(9, "lattice vectors");

    //copy all elements over
    volume = src.volume;

    for (int i = 0; i < src.numberOfAtoms; i++)
    {

        atmLabel[i] = src.atmLabel[i];
        globalNo[i] = src.globalNo[i];
        posX[i] = src.posX[i];
        posY[i] = src.posY[i];
        posZ[i] = src.posZ[i];
        mass[i] = src.mass[i];
        charge[i] = src.charge[i];
        frozen[i] = src.frozen[i];
    }

    for (int i = 0; i < 9; i++)
    {

        latVector[i] = src.latVector[i];
        rcpVector[i] = src.rcpVector[i];

    }


    return *this;

}

void Basis::reportBasisDifference(const MPIComms& mpi, Basis& basin2, const Species& spec, double detectBasinRadius, std::ofstream& outStream)
{
    int
        j;

    double
        rx, ry, rz, rsq,
        xx, yy, zz,
        radius = detectBasinRadius * detectBasinRadius;
    
    std::string
        symbol,
        atype;

    Element
        ele;

    for (int i = 0; i < numberOfAtoms; i++)
    {
        j = globalNo[i];

        rx = posX[i] - basin2.posX[j];
        ry = posY[i] - basin2.posY[j];
        rz = posZ[i] - basin2.posZ[j];

        xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
        yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
        zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        if (rsq > radius)
        {

            int lb = atmLabel[i];
            ele = spec.getSpecies(lb);
            symbol = ele.name;
            outStream << "\n " << symbol << " index " << i  << " global no " << j << std::endl;
            outStream << std::scientific << std::setprecision(5) << "  " << rx << "  " << ry << "  " << rz << "  " << sqrt(rsq) << std::endl;
            
        }

    }

}
std::vector<int> Basis::getTargetAtoms(const MPIComms& mpi, Basis& basin2, double detectBasinRadius)
{
    double
        rx, ry, rz, rsq,
        xx, yy, zz,
        radius = detectBasinRadius * detectBasinRadius;
    
    std::string
        symbol,
        atype;

    Element
        ele;

    std::vector<int>
        targetAtoms;

    //if (mpi.getRank() != 0) 
    //    return;

    for (int i = 0; i < numberOfAtoms; i++)
    {

        rx = posX[i] - basin2.posX[i];
        ry = posY[i] - basin2.posY[i];
        rz = posZ[i] - basin2.posZ[i];

        xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
        yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
        zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        if (rsq > radius)
        {

            targetAtoms.push_back(i);
            
        }

    }

    return targetAtoms;
}
std::vector<int> Basis::calculateDefects(const MPIComms& mpi, int* label, double* gridX, double* gridY, double* gridZ, 
                                         double detectBasinRadius, std::ofstream& outStream)
{
    double
        rx, ry, rz, rsq,
        xx, yy, zz,
        radius = detectBasinRadius * detectBasinRadius;

    std::vector<int>
        targetAtoms;

    bool
        found = false;

    for (int i = 0; i < numberOfAtoms; i++)
    {

        found = false;

        for (int j = 0; j < numberOfAtoms; j++)
        {

            rx = posX[j] - gridX[i];
            ry = posY[j] - gridY[i];
            rz = posZ[j] - gridZ[i];

            xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
            yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
            zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            if (rsq < radius && atmLabel[i] == label[j])
            {
                found = true;
            }

        }

        if (!found) //interstitials
        {
            targetAtoms.push_back(i);   
        }

    }

    return targetAtoms;
}
void Basis::inputBasis(int& restartIteration, double& restartTime, const Species& spec, const JobControl& job, std::ifstream& inStream, std::ofstream& outStream)
{
    int
        levCfg,
        num = spec.getNumSpecies();

    std::string
        symbol,
        dummy,
        line;

    bool
        founda;

    Element
        ele;         // temp storage for each element

    restartIteration = 0;
    restartTime = 0.0;

    std::vector<std::string> words;

    std::getline(inStream, title);

    std::getline(inStream, line);
	words = split(line);

		//keyWord = words[0];
		//transform(keyWord.begin(), keyWord.end(), keyWord.begin(), ::tolower);

    dummy = words[0];
    levCfg = std::stoi(dummy);

    dummy = words[2];
    numberOfAtoms = std::stoi(dummy);

    //dimension atom arrays
    atmLabel = alloc_ivector(numberOfAtoms, "atlabel", -1);
    posX = alloc_dvector(numberOfAtoms, "position x");
    posY = alloc_dvector(numberOfAtoms, "position y");
    posZ = alloc_dvector(numberOfAtoms, "position z");
    mass = alloc_dvector(numberOfAtoms, "mass");
    charge = alloc_dvector(numberOfAtoms, "charge");
    frozen = alloc_ivector(numberOfAtoms, "frozen");
    latVector = alloc_dvector(9, "lattice vectors");
    rcpVector = alloc_dvector(9, "lattice vectors");
    globalNo = alloc_ivector(numberOfAtoms, "copy constructor: atlabel", 0);

    std::getline(inStream, line);
	words = split(line);

    dummy = words[0];
    if (dummy == "timeStep") // the rest of the kine is restart information
    {
        dummy = words[1];
        restartIteration = std::stoi(dummy);
        dummy = words[5];
        restartTime = std::stod(dummy);

        //read in the next line ready for the lattice vectors
        std::getline(inStream, line);
	    words = split(line);
        dummy = words[0];
    }
    
    latVector[0] = std::stod(dummy);
    dummy = words[1];
    latVector[1] = std::stod(dummy);
    dummy = words[2];
    latVector[2] = std::stod(dummy);

    std::getline(inStream, line);
	words = split(line);

    dummy = words[0];
    latVector[3] = std::stod(dummy);
    dummy = words[1];
    latVector[4] = std::stod(dummy);
    dummy = words[2];
    latVector[5] = std::stod(dummy);

    std::getline(inStream, line);
	words = split(line);

    dummy = words[0];
    latVector[6] = std::stod(dummy);
    dummy = words[1];
    latVector[7] = std::stod(dummy);
    dummy = words[2];
    latVector[8] = std::stod(dummy);

    setRcpVec();
    
    for (int j = 0; j < numberOfAtoms; j++)
    {

        std::getline(inStream, line);
        
	    words = split(line);

        symbol = words[0];
        //atype = words[1];

        founda = false;
        for (int k = 0; k < num; k++)
        {

            ele = spec.getSpecies(k);

            if (symbol == ele.name)
            {
                charge[j] = ele.charge;
                mass[j] = ele.mass;
                atmLabel[j] = k;
                founda = true;
            }

        }

        if (founda == false)
        {
            outStream << "\n\n*** no atom of type " << symbol << " exists";
            outStream.flush();
            exit(EXIT_FAILURE);
        }

        std::getline(inStream, line);
	    words = split(line);
        dummy = words[0];
        posX[j] = std::stod(dummy);
        dummy = words[1];
        posY[j] = std::stod(dummy);
        dummy = words[2];
        posZ[j] = std::stod(dummy);

        //set the global number
        globalNo[j] = j;

        //dummy = line.ExtractWord(four);
        //frozen[j] = std::stoi(dummy);

        if (levCfg >= 1)
            std::getline(inStream, line);
	
        if (levCfg >= 2)
            std::getline(inStream, line);
	
    } 
/*
    for (int j = 0; j < numberOfAtoms; j++)
    {

        try
        {
            std::getline(inStream, line);                 
        }
        catch(const std::exception& e)
        {
            outStream << "\n offending line " << line << std::endl;
            outStream << e.what() << std::endl;
            outStream.flush();
        }
	    words = split(line);

        symbol = words[0];

        founda = false;
        for (int k = 0; k < num; k++)
        {

            ele = spec.getSpecies(k);

            if (symbol == ele.name)
            {
                charge[j] = ele.charge;
                mass[j] = ele.mass;
                atmLabel[j] = k;
                founda = true;
            }

        }

        if (founda == false)
        {
            outStream << "\n\n*** no atom of type " << symbol << " exists";
            outStream.flush();
            exit(EXIT_FAILURE);
        }

        dummy = words[1];
        posX[j] = std::stod(dummy);
        dummy = words[2];
        posY[j] = std::stod(dummy);
        dummy = words[3];
        posZ[j] = std::stod(dummy);

        //set the global number
        globalNo[j] = j;

        if (levCfg >= 1)
            std::getline(inStream, line);
	
        if (levCfg >= 2)
            std::getline(inStream, line);

    }
*/

    //for (j = 0; j < numberOfAtoms; j++)
    //{
    //    if (frozen[j] == 1)
    //        numFrozenAtoms++;
    //}
    
}

void Basis::deleteBasisMemory(void)
{
    if (numberOfAtoms != 0) numberOfAtoms = 0;
    //delete the memory ready for a new cluster
    free_ivector(atmLabel);
    free_dvector(posX);
    free_dvector(posY);
    free_dvector(posZ);
    free_dvector(mass);
    free_dvector(charge);
    free_ivector(frozen);
    free_dvector(latVector);
    free_dvector(rcpVector);
    free_ivector(globalNo);

}


/** writes out atoms */

void Basis::dumpBasis(const MPIComms& mpi, const Species& spec, double simulationTime, double energy, int iteration, std::ofstream& outStream)
{
    int
        lb;

    std::string
        symbol;

    Element
        ele;

    if (mpi.getRank() != 0) 
        return;    

    if (firstWrite)  // write the header
    {
        outStream << title << std::endl;
        outStream << "0    3   " << numberOfAtoms << std::endl;
    }

    firstWrite = false;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
              << std::setprecision(10) 
              << std::setw(15) << "timeStep " << iteration << "     " << numberOfAtoms << "    0      3   " << simulationTime << "     " << energy << std::endl;
              
    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(10) 
        << std::setw(15) << latVector[0]
        << std::setw(15) << latVector[1]
        << std::setw(15) << latVector[2] << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(10)  
        << std::setw(15) << latVector[3]
        << std::setw(15) << latVector[4]
        << std::setw(15) << latVector[5] << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(10) 
        << std::setw(15) << latVector[6]
        << std::setw(15) << latVector[7]
        << std::setw(15) << latVector[8] << std::endl;
    outStream.flush();

    for (int i = 0; i < numberOfAtoms; i++)
    {

        lb = atmLabel[i];
        ele = spec.getSpecies(lb);
        symbol = ele.name;
        outStream << symbol << "      " << i + 1 << std::endl;
        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                  << std::setprecision(10) << "  " << posX[i] << "  " << posY[i] << "  " << posZ[i] << std::endl;

    }

}

void Basis::writeBasisXYZ(const Species& spec, std::ofstream& outStream)
{
    int
        lb;

    std::string
        symbol;

    Element
        ele;

    outStream << numberOfAtoms << std::endl;
    outStream << " " << std::endl;

    for (int i = 0; i < numberOfAtoms; i++)
    {

        lb = atmLabel[i];
        ele = spec.getSpecies(lb);
        symbol = ele.name;
        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                  << std::setprecision(10) << symbol << "  " << posX[i] << "  " << posY[i] << "  " << posZ[i] << std::endl;

    }

}


void Basis::printBasis(const Species& spec, std::ofstream& outStream)
{
    int
        lb;

    std::string
        symbol;

    Element
        ele;


    outStream << "\n" << "\n number of atoms = " << numberOfAtoms;

    outStream << "\n";

    for (int i = 0; i < numberOfAtoms; i++)
    {

        lb = atmLabel[i];
        ele = spec.getSpecies(lb);
        symbol = ele.name;

        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
            << std::setprecision(10) << "\n" << i + 1
            << std::setw(5) << symbol
            << std::setw(15) << posX[i]
            << std::setw(15) << posY[i]
            << std::setw(15) << posZ[i];

    }

    outStream.flush();
}


void Basis::calculateCharge(std::ofstream& outStream)
{
    double chargeTotal = 0.0;

    for (int i = 0; i < numberOfAtoms; i++)
        chargeTotal += charge[i];

    if (abs(chargeTotal) > 1.0e-8)
    {
        outStream << "*** the system charge = " << charge << std::endl;
        outStream.flush();
        exit(EXIT_FAILURE);
    }

}

bool Basis::checkBasin(Basis& basin2, double detectBasinRadius, double maxBasinRadius)
{
    int
        j;
        
    bool check = false;

    double
        rx, ry, rz, rsq,
        xx, yy, zz,
        radius = detectBasinRadius * detectBasinRadius,
        maxRadius = maxBasinRadius * maxBasinRadius;

    for (int i = 0; i < numberOfAtoms; i++)
    {

        j = globalNo[i];

        rx = posX[i] - basin2.posX[j];
        ry = posY[i] - basin2.posY[j];
        rz = posZ[i] - basin2.posZ[j];

        xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
        yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
        zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;
        if (rsq >= radius)
        {
            check = true;
        }
        if (rsq > maxRadius)  //check for something horrible
        {
            check = false;
            break;   // immediately abort
        }

            

    }

    return check;
}

double Basis::getDistance(double* centre, int i)
{
    double
        rx, ry, rz,
        xx, yy, zz;

    rx = posX[i] - centre[0];
    ry = posY[i] - centre[1];
    rz = posZ[i] - centre[2];

    xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
    yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
    zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

    xx -= rint(xx);
    yy -= rint(yy);
    zz -= rint(zz);

    rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
    ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
    rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

    return sqrt(rx * rx + ry * ry + rz * rz);

}
void Basis::resetSimulationCell()
{
    double
        rx, ry, rz,
        xx, yy, zz;

    for (int i = 0; i < numberOfAtoms; i++)
    {

        rx = posX[i];
        ry = posY[i];
        rz = posZ[i];

        xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
        yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
        zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        posX[i] = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        posY[i] = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        posZ[i] = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];
       
    }

}

void Basis::fracToCart()
{
    double
        rx,
        ry,
        rz;

    for (int i = 0; i < numberOfAtoms; i++)
    {
        rx = posX[i];
        ry = posY[i];
        rz = posZ[i];

        posX[i] = latVector[0] * rx + latVector[3] * ry + latVector[6] * rz;
        posY[i] = latVector[1] * rx + latVector[4] * ry + latVector[7] * rz;
        posZ[i] = latVector[2] * rx + latVector[5] * ry + latVector[8] * rz;

    }

}


void Basis::cartToFrac()
{
    double
        rx,
        ry,
        rz;

    for (int i = 0; i < numberOfAtoms; i++)
    {
        rx = posX[i];
        ry = posY[i];
        rz = posZ[i];

        posX[i] = rcpVector[0] * rx + rcpVector[1] * ry + rcpVector[2] * rz;
        posY[i] = rcpVector[3] * rx + rcpVector[4] * ry + rcpVector[5] * rz;
        posZ[i] = rcpVector[6] * rx + rcpVector[7] * ry + rcpVector[8] * rz;

    }

}

void Basis::displaceAtoms(double* delta)
{
    int
        i,
        j = 0;

    for (i = 0; i < numberOfAtoms; i++)
    {
        if (frozen[i] == 0)  //i.e. it is free to move
        {
            posX[i] += delta[j];
            posY[i] += delta[j+1];
            posZ[i] += delta[j+2];
            j += 3;
            //std::cout << "\ndisplace i " << delta[j] << " " << delta[j+1] << " " << delta[j+2];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////
// functions that freeze atoms within the simulation cell 
////////////////////////////////////////////////////////////////////////////////////////////////
void Basis::unFreezeAtoms()
{
    numFrozenAtoms = 0;

    for (int i = 0; i < numberOfAtoms; i++)
    {
        frozen[i] = 0;   
    }
}

void Basis::freezeAtomTypes(const Species& spec, std::vector<std::string> types)
{
    int
        i, j, k,
        numFroz = types.size(),
        numSpec = spec.getNumSpecies();

    std::string
        typ;

    Element
        ele;

    numFrozenAtoms = 0;

    for (j = 0; j < numFroz; j++)
    {
        typ = types[j];

        for (k = 0; k < numSpec; k++)
        {
            ele = spec.getSpecies(k);

            if (typ == ele.name)
            {

                for (i = 0; i < numberOfAtoms; i++)
                {
                    if (atmLabel[i] == int(k))
                    {
                        frozen[i] = 1;   // freeze the atom as it is of the correct type
                        numFrozenAtoms++;
                    }
                }

            }

        }

    }
}

void Basis::freezeRegion(double maskLow, double maskHigh)
{
    for (int i = 0; i < numberOfAtoms; i++)
    {
        if (posX[i] >= maskLow and posX[i] <= maskHigh)
        {
            numFrozenAtoms++;
            frozen[i] = 1;
        }
    }
}

void Basis::freezeForces(double* __restrict__ forceX, double* __restrict__ forceY, double* __restrict__ forceZ)
{

    if (numFrozenAtoms == 0)
        return;

    for (int i = 0; i < numberOfAtoms; i++)
    {
        if (frozen[i] == 1)
        {
            forceX[i] = 0.0;
            forceY[i] = 0.0;
            forceZ[i] = 0.0;
        }
    } 
}

////////////////////////////////////////////////////////////////////////////////////////////////
// functions acting on the simulation cell 
////////////////////////////////////////////////////////////////////////////////////////////////

void Basis::setRcpVec(void)
{ 
    rcpVector[0] = (latVector[4] * latVector[8] - latVector[7] * latVector[5]);
    rcpVector[1] = (latVector[6] * latVector[5] - latVector[3] * latVector[8]);
    rcpVector[2] = (latVector[3] * latVector[7] - latVector[6] * latVector[4]);

    rcpVector[3] = (latVector[7] * latVector[2] - latVector[1] * latVector[8]);
    rcpVector[4] = (latVector[0] * latVector[8] - latVector[6] * latVector[2]);
    rcpVector[7] = (latVector[6] * latVector[1] - latVector[0] * latVector[7]);

    rcpVector[6] = (latVector[1] * latVector[5] - latVector[4] * latVector[2]);
    rcpVector[5] = (latVector[2] * latVector[3] - latVector[5] * latVector[0]);
    rcpVector[8] = (latVector[0] * latVector[4] - latVector[3] * latVector[1]);
    
    double determinant = latVector[0] * rcpVector[0] + latVector[1] * rcpVector[3] + latVector[2] * rcpVector[6];
    double fact = 1.0;

    if (abs(determinant) > 0.0)
        fact = 1.0 / determinant;

    for (int i = 0; i < 9; i++)
        rcpVector[i] *= fact;

}

void Basis::printRcpVec(std::ofstream& outStream)
{

    outStream << "\n\n reciprocal lattice vectors";

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(5) << "\n"
        << std::setw(15) << rcpVector[0]
        << std::setw(15) << rcpVector[1]
        << std::setw(15) << rcpVector[2] << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(5) 
        << std::setw(15) << rcpVector[3]
        << std::setw(15) << rcpVector[4]
        << std::setw(15) << rcpVector[5] << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(5) 
        << std::setw(15) << rcpVector[6]
        << std::setw(15) << rcpVector[7]
        << std::setw(15) << rcpVector[8] << std::endl;

}

double Basis::cellSize()
{
    double
        bc1,
        bc2,
        bc3,
        vol;

    bc1 = latVector[4] * latVector[8] - latVector[5] * latVector[7];
    bc2 = latVector[5] * latVector[6] - latVector[3] * latVector[8];
    bc3 = latVector[3] * latVector[7] - latVector[4] * latVector[6];

    vol = latVector[0] * bc1 + latVector[1] * bc2 + latVector[2] * bc3;

    if (vol < 0.0)
        vol = -vol;

    volume = vol;

    return vol;
}





/* returns the minimum dimensions of lattice vectors */
double Basis::minDimension()
{
    double
        cellx,                        // cell dimensions in x y + z
        celly,
        cellz,
        cellmin;

    cellx = pow(latVector[0], 2.0) + pow(latVector[3], 2.0) + pow(latVector[6], 2.0);

    celly = pow(latVector[1], 2.0) + pow(latVector[4], 2.0) + pow(latVector[7], 2.0);

    cellz = pow(latVector[2], 2.0) + pow(latVector[5], 2.0) + pow(latVector[8], 2.0);

    cellmin = cellx;

    if (celly < cellx)
        cellmin = celly;

    if (cellz < celly)
        cellmin = cellz;

    return sqrt(cellmin);

}


void Basis::printVectors(std::ofstream& outStream)
{

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(10) 
        << std::setw(15) << latVector[0]
        << std::setw(15) << latVector[1]
        << std::setw(15) << latVector[2] << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(10)  
        << std::setw(15) << latVector[3]
        << std::setw(15) << latVector[4]
        << std::setw(15) << latVector[5] << std::endl;

    outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << std::setprecision(10) 
        << std::setw(15) << latVector[6]
        << std::setw(15) << latVector[7]
        << std::setw(15) << latVector[8] << std::endl;

}

void Basis::inputLatticeVectors(std::ifstream& inStream)
{
    std::string
        dummy,
        line;

    std::vector<std::string> words;

   /* line.GetLine(inStream, 120, '\n');

    dummy = words[0];
    latVector[0] = std::stod(dummy);
    dummy = words[1];
    latVector[1] = std::stod(dummy);
    dummy = words[2];
    latVector[2] = std::stod(dummy);

    line.GetLine(inStream, 120, '\n');

    dummy = words[0];
    latVector[3] = std::stod(dummy);
    dummy = words[1];
    latVector[4] = std::stod(dummy);
    dummy = words[2];
    latVector[5] = std::stod(dummy);

    line.GetLine(inStream, 120, '\n');

    dummy = words[0];
    latVector[6] = std::stod(dummy);
    dummy = words[1];
    latVector[7] = std::stod(dummy);
    dummy = words[2];
    latVector[8] = std::stod(dummy); */

    std::getline(inStream, line);
	words = split(line);

    dummy = words[0];
    latVector[0] = std::stod(dummy);
    dummy = words[1];
    latVector[3] = std::stod(dummy);
    dummy = words[2];
    latVector[6] = std::stod(dummy);

    std::getline(inStream, line);
	words = split(line);

    dummy = words[0];
    latVector[1] = std::stod(dummy);
    dummy = words[1];
    latVector[4] = std::stod(dummy);
    dummy = words[2];
    latVector[7] = std::stod(dummy);

    std::getline(inStream, line);
	words = split(line);

    dummy = words[0];
    latVector[2] = std::stod(dummy);
    dummy = words[1];
    latVector[5] = std::stod(dummy);
    dummy = words[2];
    latVector[8] = std::stod(dummy);

    setRcpVec();

}

/*
void Basis::sendPositions(const MPIComms& mpi, int tag, std::ofstream& outStream)
{
    int
        rank = mpi.getRank(),
        numProcs = mpi.getNumProcs(),
        grpRank = mpi.getWorkGroupRank(),
        grp = mpi.getGroupID(rank),
        iGrp,
        //tag = 632,  //tag is used to indicate where the positions have come from
        ierr;

    double
        *buffer;

    buffer = alloc_dvector(3 * numberOfAtoms, " mpi :: broadcast basin", 0.0);

    //copy the positions to the buffer
    int j = 0;
    for (int i = 0; i < numberOfAtoms; i++)
    {
        buffer[j] = posX[i];
        buffer[j+1] = posY[i];
        buffer[j+2] = posZ[i];
        j += 3;
    }
    if (grpRank == 0)
    {
    for (int i = 0; i < numProcs; i++)
    {
        iGrp = mpi.getGroupID(i);
        if (i == rank || grp == iGrp)
            continue;

        outStream << "\n sending basin from rank, group " << rank << " " << grp << " to " << i << " " << iGrp << std::endl;
        outStream.flush();
        ierr = MPI_Send(buffer, 3 * numberOfAtoms, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }
    }
    
    free_dvector(buffer);
}

void Basis::receivePositions(const MPIComms& mpi, int tag, std::ofstream& outStream)
{
    int
        rank = mpi.getRank(),
        grpRank = mpi.getWorkGroupRank(),
        idGrp = mpi.getGroupID(rank),
        ierr;

    double
        *buffer;

    MPI_Status
        status;

    buffer = alloc_dvector(3 * numberOfAtoms, " mpi :: receive position ", 0.0);

    ierr = MPI_Recv(buffer, 3 * numberOfAtoms, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

    outStream << "\n ready to receive positions (rank, grp, grpRank) " << rank << " " << idGrp << " " << grpRank << std::endl;

    //copy the positions from the buffer
    int j = 0;
    for (int i = 0; i < numberOfAtoms; i++)
    {
        posX[i] = buffer[j];
        posY[i] = buffer[j+1];
        posZ[i] = buffer[j+2];
        j += 3;
    }    

    free_dvector(buffer);
    outStream << "\n received positions (rank, grp, grpRank) " << rank << " " << idGrp << " " << grpRank << std::endl;
} 
*/
std::vector<int> Basis::mapAtomsToVector(double* centre, double radius) const
{
    double
        rx, ry, rz,
        xx, yy, zz,
        rsq,
        radiusSqrd = radius * radius;

    std::vector<int>
        atoms;

    for (int i = 0; i < numberOfAtoms; i++)
    {

        rx = posX[i] - centre[0];
        ry = posY[i] - centre[1];
        rz = posZ[i] - centre[2];

        xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
        yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
        zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;
        if (rsq < radiusSqrd)
        {
            atoms.push_back(i);
        }
    }

    return atoms;
}

int* Basis::calculateCoordNumber(int typ1, int typ2, double coordDistance) const
{
    double
        ax, ay, az,
        bx, by, bz,
        rx, ry, rz,
        xx, yy, zz,
        rsq,
        radius = coordDistance * coordDistance;

    int* coordNumber = alloc_ivector(numberOfAtoms, " coordination number ", 0);

    for (int i = 0; i < numberOfAtoms; i++)
    {

        if (atmLabel[i] != typ1)
            continue;

        ax = posX[i];
        ay = posY[i];
        az = posZ[i];

        for (int j = 0; j < numberOfAtoms; j++)
        {
        
            if (atmLabel[j] != typ2 || i == j)
                continue;

            bx = posX[j];
            by = posY[j];
            bz = posZ[j];

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector[0] + ry * rcpVector[1] + rz * rcpVector[2];
            yy = rx * rcpVector[3] + ry * rcpVector[4] + rz * rcpVector[5];
            zz = rx * rcpVector[6] + ry * rcpVector[7] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            if (rsq <= radius)
                coordNumber[i] += 1;

        }

    }

    return coordNumber;
}
