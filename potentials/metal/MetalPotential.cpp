#include "MetalPotential.h"

std::vector<std::string> split(std::string s);

MetalPotential::MetalPotential()
{

}
MetalPotential::~MetalPotential()
{  
    free_dvector(parRep);
    free_dvector(parDens);
    free_dvector(parEmbed);
    free_dvector(distRep);
    free_dvector(distDens);
    free_dvector(distEmbed);
    free_dvector(derivRep);
    free_dvector(derivDens);
    free_dvector(derivEmbed);
}


MetalPotential::MetalPotential(const MetalPotential& src)
{
    potentialType = src.potentialType;
    numComponents = src.numComponents;
    numPtsRep = src.numPtsRep;
    numPtsDens = src.numPtsDens;
    numPtsEmbed = src.numPtsEmbed;

    startRep = src.startRep;
    finishRep = src.finishRep;
    startDens = src.startDens;
    finishDens = src.finishDens;
    startEmbed = src.startEmbed;
    finishEmbed = src.finishEmbed;
    stepEmbed = src.stepEmbed;
    stepDens = src.stepDens;
    stepRep = src.stepRep;
    potCutoff = src.potCutoff;

    param1 = src.param1;                  // short- range interaction parameters for parameterised potentials
    param2 = src.param2; 
    param3 = src.param3; 
    param4 = src.param4; 
    param5 = src.param5; 
    param6 = src.param6; 
    param7 = src.param7;  
    param8 = src.param8;
    param9 = src.param9;

    parEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
    distEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
    derivEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);

    parRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
    distRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
    derivRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);

    parDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
    distDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
    derivDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);

    for (int i = 0; i < numPtsEmbed; i++)
    {
        parEmbed[i] = src.parEmbed[i];
        derivEmbed[i] = src.derivEmbed[i];
        distEmbed[i] = src.distEmbed[i];
    }

    for (int i = 0; i < numPtsDens; i++)
    {
        parDens[i] = src.parDens[i];
        derivDens[i] = src.derivDens[i];
        distDens[i] = src.distDens[i];
    }
    
    for (int i = 0; i < numPtsRep; i++)
    {
        parRep[i] = src.parRep[i];
        distRep[i] = src.distRep[i];
        derivRep[i] = src.derivRep[i];
    }

    atA = src.atA;
    atB = src.atB;
}

MetalPotential& MetalPotential::operator=(const MetalPotential& src)
{
    if (this == &src) //ie self assignment
    {
        return *this;
    }

    potentialType = src.potentialType;
    numComponents = src.numComponents;
    numPtsRep = src.numPtsRep;
    numPtsDens = src.numPtsDens;
    numPtsEmbed = src.numPtsEmbed;

    startRep = src.startRep;
    finishRep = src.finishRep;
    startDens = src.startDens;
    finishDens = src.finishDens;
    startEmbed = src.startEmbed;
    finishEmbed = src.finishEmbed;
    stepEmbed = src.stepEmbed;
    stepDens = src.stepDens;
    stepRep = src.stepRep;
    potCutoff = src.potCutoff;

    param1 = src.param1;                  // short- range interaction parameters for parameterised potentials
    param2 = src.param2; 
    param3 = src.param3; 
    param4 = src.param4; 
    param5 = src.param5; 
    param6 = src.param6; 
    param7 = src.param7;  
    param8 = src.param8;
    param9 = src.param9;

    if (parEmbed == nullptr) 
        parEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
    if (distEmbed == nullptr) 
        distEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
    if (derivEmbed == nullptr) 
        derivEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);

    if (parRep == nullptr) 
        parRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
    if (distRep == nullptr) 
        distRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
    if (derivRep == nullptr) 
        derivRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);

    if (parDens == nullptr) 
        parDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
    if (distDens == nullptr) 
        distDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
    if (derivDens == nullptr) 
        derivDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);

    for (int i = 0; i < numPtsEmbed; i++)
    {
        parEmbed[i] = src.parEmbed[i];
        derivEmbed[i] = src.derivEmbed[i];
        distEmbed[i] = src.distEmbed[i];
    }

    for (int i = 0; i < numPtsDens; i++)
    {
        parDens[i] = src.parDens[i];
        derivDens[i] = src.derivDens[i];
        distDens[i] = src.distDens[i];
    }
    
    for (int i = 0; i < numPtsRep; i++)
    {
        parRep[i] = src.parRep[i];
        distRep[i] = src.distRep[i];
        derivRep[i] = src.derivRep[i];
    }

    atA = src.atA;
    atB = src.atB;
    
    return *this;
}


void MetalPotential::readEAM(std::string infile, std::ofstream& outStream)
{
    std::string
        dummy,
        keyWord,
        line;
    std::vector<std::string> words;
    std::ifstream inStream;

    inStream.open(infile, std::ios::in);
    if (inStream.fail())              // check to see file is there
    {
        outStream << "\n*** could not find EAM file" << std::endl;
        outStream.flush();
        outStream.close();
        exit(EXIT_FAILURE);
    }

    std::getline(inStream, line); // first lines are dummy

    std::getline(inStream, line);
	words = split(line);
    numComponents = std::stoi(words[0]);

    for (int nc = 0; nc < numComponents; nc++)
    {
                
        std::getline(inStream, line);
	    words = split(line);
        keyWord = words[0];
        transform(keyWord.begin(), keyWord.end(), keyWord.begin(), ::tolower);

        if (keyWord == "pair")
        {
            numPtsRep = std::stoi(words[3]);
            startRep = std::stoi(words[4]);
            finishRep = std::stoi(words[5]);
            stepRep = 1 + (finishRep - startRep) / numPtsRep;

            parRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
            distRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
            derivRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);

            readParameters(parRep, numPtsRep, inStream);

            for (int i = 0; i < numPtsRep; i++)
            {
                distRep[i] = startRep + i * stepRep;
            }
            calculateDerivativeMesh(parRep, derivRep, stepRep, numPtsRep);

            potCutoff = finishRep;
        }
        else if (keyWord == "dens")
        {
            numPtsDens = std::stoi(words[2]);
            startDens = std::stoi(words[3]);
            finishDens = std::stoi(words[4]);
            stepDens = (finishDens - startDens) / numPtsDens;

            parDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
            distDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
            derivDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);

            readParameters(parDens, numPtsDens, inStream);

            for (int i = 0; i < numPtsDens; i++)
            {
                distDens[i] = startDens + i * stepDens;
            }
            calculateDerivativeMesh(parDens, derivDens, stepDens, numPtsDens);

            if (potCutoff < finishDens)
                potCutoff = finishDens;
        }
        else if (keyWord == "embed")
        {
            numPtsEmbed = std::stoi(words[2]);
            startEmbed = std::stoi(words[3]);
            finishEmbed = std::stoi(words[4]);
            stepEmbed = (finishEmbed - startEmbed) / numPtsEmbed;

            parEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
            distEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
            derivEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);

            readParameters(parEmbed, numPtsEmbed, inStream);

            for (int i = 0; i < numPtsEmbed; i++)
            {
                distEmbed[i] = startEmbed + i * stepEmbed;
            }
            calculateDerivativeMesh(parEmbed, derivEmbed, stepEmbed, numPtsEmbed);
        }
        else
        {
            outStream << " wrong line in readEAM " << line << std::endl;
            outStream.flush();
            exit(EXIT_FAILURE);
        }
    }
    
}

void MetalPotential::readLEAM(std::string infile, std::ofstream& outStream)
{
    std::string
        dummy,
        keyWord,
        line;
    std::vector<std::string> words;
    std::ifstream inStream;

    inStream.open(infile, std::ios::in);
    if (inStream.fail())              // check to see file is there
    {
        outStream << "\n*** could not find LAMMPS style EAM file " << infile << std::endl;
        outStream.flush();
        outStream.close();
        exit(EXIT_FAILURE);
    }

    std::getline(inStream, line); // first 4 lines are dummy (3 appear to be for comments)
    std::getline(inStream, line); 
    std::getline(inStream, line); 
    std::getline(inStream, line); 

    std::getline(inStream, line);

	words = split(line);

    numPtsEmbed = std::stoi(words[0]);
    stepEmbed = std::stod(words[1]);
    numPtsRep = numPtsDens = std::stoi(words[2]);
    stepRep = stepDens = std::stod(words[3]);
    potCutoff = std::stod(words[4]);
    

    parEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
    distEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);
    derivEmbed = alloc_dvector(numPtsEmbed, " Repulsive for metal potential", 0.0);

    parRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
    distRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);
    derivRep = alloc_dvector(numPtsRep, " Repulsive for metal potential", 0.0);

    parDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
    distDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);
    derivDens = alloc_dvector(numPtsDens, " Repulsive for metal potential", 0.0);

    std::getline(inStream, line); // dummy line

    //embedding or glue
    readParameters(parEmbed, numPtsEmbed, inStream);
    for (int i = 0; i < numPtsEmbed; i++)
    {
        //std::getline(inStream, line);
	    //words = split(line);
        //parEmbed[i] = std::stod(words[0]);
        distEmbed[i] = i * stepEmbed;
    }

    calculateDerivativeMesh(parEmbed, derivEmbed, stepEmbed, numPtsEmbed);
    //outStream << " finished embed " << line << std::endl;

    //density for calculating rho
    readParameters(parDens, numPtsDens, inStream);
    for (int i = 0; i < numPtsDens; i++)
    {
        //std::getline(inStream, line);
	    //words = split(line);
        //parDens[i] = std::stod(words[0]);
        distDens[i] = i * stepDens;
    }
    calculateDerivativeMesh(parDens, derivDens, stepDens, numPtsDens);
    //outStream << " finished dens " << line << std::endl;

    //pair potential
    readParameters(parRep, numPtsRep, inStream);
    for (int i = 0; i < numPtsRep; i++)
    {
        //std::getline(inStream, line);
	    //words = split(line);
        distRep[i] = i * stepRep;
        //parRep[i] = std::stod(words[0]) / distRep[i];
        parRep[i] = parRep[i] / distRep[i];
    }
    calculateDerivativeMesh(parRep, derivRep, stepRep, numPtsRep);
    //outStream << " finished rep " << line << std::endl;
}

void MetalPotential::readParameters(double* par, int numPts, std::ifstream& inStream)
{
    int
        id = 0;

    std::string
        line;

    std::vector<std::string> words;

    while (id < numPts)
    {
        std::getline(inStream, line);
        words = split(line);

        for (unsigned int i = 0; i < words.size(); i++)
        {
            par[id] = std::stod(words[i]);
            id++;
        }

    }

}

void MetalPotential::calculateDerivativeMesh(double* pot, double* deriv, double dr, int npts)
{
    for(int i = 0; i < npts - 1;i++)
        deriv[i] = (pot[i+1] - pot[i]) / dr;
}

