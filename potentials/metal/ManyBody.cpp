#include "ManyBody.h"
//////////////////////////////////////////////////////////////////////
// construction/destruction
//////////////////////////////////////////////////////////////////////
std::vector<std::string> split(std::string s);

ManyBody::ManyBody()
{

}

ManyBody::~ManyBody()
{  
    /*free_dvector(distRepEng);
    free_dvector(distDensEng);
    free_dvector(distEmbedEng);
    free_dvector(meshRepEng);
    free_dvector(meshDensEng);
    free_dvector(meshEmbedEng);
    free_dvector(meshRepFrc);
    free_dvector(meshDensFrc);
    free_dvector(meshEmbedFrc); */
    free_ivector(meshTable);  
    free_ivector(meshType);  
    free_dvector(mesh);             
}

ManyBody::ManyBody(const ManyBody& src)
{
    numberPotentials = src.numberPotentials;
    numSpec = src.numSpec;
    maxInteraction = src.maxInteraction;

    if (meshTable == nullptr) meshTable = alloc_ivector(maxInteraction, "two body table", 0);
    if (meshType == nullptr) meshType = alloc_ivector(maxInteraction, "two body table", 0);
    if (mesh == nullptr) mesh = alloc_dvector(maxInteraction * maxParameter, "two body table", 0.0);

    for (unsigned int i = 0; i < maxInteraction; i++)
    {
        meshTable[i] = src.meshTable[i];
        meshType[i] = src.meshType[i];
    }

    for (unsigned int i = 0; i < maxInteraction * maxParameter; i++)
        mesh[i] = src.mesh[i];

    for (unsigned int i = 0; i < metPots.size(); i++)
        metPots.push_back(src.metPots[i]);
    
}
ManyBody& ManyBody::operator=(const ManyBody& src)
{

    if (this == &src) //ie self assignment
    {
        return *this;
    }

    numberPotentials = src.numberPotentials;
    numSpec = src.numSpec;
    maxInteraction = src.maxInteraction;

    if (meshTable == nullptr) 
        meshTable = alloc_ivector(maxInteraction, "two body table", 0);
    if (meshType == nullptr)
        meshType = alloc_ivector(maxInteraction, "two body table", 0);
    if (mesh == nullptr)
        mesh = alloc_dvector(maxInteraction * maxParameter, "two body table", 0.0);

    for (unsigned int i = 0; i < maxInteraction; i++)
    {
        meshTable[i] = src.meshTable[i];
        meshType[i] = src.meshType[i];
    }

    for (unsigned int i = 0; i < maxInteraction * maxParameter; i++)
        mesh[i] = src.mesh[i];

    for (unsigned int i = 0; i < metPots.size(); i++)
        metPots.push_back(src.metPots[i]);
    
    return *this;
}

void ManyBody::loadPotential(int num, std::ifstream& inStream, std::ofstream& outStream)
{
    int
        keyval;

    std::string
        fileName,
        finished = " ",
        line,
        dummy,
        keyWord;            // directive read as temp

    MetalPotential
        tpot;

    

    numberPotentials = num;

    //metPots = std::vector<MetalPotential>(numberPotentials);  // dimension the array
    std::vector<std::string> words;

    for (int pot = 0; pot < numberPotentials; pot++)
    {

        std::getline(inStream, line);

	    words = split(line);

	    keyWord = words[0];
	    transform(keyWord.begin(), keyWord.end(), keyWord.begin(), ::tolower);

        if (keyWord == "close")
        {
            outStream << "\n it looks like the wrong number of potentials were given!";
            outStream.flush();
            exit(EXIT_FAILURE);
        }

        keyval = potenKey(keyWord);                  // change keyWord to int

        if (keyval == 0)
        {                          // check to see if                                          // valid keyWord
            outStream << "\n\n***potential " << keyWord << " not found";
            outStream.flush();
            exit(EXIT_FAILURE);
        }


        switch (keyval)
        {

            case sutchen:
                tpot.potentialType = keyval;

                potSutChen = true;

                if (potFinSin == true || potGupta == true || potEAM == true)
                {
                    outStream << "\n many-body potential types can not be mixed" << std::endl;
                    exit(EXIT_FAILURE);
                }

                std::getline(inStream, line);
		        words = split(line);

                tpot.atA = words[0];
                tpot.atB = words[1];

                dummy = words[2];
                tpot.param1 = std::stod(dummy);
                dummy = words[3];
                tpot.param2 = std::stod(dummy);
                dummy = words[4];
                tpot.param3 = std::stod(dummy);
                dummy = words[5];
                tpot.param4 = std::stod(dummy);
                dummy = words[6];
                tpot.param5 = std::stod(dummy);

                break;

            case fnsc:
                tpot.potentialType = keyval;

                potFinSin = true;

                if (potSutChen == true || potGupta == true || potEAM == true)
                {
                    outStream << "\n many-body potential types can not be mixed" << std::endl;
                    exit(EXIT_FAILURE);
                }

                std::getline(inStream, line);
		        words = split(line);

                tpot.atA = words[0];
                tpot.atB = words[1];

                dummy = words[2];
                tpot.param1 = std::stod(dummy);
                dummy = words[3];
                tpot.param2 = std::stod(dummy);
                dummy = words[4];
                tpot.param3 = std::stod(dummy);
                dummy = words[5];
                tpot.param4 = std::stod(dummy);
                dummy = words[6];
                tpot.param5 = std::stod(dummy);
                dummy = words[7];
                tpot.param6 = std::stod(dummy);
                dummy = words[8];
                tpot.param7 = std::stod(dummy);

                break;

            case xfnsc:
                tpot.potentialType = keyval;

                potFinSin = true;

                if (potSutChen == true || potGupta == true || potEAM == true)
                {
                    outStream << "\n many-body potential types can not be mixed" << std::endl;
                    exit(EXIT_FAILURE);
                }

                std::getline(inStream, line);
		        words = split(line);

                tpot.atA = words[0];
                tpot.atB = words[1];

                dummy = words[2];
                tpot.param1 = std::stod(dummy);
                dummy = words[3];
                tpot.param2 = std::stod(dummy);
                dummy = words[4];
                tpot.param3 = std::stod(dummy);
                dummy = words[5];
                tpot.param4 = std::stod(dummy);
                dummy = words[6];
                tpot.param5 = std::stod(dummy);
                dummy = words[7];
                tpot.param6 = std::stod(dummy);
                dummy = words[8];
                tpot.param7 = std::stod(dummy);
                dummy = words[9];
                tpot.param8 = std::stod(dummy);
                dummy = words[10];
                tpot.param9 = std::stod(dummy);

                break;

            case gupt:
                tpot.potentialType = keyval;

                potGupta = true;

                if (potFinSin == true || potSutChen == true || potEAM == true)
                {
                    outStream << "\n many-body potential types can not be mixed" << std::endl;
                    exit(EXIT_FAILURE);
                }

                std::getline(inStream, line);
		        words = split(line);

                tpot.atA = words[0];
                tpot.atB = words[1];

                dummy = words[2];
                tpot.param1 = std::stod(dummy);
                dummy = words[3];
                tpot.param2 = std::stod(dummy);
                dummy = words[4];
                tpot.param3 = std::stod(dummy);
                dummy = words[5];
                tpot.param4 = std::stod(dummy);
                dummy = words[6];
                tpot.param5 = std::stod(dummy);

                break;

            case eam:
                tpot.potentialType = keyval;

                potEAM = true;

                if (potFinSin == true || potSutChen == true || potGupta == true)
                {
                    outStream << "\n many-body potential types can not be mixed" << std::endl;
                    exit(EXIT_FAILURE);
                }

                std::getline(inStream, line);
		        words = split(line);
                outStream << "\n eam line: " << line << std::endl;
                tpot.atA = words[0];
                tpot.atB = words[1];
                fileName = words[2];

                dummy = "lammps";

                if (words.size() > 3)
                    dummy = words[3];

                if (dummy == "lammps")
                {
                    outStream << "lammps " << fileName << std::endl;
                    tpot.readLEAM(fileName, outStream);
                }
                else if (dummy == "dlpoly")
                {
                    tpot.readEAM(fileName, outStream);
                }
                else
                {
                    outStream << "\n unrecognised EAM format" << std::endl;
                    exit(EXIT_FAILURE);
                }

                break;

            default:
                outStream << "\n unrecognised many-body potential type" << std::endl;
                exit(EXIT_FAILURE);
            
        }

        metPots.push_back(tpot);
        
    }
	
}
/********************************************************
 convert the potentials to internal energy units
 ********************************************************/
void ManyBody::convertPotential()
{
    MetalPotential
        tpot;        // temp storage for a potential

    double
        unit = 1.0;
/*
    for (int pot = 0; pot < numberPotentials; pot++)
    {
        tpot = metPots[pot];

        

        metPots[pot] = tpot;  // change stored value

    }
*/
}

/********************************************************
 print out all the potentials to a file
 ********************************************************/
void ManyBody::printPotential(std::ofstream& outStream)
{
    MetalPotential
        tpot;        // temp storage for a potential


    outStream << "\n\n******************************************"
              << "\n                many-body potential data"
              << "\n\n******************************************" << std::endl;
/*
    for (int pot = 0; pot < numberPotentials; pot++)
    {
        tpot = metPots[pot];

        outStream << "\n\n      potential          atom a         atom b" << std::endl;


        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
            << "\n"
            << std::setw(13) << tpot.atA
            << std::setw(13) << tpot.atB << std::endl;


    }
*/

    outStream << "\n leaving print potentials" << std::endl;

}

void ManyBody::calculateManyBodyPairEnergy(int ltypea, int ltypeb, double r, double rsq, double& vPair)
{
    int
        pot,
        idx,
        potnum = -1;

    double
        dr,
        arij;

    vPair = 0.0;
        
    if (r < 1.0e-8)
        return;

        
    pot = ltypea * numSpec + ltypeb;

    potnum = meshTable[pot];

    if (potnum < 0)
        return;

    //screeening function
    
    //both potentials are the same here
    switch (meshType[potnum])
    {
        

        case sutchen:
            potnum = potnum * maxParameter;
            arij = mesh[potnum+1] / r;
            //scrn = pow(mesh[potnum]*mesh[potnum+4],2.0);
            vPair = mesh[potnum] * pow(arij, mesh[potnum+2]);
            
            break;

        case fnsc:
            potnum = potnum * maxParameter;

            if (r <= mesh[potnum+3])
                vPair = (mesh[potnum] + mesh[potnum+1] * r + mesh[potnum+2] * rsq) * pow(r - mesh[potnum+3], 2.0);
            break;

        case xfnsc:
            potnum = potnum * maxParameter;

            if (r <= mesh[potnum+5])
            {
                arij = r - mesh[potnum+5];
                vPair = (arij*arij) * (mesh[potnum] + mesh[potnum+1] * r + (mesh[potnum+2] * rsq) + (mesh[potnum+3] * pow(r, 3)) + (mesh[potnum+4] * pow(r, 4)));
        
                //std::cout << "\n pair " << r << " " << vPair << " " << (arij*arij) << " " << mesh[potnum] << " " 
                //<< mesh[potnum+1] * r << " " << mesh[potnum+5];
            }

            break;

        case gupt:
                    
            potnum = potnum * maxParameter;
            arij = (r - mesh[potnum+1]) / mesh[potnum+1];

            vPair = 2.0 * mesh[potnum] * exp(-mesh[potnum+2] * arij);
            break;

        case eam:
            if (r < metPots[potnum].potCutoff)
            {
                dr = metPots[potnum].stepRep;
                idx =  int((r - metPots[potnum].startRep) / dr);
                vPair = eamMeshEnergy(metPots[potnum].parRep , metPots[potnum].distRep , r, dr, idx) * 2.0;
            }
            break;
    }
}

double ManyBody::eamMeshEnergy(double* energyMesh, double* distMesh, double r, double dr, int index)
{
    //linear interpolation of the energy
    double energy = energyMesh[index] + ((energyMesh[index + 1] - energyMesh[index]) * (r - distMesh[index]) / dr);
    return energy;
}

double ManyBody::calculateManyBodyDensEnergy(int ltypea, int ltypeb, double r, double rsq)
{
    int
        pot,
        idx,
        potnum = -1;

    double
        vDens = 0.0,
        dr,
        arij,
        temp1;    
        
    if (r < 1.0e-8)
        return vDens;


    //screeening function
    if (potEAM)
    {
        pot = ltypeb * numSpec + ltypeb;  //density is type b in EAM model
    }
    else
    {
        pot = ltypea * numSpec + ltypeb;
    }

    potnum = meshTable[pot];

    if (potnum < 0)
        return vDens;

    //both potentials are the same here
    switch (meshType[potnum])
    {

        case sutchen:

            potnum = potnum * maxParameter;
            arij = mesh[potnum+1] / r;
            temp1 = mesh[potnum]*mesh[potnum+4];
            vDens = temp1 * temp1 * pow(arij, mesh[potnum+3]);
            break;

        case fnsc:

            potnum = potnum * maxParameter;
            
            if (r <= mesh[potnum+5])
            {
                arij = r - mesh[potnum+5];
                temp1 = mesh[potnum+4] * mesh[potnum+4];
                vDens = temp1 * (pow(arij, 2.0) + mesh[potnum+6] * pow(arij, 3.0) / mesh[potnum+5]);
            }
            break;

        case xfnsc:

            potnum = potnum * maxParameter;
            if (r <= mesh[potnum+7])
            {
                arij = r - mesh[potnum+7];
                temp1 = mesh[potnum+6] * mesh[potnum+6];
                vDens = temp1 * (pow(arij, 2.0) + pow(mesh[potnum+8],2.0) * pow(arij, 4.0));
                //std::cout << "\n dens " << r << " " << vDens;
            }
            break;

        case gupt:

            potnum = potnum * maxParameter;
            
            temp1 = mesh[potnum+3] * mesh[potnum+3];
            arij = (r - mesh[potnum+1]) / mesh[potnum+1];
            vDens = temp1 * exp(-2.0 * mesh[potnum+4] * arij);
            break;

        case eam:

            if (r < metPots[potnum].potCutoff)
            {
                dr = metPots[potnum].stepDens;
                idx =  int(r / dr); 
                vDens = eamMeshEnergy(metPots[potnum].parDens , metPots[potnum].distDens , r, dr, idx);
                //std::cout << "\n dens " << r << " " << dr << " " << idx  << " " << vDens;
                //exit(-1);
            }
            break;
    }

    return vDens;
}
double ManyBody::embed(int ltypea, int ltypeb, double rho)
{
    int
        pot,
        idx,
        potnum = -1;

    double
        dr,
        embed_rho = 0.0;

    if (potEAM == false)
    {
        embed_rho -= sqrt(rho);
        return embed_rho;
    }
        

    pot = ltypea * numSpec + ltypeb;

    potnum = meshTable[pot];

    if (potnum < 0)
        return embed_rho;

    switch (meshType[potnum])
    {
        
        case sutchen:
            embed_rho = -sqrt(rho);
            break;
        
        case gupt:
            embed_rho -= sqrt(rho);
            break;

        case fnsc:
            embed_rho = -sqrt(rho);
            break;

        case xfnsc:
            embed_rho = -sqrt(rho);
            break;
        
        case eam:
            dr = metPots[potnum].stepEmbed;
            idx =  int(rho / dr); 

            //make sure the index is within bounds
            if (idx < 0 || idx >= metPots[potnum].numPtsEmbed)
            {
                std::cout << "\n the rho for embedding is out of bounds " << idx << " " << rho << " " << dr << std::endl;
                std::cerr << "\n the rho for embedding is out of bounds " << idx << " " << rho << " " << dr << std::endl;
                exit(EXIT_FAILURE);
            }

            //std::cout << "\n emedding " << dr << " " << rho << " " << idx << " " << metPots[potnum].parEmbed[idx] 
            //                                                              << " " << metPots[potnum].distEmbed[idx] << std::endl;
            //embed_rho = metPots[potnum].parEmbed[idx] + ((metPots[potnum].parEmbed[idx + 1] - metPots[potnum].parEmbed[idx]) 
            //          * (rho - metPots[potnum].distEmbed[idx]) / dr);
            embed_rho = eamMeshEnergy(metPots[potnum].parEmbed , metPots[potnum].distEmbed , rho, dr, idx);
            break;
    }

    return embed_rho;    
}

void ManyBody::calculateManyBodyPairForce(int ltypea, int ltypeb, double r, double rsq, double& vPair, double& fPair)
{
    int
        pot,
        idx,
        potnum = -1;

    double
        dr,
        arij,
        temp1,
        temp2;

    vPair = fPair = 0.0;
        
    if (r < 1.0e-8)
        return;

        
    pot = ltypea * numSpec + ltypeb;

    potnum = meshTable[pot];

    if (potnum < 0)
        return;

    //screeening function
    
    //both potentials are the same here
    switch (meshType[potnum])
    {
        

        case sutchen:
            potnum = potnum * maxParameter;
            arij = mesh[potnum+1] / r;
            vPair = mesh[potnum] * pow(arij, mesh[potnum+2]);
            fPair = mesh[potnum] * mesh[potnum+2] * pow(arij, mesh[potnum+2]);
            break;

        case fnsc:
            potnum = potnum * maxParameter;

            if (r <= mesh[potnum+3])
            {
                arij = r - mesh[potnum+3];
                vPair = (mesh[potnum] + mesh[potnum+1] * r + mesh[potnum+2] * rsq) * pow(r - mesh[potnum+3], 2.0);
                temp1 = 2.0 * arij * (mesh[potnum] + mesh[potnum+1] * r + mesh[potnum+2] * rsq);
                temp2 = arij * arij * (mesh[potnum+1] + 2.0 * mesh[potnum+2] * r);
                fPair -= r * (temp1 + temp2);
            }
            break;

        case xfnsc:
            potnum = potnum * maxParameter;

            if (r <= mesh[potnum+5])
            {
                arij = r - mesh[potnum+5];
                vPair = (arij*arij) * (mesh[potnum] + mesh[potnum+1] * r + (mesh[potnum+2] * rsq) + (mesh[potnum+3] * pow(r, 3)) + (mesh[potnum+4] * pow(r, 4)));
        

                temp1 = 2.0 * arij * (mesh[potnum] + (mesh[potnum+1] * r) + (mesh[potnum+2] * rsq) + 
                        (mesh[potnum+3] * pow(r, 3)) + (mesh[potnum+4] * pow(r, 4)));
                temp2 = arij * arij * (mesh[potnum+1] + (2.0 * mesh[potnum+2] * r) + (3.0 * mesh[potnum+3] * pow(r, 2.0)) + 
                        (4.0 * mesh[potnum+4] * pow(r, 3.0)));
                fPair -= r * (temp1 + temp2);
            }
            break;

        case gupt:
                    
            potnum = potnum * maxParameter;
            arij = (r - mesh[potnum+1]) / mesh[potnum+1];

            vPair = mesh[potnum] * exp(-mesh[potnum+2] * arij);
            fPair = ((mesh[potnum] * mesh[potnum+2]) / mesh[potnum+1]) * exp(-mesh[potnum+2] * arij) * r;
            break;

        case eam:
            if (r < metPots[potnum].potCutoff)
            {
                dr = metPots[potnum].stepRep;
                idx =  int(r / dr);
                vPair = eamMeshEnergy(metPots[potnum].parRep , metPots[potnum].distRep , r, dr, idx);
                fPair -= r * eamMeshEnergy(metPots[potnum].derivRep , metPots[potnum].distRep , r, dr, idx);
            }
            break;
    }
}
void ManyBody::calculateManyBodyForce(int ltypea, int ltypeb, double r, double rsq, double rhoi, double rhoj, double& fDens)
{
    int
        idx,
        pot,
        potnum = -1;

    double
        dr,
        arij,
        fact,
        temp1,
        temp2,
        drho1,
        drho2,
        tiny = 1.0e-8;

    fDens = 0.0;
        
    if (r < tiny || rhoi < tiny || rhoj < tiny)
        return;

        
    pot = ltypea * numSpec + ltypeb;

    potnum = meshTable[pot];

    if (potnum < 0)
        return;

    //screeening function
    

    //both potentials are the same here
    switch (meshType[potnum])
    {
        

        case sutchen:
            potnum = potnum * maxParameter;
            arij = mesh[potnum+1] / r;
            fact = mesh[potnum]*mesh[potnum+4];
            temp1 = (fact / sqrt(rhoi)) + (fact / sqrt(rhoj));
            temp2 = 0.5 * mesh[potnum] * mesh[potnum+3] * mesh[potnum+4];
            //std::cout << "\n parameters " << mesh[potnum] << " " << mesh[potnum+3]  << " " <<  mesh[potnum+4] << " " << pow(arij, mesh[potnum+3]) << " " << temp1;
            fDens = temp1 * temp2 * pow(arij, mesh[potnum+3]);
            break;

        case fnsc:
            potnum = potnum * maxParameter;
            if (r < mesh[potnum+5])
            {
                arij = r - mesh[potnum+5];
                temp1 = mesh[potnum+4] / sqrt(rhoi) + mesh[potnum+4] / sqrt(rhoj);
                temp2 = 2.0 * arij + 3.0 * mesh[potnum+6] * (arij * arij / mesh[potnum+5]);
                fDens -= r * 0.5 * mesh[potnum+4] * temp1 * temp2;
            }
            break;

        case xfnsc:
            potnum = potnum * maxParameter;
            if (r < mesh[potnum+7])
            {
                arij = r - mesh[potnum+7];
                temp1 = mesh[potnum+6] / sqrt(rhoi) + mesh[potnum+6] / sqrt(rhoj);
                temp2 = 2.0 * arij + 4.0 * pow(mesh[potnum+8], 2.0) * pow(arij, 3.0);
                fDens -= r * 0.5 * mesh[potnum+6] * temp1 * temp2;
            }
            break;

        case gupt:
            potnum = potnum * maxParameter;
            arij = (r - mesh[potnum+1]) / mesh[potnum+1];
            fact = mesh[potnum+3];
            temp1 = (fact / sqrt(rhoi)) + (fact / sqrt(rhoj));
            temp2 = exp(-2.0 * mesh[potnum+4] * arij);
            fDens = r * ((mesh[potnum+3] * mesh[potnum+4]) / mesh[potnum+1]) * temp1 * temp2;
            break;

        case eam:
            if (r < metPots[potnum].potCutoff)
            {
                //the d_embed / d_rho for i and j
                dr = metPots[potnum].stepEmbed;
                idx =  int(rhoi / dr); 

                //make sure the index is within bounds
                if (idx < 0 || idx >= metPots[potnum].numPtsEmbed)
                {
                    std::cout << "\n the rho_i for embedding is out of bounds " << idx << " " << rhoi << std::endl;
                    std::cerr << "\n the rho_i for embedding is out of bounds " << idx << " " << rhoi << std::endl;
                    exit(EXIT_FAILURE);
                }
                drho1 = eamMeshEnergy(metPots[potnum].derivEmbed , metPots[potnum].distEmbed , rhoi, dr, idx);
                //std::cout << "\n d_embed 1 " << dr << " " << idx << " " << drho1;
                dr = metPots[potnum].stepEmbed;
                idx =  int(rhoj / dr); 

                //make sure the index is within bounds
                if (idx < 0 || idx >= metPots[potnum].numPtsEmbed)
                {
                    std::cout << "\n the rho_j for embedding is out of bounds " << idx << " " << rhoj << std::endl;
                    std::cerr << "\n the rho_j for embedding is out of bounds " << idx << " " << rhoj << std::endl;
                    exit(EXIT_FAILURE);
                }
                drho2 = eamMeshEnergy(metPots[potnum].derivEmbed , metPots[potnum].distEmbed , rhoj, dr, idx);
                //std::cout << "\n d_embed 2 " << dr << " " << idx << " " << drho2;
                //density term
                dr = metPots[potnum].stepDens;
                idx =  int(r / dr);
                temp1 = eamMeshEnergy(metPots[potnum].derivDens , metPots[potnum].distDens , r, dr, idx);
                //std::cout << "\n d_dens " << dr << " " << idx << " " << temp1;

                fDens = ((drho1 + drho2) * temp1) * r;
            }
            break;


    }

}
/*******************************************************************
 a mesh is calculated for the real space sum and during energy calculation
 this is interpolated
 *******************************************************************/
void ManyBody::calculateEnergyMesh(const Species& ele, double shortrangecut, std::ofstream& outStream)

{
    int
        pot = 999,
        num = ele.getNumSpecies();

    Element
        elei,
        elej;

    MetalPotential
        poten;

    bool
        founda,
        foundb;

    // set up mesh
    // first dimension arrays
    numSpec = num;

    outStream << "\n calculating energy mesh " << std::endl;

    if (potEAM == true)
    {
        calculateEAMEnergyMesh(ele, shortrangecut, outStream);
    }
    else
    {
        /* maxInteraction = numSpec*NumSpec. This allows for a table to be set up incase i > or j > i which takes a bit of memory
        but hopefully speeds up the calculation of the forces */
        maxInteraction = ele.getMaxNumInteractions();
        meshTable = alloc_ivector(maxInteraction, "many body table", -1);
        meshType = alloc_ivector(maxInteraction, "many body table", -1);
        mesh = alloc_dvector(maxInteraction * maxParameter, "two body table", 0.0);


        cutoff = shortrangecut;

        for (int i = 0; i < num; i++)
        {

            elei = ele.getSpecies(i);

            for (int j = 0; j < num; j++)
            {

                elej = ele.getSpecies(j);

                for (int k = 0; k < numberPotentials; k++)
                {

                    poten =metPots[k];

                    founda = foundb = false;
                    
                    if (elei.name == poten.atA) 
                        founda = true;

                    if (elej.name == poten.atB) 
                        foundb = true;

                    if (founda && foundb) 
                    {

                        pot = j * numSpec + i;
                        meshTable[pot] = k;
                        pot = i * numSpec + j;
                        meshTable[pot] = k;   // symmetrical part
                        outStream << "\n potential found between " << elei.name << "  and  " << elej.name;
                        //outStream << "\n poten " << k << " " << i << " " << j << " " << pot  << " " << k * maxParameter;

                        //copy potential parameters over to the work array
                        meshType[k] = poten.potentialType;
                        mesh[maxParameter*k+0] = poten.param1;
                        mesh[maxParameter*k+1] = poten.param2;
                        mesh[maxParameter*k+2] = poten.param3;
                        mesh[maxParameter*k+3] = poten.param4;
                        mesh[maxParameter*k+4] = poten.param5;
                        mesh[maxParameter*k+5] = poten.param6;
                        mesh[maxParameter*k+6] = poten.param7;
                        mesh[maxParameter*k+7] = poten.param8;
                        mesh[maxParameter*k+8] = poten.param9;
                        //outStream << "\n mesh " << mesh[maxParameter*k+0] << " " << mesh[maxParameter*k+1] << " " << mesh[maxParameter*k+2];

                    } // end if for found both
                
                } // end loop over potentials

                //if (!founda || !foundb)
                //{
                //    outStream << "\n WARNING: potential not found between " << elei.name << "  and  " << elej.name;
                //}
                
            } // end loop over species j

        } // end loop over species i
    }

}
void ManyBody::calculateEAMEnergyMesh(const Species& ele, double shortrangecut, std::ofstream& outstream)

{
    int
        numSpec = ele.getNumSpecies(),
        pot;

    Element
        elei,
        elej;

    bool
        founda,
        foundb;
    
    // dimension arrays

    int maxInteraction = ele.getMaxNumInteractions();
    meshTable = alloc_ivector(maxInteraction, "many-body table", -1);
    meshType = alloc_ivector(maxInteraction, "many body table", 0);
    mesh = alloc_dvector_aligned(maxInteraction * maxParameter, "two body table", 0.0);
    
    for (int k = 0; k < numberPotentials; k++)
    {
        //poten = metPots[k];

        for (int i = 0; i < numSpec; i++)
        {

            elei = ele.getSpecies(i);

            for (int j = 0; j < numSpec; j++)
            {

                elej = ele.getSpecies(j);

                founda = foundb = false;

                if (elei.name == metPots[k].atA)
                    founda = true;

                if (elej.name == metPots[k].atB) 
                    foundb = true;

                if (founda && foundb) 
                {

                    pot = j * numSpec + i;
                    meshTable[pot] = k;
                    pot = i * numSpec + j;
                    meshTable[pot] = k;   // symmetrical part
                    meshType[k] = metPots[k].potentialType;
                    outstream << "\n potential between " << elei.name << "  and  " << elej.name << " " << k << std::endl;
                    //outstream << "\n poten " << k << " " << i << " " << j << " " << pot  << " " << k * maxParameter;

                } // end if for found both

            } // end loop over species j

        } // end loop over species i

    } // end loop over potentials

}

void ManyBody::checkMesh(void)
{
    std::cout << "\n vpair " << metPots[0].numPtsRep;
    for (int jj = 0; jj < metPots[0].numPtsRep; jj++)
    {
        //std::cout << "\n v-dens " << jj << " " << metPots[0].distDens[jj] << " " << metPots[0].parDens[jj];
        std::cout << "\n v-pair " << jj << " " << metPots[0].distRep[jj] << " " << metPots[0].parRep[jj];
    }
}
int ManyBody::potenKey(std::string keyWord)
{
    int
        ident = 0;

    if (keyWord == "suttonchen")
        ident = sutchen;

    else if (keyWord == "finnis")
        ident = fnsc;

    else if (keyWord == "xfnsc")
        ident = xfnsc;

    else if (keyWord == "gupt")
        ident = gupt;

    else if (keyWord == "eam")
        ident = eam;

    else if (keyWord == "meam")
        ident = meam;

    else if (keyWord == "end")
        ident = metalend;

    return ident;

}
