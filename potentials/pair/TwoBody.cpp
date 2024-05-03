#include "TwoBody.h"

std::vector<std::string> split(std::string s);

//////////////////////////////////////////////////////////////////////
// construction/destruction
//////////////////////////////////////////////////////////////////////

TwoBody::TwoBody()
{

    numberPotentials = 0;  // number of potential parameters  

}

TwoBody::~TwoBody()
{
    free_ivector(meshType);          
    free_ivector(meshTable);                 

    free_dvector_aligned(mesh);
}

TwoBody& TwoBody::operator=(const TwoBody& src)
{
    int
        i;

    if (this == &src) //ie self assignment
    {
        return *this;
    }

    numberPotentials = src.numberPotentials;
    numSpec = src.numSpec;
    maxInteraction = src.maxInteraction;

    meshTable = alloc_ivector(maxInteraction, "two body table", 0);
    meshType = alloc_ivector(maxInteraction, "two body table", 0);
    mesh = alloc_dvector(maxInteraction * MAXPAR, "two body table", 0.0);

    for (i = 0; i < maxInteraction; i++)
    {
        meshTable[i] = src.meshTable[i];
        meshType[i] = src.meshType[i];
    }

    for (i = 0; i < maxInteraction * MAXPAR; i++)
        mesh[i] = src.mesh[i];
    
    return *this;
}

/* member function to dimension and load short-range potentail parameters
 */
void TwoBody::loadPotential(std::ifstream& inStream, std::ofstream& outStream, int num)
{
    int
        keyval;

    std::string
        finished = " ",
        line,
        dummy,
        keyWord;            // directive read as temp

    Potential
        tpot;

    numberPotentials = num;

    vdwPots = std::vector<Potential>(numberPotentials);  // dimension the array
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

        case buckingham:                         // buckingham potential

            tpot.potentialType = keyval;

            std::getline(inStream, line);
	        words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

            dummy = words[2];
            tpot.param1 = std::stod(dummy);
            dummy = words[3];
            tpot.param2 = std::stod(dummy);
            dummy = words[4];
            tpot.param3 = std::stod(dummy);

            break;

        case morse:                            // morse potential

            tpot.potentialType = keyval;

            std::getline(inStream, line);
	        words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

            dummy = words[2];
            tpot.param1 = std::stod(dummy);
            dummy = words[3];
            tpot.param2 = std::stod(dummy);
            dummy = words[4];
            tpot.param3 = std::stod(dummy);
            dummy = words[5];
            tpot.param4 = std::stod(dummy);
            
            break;

        case ljones:     // len-jones parameters                    

            tpot.potentialType = keyval;

            std::getline(inStream, line);
		    words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

            dummy = words[2];
            tpot.param1 = std::stod(dummy);
            dummy = words[3];
            tpot.param2 = std::stod(dummy);
            
            break;

        case bhm:                            // born huggins mayer

            tpot.potentialType = keyval;

            std::getline(inStream, line);
		    words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

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

        case wca:

            tpot.potentialType = keyval;

            std::getline(inStream, line);
		    words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

            dummy = words[2];  // eta
            tpot.param1 = std::stod(dummy);  // sigma
            dummy = words[3];
            tpot.param2 = std::stod(dummy);
            
            break;

        case ew:     // espanol and warren potential

            tpot.potentialType = keyval;
            
            std::getline(inStream, line);
		    words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

            dummy = words[2];  // A  
            tpot.param1 = std::stod(dummy);

            break;

        case bond:     // espanol and warren potential

            tpot.potentialType = keyval;
            
            std::getline(inStream, line);
		    words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

            dummy = words[2];  // A  
            tpot.param1 = std::stod(dummy);

            break;

        case sevnsix:   // severn six

            tpot.potentialType = keyval;
            
            std::getline(inStream, line);
		    words = split(line);

            tpot.iona = words[0];
            tpot.ionb = words[1];

            dummy = words[2];  // eta
            tpot.param1 = std::stod(dummy);  // sigma
            dummy = words[3];
            tpot.param2 = std::stod(dummy);

            break;

        case potends:

            break;


        default:
            std::cerr << "\n\n*** program error : potential type not recognised ***";
            std::cerr.flush();

        }                  // end of case statements

        vdwPots[pot] = tpot;


    }


}

/********************************************************
 convert the potentials to internal energy units
 ********************************************************/
void TwoBody::convertPotential(double unit)
{
    Potential
        tpot;        // temp storage for a potential


    for (int pot = 0; pot < numberPotentials; pot++)
    {
        tpot = vdwPots[pot];




        if (tpot.potentialType == buckingham)
        {
            tpot.param1 *= unit;
            tpot.param3 *= unit;

        }

        if (tpot.potentialType == morse)
        {
            tpot.param1 *= unit;
        }

        if (tpot.potentialType == ljones)
        {
            tpot.param1 *= unit;
            //tpot.param2 *= unit;
        }

        if (tpot.potentialType == bhm)
        {
            tpot.param1 *= unit;
            tpot.param4 *= unit;
            tpot.param5 *= unit;
        }

        if (tpot.potentialType == wca)     // weeks-chandler-anderson potential
        {
            tpot.param1 *= unit;
        }

        if (tpot.potentialType == ew)     // espanol and warren potential
        {
            tpot.param1 *= unit;
        }

        if (tpot.potentialType == bond)     // espanol and warren potential
        {
            tpot.param1 *= unit;
        }

        if (tpot.potentialType == sevnsix)
        {
            tpot.param1 *= unit;
        }

        vdwPots[pot] = tpot;  // change stored value

    }

}

/********************************************************
 print out all the potentials to a file
 ********************************************************/
void TwoBody::printPotential(std::ofstream& outStream)
{
    Potential
        tpot;        // temp storage for a potential


   
    outStream << "\n \n *********************************************************************" << std::endl;
    outStream << "            two-body potential data";
    outStream << "\n \n *********************************************************************" << std::endl;
    for (int pot = 0; pot < numberPotentials; pot++)
    {
        tpot = vdwPots[pot];


        outStream << "\n\n      potential          ion         ion ";


        outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
            << "\n"
            << std::setw(15) << tpot.potentialType
            << std::setw(13) << tpot.iona  << std::setw(13) << tpot.ionb << std::endl;

        if (tpot.potentialType == buckingham)
        {
            outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                << std::setprecision(5) << "\n"
                << std::setw(6) << "a"
                << std::setw(15) << tpot.param1
                << std::setw(6) << "rho"
                << std::setw(15) << tpot.param2
                << std::setw(6) << "c"
                << std::setw(15) << tpot.param3 << std::endl;
        }

        if (tpot.potentialType == morse)
        {
            outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                << std::setprecision(5) << "\n"
                << std::setw(6) << "d"
                << std::setw(15) << tpot.param1
                << std::setw(6) << "beta"
                << std::setw(15) << tpot.param2
                << std::setw(6) << "req"
                << std::setw(15) << tpot.param3 << std::endl;
        }

        if (tpot.potentialType == wca)
        {
            outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                << std::setprecision(5) << "\n"
                << std::setw(6) << "eta"
                << std::setw(15) << tpot.param1
                << std::setw(6) << "sigma"
                << std::setw(15) << tpot.param2 << std::endl;
        }

        if (tpot.potentialType == bhm)
        {
            outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                << std::setprecision(5) << "\n"
                << std::setw(6) << "a"
                << std::setw(15) << tpot.param1
                << std::setw(6) << "rho"
                << std::setw(15) << tpot.param2
                << std::setw(6) << "req"
                << std::setw(15) << tpot.param3
                << std::setw(6) << "ddd"
                << std::setw(15) << tpot.param4
                << std::setw(6) << "eee"
                << std::setw(15) << tpot.param5 << std::endl;
        }

        if (tpot.potentialType == bond)
        {
            outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                << std::setprecision(5) << "\n"
                << std::setw(6) << "bond energy"
                << std::setw(15) << tpot.param1 << std::endl;
        }

        if (tpot.potentialType == bond)
        {
            outStream << std::setiosflags(std::ios::showpoint | std::ios::fixed)
                << std::setprecision(5) << "\n"
                << std::setw(6) << "eta"
                << std::setw(15) << tpot.param1
                << std::setw(6) << "sigma"
                << std::setw(15) << tpot.param2 << std::endl;
        }

    }

}



/* returns an integer value that corresponds to a given potential keyWord
 * if the keyWord is not found a value of zero is returned
 */
int TwoBody::potenKey(std::string keyWord)
{
    int
        ident = 0;

    if (keyWord == "buck")
        ident = buckingham;

    else if (keyWord == "morse")
        ident = morse;

    else if (keyWord == "ljones")
        ident = ljones;

    else if (keyWord == "bhm")
        ident = bhm;

    else if (keyWord == "ew")
        ident = ew;

    else if (keyWord == "wca")
        ident = wca;

    else if (keyWord == "bond")
        ident = bond;

    else if (keyWord == "severnsix")
        ident = sevnsix;

    else if (keyWord == "end")
        ident = potends;

    return ident;

}

/*******************************************************************
 a mesh is calculated for the real space sum and during energy calculation
 this is interpolated
 
 *******************************************************************/

void TwoBody::calculateEnergyMesh(const Species& ele, double shortrangecut, std::ofstream& outStream)

{
    int
        pot = 999,
        num = ele.getNumSpecies();
        

    double
        radius;

    Element
        elei,
        elej;

    Potential
        poten;

    bool
        founda,
        foundb;

    // set up mesh
    // first dimension arrays
    numSpec = num;

    /* maxInteraction = numSpec*NumSpec. This allows for a table to be set up incase i > or j > i which takes a bit of memory
    but hopefully speeds up the calculation of the forces */
    maxInteraction = ele.getMaxNumInteractions();
    meshTable = alloc_ivector(maxInteraction, "two body table", 0);
    meshType = alloc_ivector(maxInteraction, "two body table", 0);
    mesh = alloc_dvector_aligned(maxInteraction * MAXPAR, "two body table", 0.0);


    cutoff = shortrangecut;
    radius = shortrangecut;

    delMesh = radius / (MAXMESH - 4);      // spacing between mesh points


    for (int i = 0; i < num; i++)
    {

        elei = ele.getSpecies(i);

        for (int j = 0; j < num; j++)
        {

            elej = ele.getSpecies(j);

            for (int k = 0; k < numberPotentials; k++)
            {

                poten =vdwPots[k];

                founda = foundb = false;

                if (elei.name == poten.iona) 
                    founda = true;

                if (elej.name == poten.ionb) 
                    foundb = true;

                if (founda && foundb) 
                {

                    pot = j * numSpec + i;
                    meshTable[pot] = k;
                    pot = i * numSpec + j;
                    meshTable[pot] = k;   // symmetrical part
                    outStream << "\n potential found between " << elei.name << "  and  " << elej.name;
                    //outStream << "\n poten " << k << " " << i << " " << j << " " << pot  << " " << k * MAXPAR;

                    //copy potential parameters over to the work array
                    meshType[k] = poten.potentialType;
                    mesh[MAXPAR*k+0] = poten.param1;
                    mesh[MAXPAR*k+1] = poten.param2;
                    mesh[MAXPAR*k+2] = poten.param3;
                    mesh[MAXPAR*k+3] = poten.param4;
                    mesh[MAXPAR*k+4] = poten.param5;
                    mesh[MAXPAR*k+5] = poten.param6;
                    //outStream << "\n mesh " << mesh[MAXPAR*k+0] << " " << mesh[MAXPAR*k+1] << " " << mesh[MAXPAR*k+2];

                } // end if for found both
                
            } // end loop over potentials

            //if (!founda || !foundb)
            //{
            //    outStream << "\n WARNING: potential not found between " << elei.name << "  and  " << elej.name;
            //}
                
        } // end loop over species j

    } // end loop over species i

}
