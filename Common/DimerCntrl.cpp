#include "DimerCntrl.h"

std::vector<std::string> split(std::string s);


void DimerCntrl::readDimer(std::ifstream& inStream, std::ofstream& outStream)
{
    std::string
        dummy,
        keyWord,
        subWord,
        line;

    bool
        readMore = true;

    std::vector<std::string> words;

	//start reading the job control parameters
	while (readMore)
	{
		std::getline(inStream, line);

        if (line[0] == '#')
            continue;

        transform(line.begin(), line.end(), line.begin(), ::tolower);

        words = split(line);

        if (words.size() == 0)
            continue;

        keyWord = words[0];

        if (keyWord == "}")
        {
            readMore = false;
            break;
        }
        else if (keyWord == "target")
        {
            subWord = words[1];

            if (subWord == "global")
            {
                regionStyle = 1;
            }
            else if (subWord == "local")
            {
                regionStyle = 2;
                dummy = words[2];
                regionRadius = std::stod(dummy);
            }
            else if (subWord == "globaltypes")  // all the atoms of these types are displaced
            {
                regionStyle = 3;

                dummy = words[2];
                int num = std::stoi(dummy);
                for (int i = 0; i < num; i++)
                {
                    std::getline(inStream, line);
					words = split(line);
                    dummy = words[0];
                    dimerTypes1.push_back(dummy);
                }
            }
            else if (subWord == "localtypes")  // target atom of type is chosen and only atoms within a given radius are ddisplaced
            {
                regionStyle = 4;

                dummy = words[2];
                int num = std::stoi(dummy);

                dummy = words[3];
                regionRadius = std::stod(dummy);

                for (int i = 0; i < num; i++)
                {
                    std::getline(inStream, line);
					words = split(line);
                    dummy = words[0];
                    dimerTypes1.push_back(dummy);
                }
            }
            else if (subWord == "coordination")  // target atom of type is chosen and only atoms within a given radius are ddisplaced
            {
                regionStyle = 5;

                dummy = words[2];
                int num = std::stoi(dummy);

                dummy = words[3];
                regionRadius = std::stod(dummy);

                for (int i = 0; i < num; i++)
                {
                    std::getline(inStream, line);
					words = split(line);
                    dummy = words[0];
                    dimerTypes1.push_back(dummy);
                    dummy = words[1];
                    dimerTypes2.push_back(dummy);
                    dummy = words[2];
                    coordDistance.push_back(std::stod(dummy));
                    dummy = words[3];
                    idealCoordNumber.push_back(std::stoi(dummy));
                }
            }
        }
        else if (keyWord == "maxrotations")
        {
            dummy = words[1];
            maxRotations = std::stoi(dummy);
        }
        else if (keyWord == "dimerangle")
        {
            dummy = words[1];
            deltaTheta = std::stod(dummy);
        }
        else if (keyWord == "dimerdisplacement")
        {
            dummy = words[1];
            dimerDisplacement = std::stod(dummy);
        }
        else if (keyWord == "linedisplacement")
        {
            dummy = words[1];
            lineDisplacement = std::stod(dummy);
        }
        else if (keyWord == "maxsteplength")
        {
            dummy = words[1];
            maxDeltaX = std::stod(dummy);
        }
        else if (keyWord == "rotforcetolerance")
        {
            dummy = words[1];
            forceTol = std::stod(dummy);
        }
		else
		{
            outStream << "\n DIMER: keyword not found " << words[0];
			outStream.flush();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
		}
    }
}
								
