#include "NbrListPBC.h"
NbrListPBC::NbrListPBC()
{
    maxNbrs = 0;

    numNbrs = nullptr;

    map = nullptr;

    rZeroX = nullptr;
    rZeroY = nullptr;
    rZeroZ = nullptr;
}

NbrListPBC::~NbrListPBC()
{

    free_dvector(rZeroX);
    free_dvector(rZeroY);
    free_dvector(rZeroZ);
}

void NbrListPBC::allocNbrListArrays(int natoms, int num)
{

    maxNbrs = num;

    numNbrs = alloc_ivector(natoms, "nbrlist allocation", 0);
    map = alloc_ivector(natoms * num, "nbrlist allocation", 0);

    rZeroX = alloc_dvector(natoms,  "nbrlist allocation", 0.0);
    rZeroY = alloc_dvector(natoms,  "nbrlist allocation", 0.0);
    rZeroZ = alloc_dvector(natoms,  "nbrlist allocation", 0.0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calculation of neighbour lists
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void NbrListPBC::calculateNbrList(double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int* frozen,
                                  int natoms, double cutoff)
{
    int
        aatom,
        batom,
        idx,
        num;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        rx, ry, rz,
        radius, rsq;

    radius = cutoff * cutoff;
    
    //for (aatom = rank; aatom < natoms - 1; aatom = aatom + numProcs)
    for (aatom = 0; aatom < natoms; aatom++)
    {

        num = 0;
        idx = maxNbrs * aatom;
        //std::cout << "\n atom nbrlist " << aatom;
        //set new start for position check accummulators
        rZeroX[aatom] = 0.0;
        rZeroY[aatom] = 0.0;
        rZeroZ[aatom] = 0.0;

        ax = posX[aatom];
        ay = posY[aatom];
        az = posZ[aatom];

        for (batom = 0; batom < natoms; batom++)
        {

            if (aatom == batom)
                continue;

            bx = posX[batom];
            by = posY[batom];
            bz = posZ[batom];

            rx = ax - bx;
            ry = ay - by;
            rz = az - bz;

            xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
            yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
            zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

            xx -= rint(xx);
            yy -= rint(yy);
            zz -= rint(zz);

            rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
            ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
            rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

            rsq = rx * rx + ry * ry + rz * rz;

            if (rsq <= radius)
            {

                map[idx+num] = batom;
                num++;

                if (num >= maxNbrs)
                {
                    std::cout << "\n exceeded max nbrs " << maxNbrs << " " << num << std::endl;
                    std::cout.flush();
                    exit(EXIT_FAILURE);

                }

            }

        }

        numNbrs[aatom] = num;
        //std::cout << "\n PBC nbrlist " << aatom << " " << num;

    }

}

void NbrListPBC::calculateAtomNbrList(int aatom, double* posX, double* posY, double* posZ, double* latVector, double* rcpVector, int* frozen,
                                  int natoms, double cutoff)
{
    int
        batom,
        idx,
        num;

    double
        ax, ay, az,
        bx, by, bz,
        xx, yy, zz,
        rx, ry, rz,
        radius, rsq;

    radius = cutoff * cutoff;
    
    num = 0;
    idx = maxNbrs * aatom;
    //std::cout << "\n atom nbrlist " << aatom;
    //set new start for position check accummulators
    rZeroX[aatom] = 0.0;
    rZeroY[aatom] = 0.0;
    rZeroZ[aatom] = 0.0;

    ax = posX[aatom];
    ay = posY[aatom];
    az = posZ[aatom];

    for (batom = 0; batom < natoms; batom++)
    {

        if (aatom == batom)
            continue;

        bx = posX[batom];
        by = posY[batom];
        bz = posZ[batom];

        rx = ax - bx;
        ry = ay - by;
        rz = az - bz;

        xx = rx * rcpVector[0] + ry * rcpVector[3] + rz * rcpVector[6];
        yy = rx * rcpVector[1] + ry * rcpVector[4] + rz * rcpVector[7];
        zz = rx * rcpVector[2] + ry * rcpVector[5] + rz * rcpVector[8];

        xx -= rint(xx);
        yy -= rint(yy);
        zz -= rint(zz);

        rx = xx * latVector[0] + yy * latVector[3] + zz * latVector[6];
        ry = xx * latVector[1] + yy * latVector[4] + zz * latVector[7];
        rz = xx * latVector[2] + yy * latVector[5] + zz * latVector[8];

        rsq = rx * rx + ry * ry + rz * rz;

        if (rsq <= radius)
        {

            map[idx+num] = batom;
            num++;

            if (num >= maxNbrs)
            {
                std::cout << "\n exceeded max nbrs " << maxNbrs << " " << num << std::endl;
                std::cout.flush();
                exit(EXIT_FAILURE);

            }

        }

    }

    numNbrs[aatom] = num;
    //std::cout << "\n PBC nbrlist " << aatom << " " << num;

}

bool NbrListPBC::checkNbrList(int i, double deltaX, double deltaY, double deltaZ, double verletShell)
{
    double
        dr,
        rMax = pow(verletShell, 2.0);

    bool
        check = false;

    rZeroX[i] += deltaX;
    rZeroY[i] += deltaY;
    rZeroZ[i] += deltaZ;

    dr = rZeroX[i] * rZeroX[i] + rZeroY[i] * rZeroY[i] + rZeroZ[i] * rZeroZ[i];

    if (dr > rMax)
        check = true;
    
    return check;
}
