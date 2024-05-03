#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>

std::vector<std::string> split(std::string s)
{
    std::vector<std::string> words;

    std::stringstream ss(s);
    std::string word;
    
    while (ss >> word) 
    {
        words.push_back(word);
    }

    return words;
}

void dcell(double *aaa, double bbb[10])

/*
!
! dl_poly_3 subroutine to calculate the dimensional properies of a
! simulation cell specified by the input 3x3 matrix aaa [cell vectors in
! rows, the matrix is in the form of one dimensional reading
! [row1,row2,row3].
!
! The results are returned in the array bbb, with:
!
! bbb[1 to 3] - lengths of cell vectors
! bbb[4 to 6] - cosines of cell angles
! bbb[7 to 9] - perpendicular cell widths
! bbb[10]     - cell volume
!
! copyright - daresbury laboratory
! author    - w.smith july 1992
! amended   - i.t.todorov may 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
{
  double
      axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3, 
      x[3],y[3],z[3],d[3];

//    calculate lengths of cell vectors

      bbb[0]=sqrt(aaa[0]*aaa[0]+aaa[1]*aaa[1]+aaa[2]*aaa[2]);
      bbb[1]=sqrt(aaa[3]*aaa[3]+aaa[4]*aaa[4]+aaa[5]*aaa[5]);
      bbb[2]=sqrt(aaa[6]*aaa[6]+aaa[7]*aaa[7]+aaa[8]*aaa[8]);

//    calculate cosines of cell angles

      bbb[3]=(aaa[0]*aaa[3]+aaa[1]*aaa[4]+aaa[2]*aaa[5])/(bbb[0]*bbb[1]);
      bbb[4]=(aaa[0]*aaa[6]+aaa[1]*aaa[7]+aaa[2]*aaa[8])/(bbb[0]*bbb[2]);
      bbb[5]=(aaa[3]*aaa[6]+aaa[4]*aaa[7]+aaa[5]*aaa[8])/(bbb[1]*bbb[2]);

//    calculate vector products of cell vectors

      axb1=aaa[1]*aaa[5]-aaa[2]*aaa[4];
      axb2=aaa[2]*aaa[3]-aaa[0]*aaa[5];
      axb3=aaa[0]*aaa[4]-aaa[1]*aaa[3];
      bxc1=aaa[4]*aaa[8]-aaa[5]*aaa[7];
      bxc2=aaa[5]*aaa[6]-aaa[3]*aaa[8];
      bxc3=aaa[3]*aaa[7]-aaa[4]*aaa[6];
      cxa1=aaa[7]*aaa[2]-aaa[1]*aaa[8];
      cxa2=aaa[0]*aaa[8]-aaa[2]*aaa[6];
      cxa3=aaa[1]*aaa[6]-aaa[0]*aaa[7];

//    calculate volume of cell

      bbb[9]=abs(aaa[0]*bxc1+aaa[1]*bxc2+aaa[2]*bxc3);

//    calculate cell perpendicular widths

      bbb[6]=bbb[9]/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3);
      bbb[7]=bbb[9]/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3);
      bbb[8]=bbb[9]/sqrt(axb1*axb1+axb2*axb2+axb3*axb3);

}
