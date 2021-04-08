/**************************************************************************
 *
 * File name: omputils.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 11.4.2009
 *
 *************************************************************************/

#include "omputils.h"
#include <math.h>


const char FULL_GAMMA_STR[] = "full";
const char SPARSE_GAMMA_STR[] = "sparse";


/* verify that the mxArray contains a double matrix */

void checkmatrix(const mxArray *param, char *name)
{
  char errmsg[100];
  sprintf(errmsg, "OMP requires that %.20s be a double matrix.", name);
  if (!mxIsDouble(param) || mxIsComplex(param) || mxGetNumberOfDimensions(param)>2) {
    mexErrMsgTxt(errmsg);
  }
}


/* verify that the mxArray contains a double scalar */

void checkscalar(const mxArray *param, char *name)
{
  char errmsg[100];
  sprintf(errmsg, "OMP requires that %.20s be a double scalar.", name);
  if (!mxIsDouble(param) || mxIsComplex(param) || mxGetNumberOfDimensions(param)>2 || 
      mxGetM(param)!=1 || mxGetN(param)!=1) 
  {
    mexErrMsgTxt(errmsg);
  }
}


/* convert seconds to hours, minutes and seconds */

void secs2hms(double sectot, int *hrs, int *mins, double *secs)
{
  *hrs = (int)(floor(sectot/3600)+1e-2);
  sectot = sectot - 3600*(*hrs);
  *mins = (int)(floor(sectot/60)+1e-2);
  *secs = sectot - 60*(*mins);
}

