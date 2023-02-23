/**************************************************************************
 *
 * File name: omp2.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 11.4.2009
 *
 *************************************************************************/

#include "ompcore.h"
#include "omputils.h"


/* Input Arguments */

#define	IN_D	        prhs[0]
#define IN_X          prhs[1]
#define IN_DtX        prhs[2]
#define IN_XtX        prhs[3]
#define IN_G          prhs[4]
#define IN_EPS        prhs[5]
#define IN_SPARSE_G   prhs[6]
#define IN_MSGDELTA   prhs[7]
#define IN_MAXATOMS   prhs[8]
#define IN_PROFILE    prhs[9]


/* Output Arguments */

#define	GAMMA_OUT     plhs[0]


/***************************************************************************************/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

{
  double *D, *x, *DtX, *XtX, *G, eps, msgdelta;
  int gmode, maxatoms, profile;
  mwSize m, n, L;    /* D is n x m , X is n x L, DtX is m x L */

  
  /* check parameters */
  
  checkmatrix(IN_D, "D");
  checkmatrix(IN_X, "X");
  checkmatrix(IN_DtX, "DtX");
  checkmatrix(IN_XtX, "XtX");
  checkmatrix(IN_G, "G");
  
  checkscalar(IN_EPS, "EPSILON");
  checkscalar(IN_SPARSE_G, "sparse_g");
  checkscalar(IN_MSGDELTA, "msgdelta");
  checkscalar(IN_MAXATOMS, "maxatoms");
  checkscalar(IN_PROFILE, "profile");
  
  
  /* get parameters */
  
  x = D = DtX = XtX = G = 0;
  
  if (!mxIsEmpty(IN_D))
    D = mxGetPr(IN_D);
  
  if (!mxIsEmpty(IN_X))
    x = mxGetPr(IN_X);
  
  if (!mxIsEmpty(IN_DtX))
    DtX = mxGetPr(IN_DtX);
  
  if (!mxIsEmpty(IN_XtX))
    XtX = mxGetPr(IN_XtX);
  
  if (!mxIsEmpty(IN_G))
    G = mxGetPr(IN_G);
  
  eps = mxGetScalar(IN_EPS);
  if ((int)(mxGetScalar(IN_SPARSE_G)+1e-2)) {
    gmode = SPARSE_GAMMA;
  }
  else {
    gmode = FULL_GAMMA;
  }
  msgdelta = mxGetScalar(IN_MSGDELTA);
  if (mxGetScalar(IN_MAXATOMS) < -1e-5) {
    maxatoms = -1;
  }
  else {
    maxatoms = (int)(mxGetScalar(IN_MAXATOMS)+1e-2);
  }
  profile = (int)(mxGetScalar(IN_PROFILE)+1e-2);
  
  
  /* check sizes */
  
  if (D && x) {
    n = mxGetM(IN_D);
    m = mxGetN(IN_D);
    L = mxGetN(IN_X);
    
    if (mxGetM(IN_X) != n) {
      mexErrMsgTxt("D and X have incompatible sizes.");
    }
    
    if (G) {
      if (mxGetN(IN_G)!=mxGetM(IN_G)) {
        mexErrMsgTxt("G must be a square matrix.");
      }
      if (mxGetN(IN_G) != m) {
        mexErrMsgTxt("D and G have incompatible sizes.");
      }
    }
  }
  
  else if (DtX && XtX) {
    m = mxGetM(IN_DtX);
    L = mxGetN(IN_DtX);
    
    /* set n to an arbitrary value that is at least the max possible number of selected atoms */
    
    if (maxatoms>0) {
      n = maxatoms;
    }
    else {
      n = m;
    }
    
    if ( !(mxGetM(IN_XtX)==L && mxGetN(IN_XtX)==1) && !(mxGetM(IN_XtX)==1 && mxGetN(IN_XtX)==L) ) {
      mexErrMsgTxt("DtX and XtX have incompatible sizes.");
    }
    
    if (mxGetN(IN_G)!=mxGetM(IN_G)) {
      mexErrMsgTxt("G must be a square matrix.");
    }
    if (mxGetN(IN_G) != m) {
      mexErrMsgTxt("DtX and G have incompatible sizes.");
    }
  }
  
  else {
    mexErrMsgTxt("Either D and X, or DtX and XtX, must be specified.");
  }
  
  
  /* Do OMP! */
  
  GAMMA_OUT = ompcore(D, x, DtX, XtX, G, n, m, L, maxatoms, eps, gmode, profile, msgdelta, 1);
  
  return;
}
