/**************************************************************************
 *
 * File name: sprow.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 24.8.2009
 *
 *************************************************************************/


#include "mex.h"
#include "mexutils.h"


/* Input Arguments */

#define	A_IN	prhs[0]
#define J_IN  prhs[1]


/* Output Arguments */

#define	X_OUT	  plhs[0]
#define	ID_OUT	plhs[1]


void mexFunction(int nlhs, mxArray *plhs[], 
		             int nrhs, const mxArray*prhs[])
     
{ 
  double *pr, *x, *id, rowid;
  mwIndex *ir, *jc;
  mwSize m, n;
  mwIndex i, j, k, l, rowlen;
  
  if (nrhs != 2) {
    mexErrMsgTxt("GETSPROW requires two input arguments.");
  } else if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  checkmatrix(A_IN, "GETSPROW", "A");
  checksparse(A_IN, "GETSPROW", "A");
  checkscalar(J_IN, "GETSPROW", "J");
  
  m = mxGetM(A_IN);
  n = mxGetN(A_IN);
  
  rowid = mxGetScalar(J_IN);
  if (rowid < 0) {
    mexErrMsgTxt("Invalid row index.");
  }
  j = (mwIndex)(rowid + 1e-2);
  if (j<1 || j>m) {
    mexErrMsgTxt("Row index is out of range.");
  }
  j--;
  
  pr = mxGetPr(A_IN);
  ir = mxGetIr(A_IN);
  jc = mxGetJc(A_IN);
  
  /* Determine length of row */
  rowlen = 0;
  for (i=0; i<jc[n]; ++i) {
    rowlen += (ir[i]==j);
  }
  
  /* Allocate output parameters */
  X_OUT = mxCreateDoubleMatrix(1, rowlen, mxREAL);
  ID_OUT = mxCreateDoubleMatrix(1, rowlen, mxREAL);
  
  x = mxGetPr(X_OUT);
  id = mxGetPr(ID_OUT);
  
  /* Compute j-th row */
  k=0;
  for (l=1; l<=n; ++l) {
    i = jc[l-1];
    while (i<jc[l] && ir[i]<j) {
      i++;
    }
    if (i<jc[l] && ir[i]==j) {
      x[k] = pr[i];
      id[k] = l;
      k++;
    }
  }
  
}
