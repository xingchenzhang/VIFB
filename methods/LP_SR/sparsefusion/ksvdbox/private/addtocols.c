/**************************************************************************
 *
 * File name: addtocols.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 19.4.2009
 *
 *************************************************************************/


#include "mex.h"


/* Input Arguments */

#define	X_IN	prhs[0]
#define V_IN  prhs[1]


/* Output Arguments */

#define	Y_OUT	plhs[0]


void mexFunction(int nlhs, mxArray *plhs[], 
		             int nrhs, const mxArray*prhs[])
     
{ 
    double *x, *y, *v, *xend;  
    mwSize m,n,m1,n1;
    mwIndex counter; 
    
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) {
      mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
      mexErrMsgTxt("Too many output arguments."); 
    } 
    
    
    /* Check the the input dimensions */ 
    
    m = mxGetM(X_IN);
    n = mxGetN(X_IN);
    if (!mxIsDouble(X_IN) || mxIsComplex(X_IN) || mxGetNumberOfDimensions(X_IN)>2) {
      mexErrMsgTxt("ADDTOCOLS requires that X be a double matrix.");
    }
    m1 = mxGetM(V_IN);
    n1 = mxGetN(V_IN);
    if (!mxIsDouble(V_IN) || mxIsComplex(V_IN) || (m1!=1 && n1!=1)) {
      mexErrMsgTxt("ADDTOCOLS requires that V be a double vector.");
    }
    if ((m1==1 && n1!=n) || (n1==1 && m1!=n)) {
      mexErrMsgTxt("Error in ADDTOCOLS: dimensions of V and X must agree.");
    }
    
    
    /* Create a matrix for the return argument */ 
    Y_OUT = mxCreateDoubleMatrix(m, n, mxREAL);
         
    
    /* Assign pointers to the various parameters */ 
    x = mxGetPr(X_IN); 
    v = mxGetPr(V_IN);
    y = mxGetPr(Y_OUT);
            
    
    /* Do the actual computation */
    
    xend = x+(m*n);
    counter = 0;
    while (x<xend) {
      (*y) = (*x) + (*v);
      y++; x++; counter++;
      if (counter==m) {v++; counter=0;}
    }
    
    return;
}
