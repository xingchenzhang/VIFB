/**************************************************************************
 *
 * File name: collincomb.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 21.5.2009
 *
 *************************************************************************/


#include "mex.h"


/* Input Arguments */

#define	A_IN	   prhs[0]
#define ROWS_IN  prhs[1]
#define COLS_IN1 prhs[1]
#define COLS_IN2 prhs[2]
#define X_IN1    prhs[2]
#define X_IN2    prhs[3]


/* Output Arguments */

#define	Y_OUT	plhs[0]


void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray*prhs[])

{
  double *A, *x, *y, *rows, *cols;
  mwSize m,n,m1,n1,m2,n2,rownum,colnum;
  mwIndex col_id,*row_ids,i,j;
  int rownumspecified=0;
  
  
  /* Check for proper number of arguments */
  
  if (nrhs!=3 && nrhs!=4) {
    mexErrMsgTxt("Invalid number of arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  
  /* Check the input dimensions */
  
  m = mxGetM(A_IN);
  n = mxGetN(A_IN);
  if (!mxIsDouble(A_IN) || mxIsComplex(A_IN) || mxGetNumberOfDimensions(A_IN)>2) {
    mexErrMsgTxt("COLLINCOMB requires that A be a double matrix.");
  }
  
  if (nrhs==3) {
    
    m1 = mxGetM(COLS_IN1);
    n1 = mxGetN(COLS_IN1);
    if (!mxIsDouble(COLS_IN1) || mxIsComplex(COLS_IN1) || (m1!=1 && n1!=1)) {
      mexErrMsgTxt("COLLINCOMB requires that COLS be an index vector of type double.");
    }
    colnum = (m1 > n1) ? m1 : n1;   /* the number of columns in the linear combination */
    
    m2 = mxGetM(X_IN1);
    n2 = mxGetN(X_IN1);
    if (!mxIsDouble(X_IN1) || mxIsComplex(X_IN1) || (m2!=1 && n2!=1)) {
      mexErrMsgTxt("COLLINCOMB requires that X be a double vector.");
    }
    
    if (m2!=colnum && n2!=colnum) {
      mexErrMsgTxt("The length of X does not match the number of columns in COLS.");
    }
    
    rows = 0;
    Y_OUT = mxCreateDoubleMatrix(m, 1, mxREAL);
    cols = mxGetPr(COLS_IN1);
    x = mxGetPr(X_IN1);
  }
  else {
    
    m1 = mxGetM(ROWS_IN);
    n1 = mxGetN(ROWS_IN);
    if (!mxIsDouble(ROWS_IN) || mxIsComplex(ROWS_IN) || (m1!=1 && n1!=1)) {
      mexErrMsgTxt("COLLINCOMB requires that ROWS be an index vector of type double.");
    }
    rownum = (m1 > n1) ? m1 : n1;   /* the number of rows in the linear combination */
    rownumspecified = 1;
    rows = mxGetPr(ROWS_IN);
    
    m1 = mxGetM(COLS_IN2);
    n1 = mxGetN(COLS_IN2);
    if (!mxIsDouble(COLS_IN2) || mxIsComplex(COLS_IN2) || (m1!=1 && n1!=1)) {
      mexErrMsgTxt("COLLINCOMB requires that COLS be an index vector of type double.");
    }
    colnum = (m1 > n1) ? m1 : n1;   /* the number of columns in the linear combination */
    
    m2 = mxGetM(X_IN2);
    n2 = mxGetN(X_IN2);
    if (!mxIsDouble(X_IN2) || mxIsComplex(X_IN2) || (m2!=1 && n2!=1)) {
      mexErrMsgTxt("COLLINCOMB requires that X be a double vector.");
    }
    
    if (m2!=colnum && n2!=colnum) {
      mexErrMsgTxt("The length of X does not match the number of columns in COLS.");
    }
    
    Y_OUT = mxCreateDoubleMatrix(rownum, 1, mxREAL);
    cols = mxGetPr(COLS_IN2);
    x = mxGetPr(X_IN2);
  }
  
  
  /* Assign pointers to the various parameters */
  A = mxGetPr(A_IN);
  y = mxGetPr(Y_OUT);
  
  
  if (rownumspecified) {
    
     /* check row indices */
    
    row_ids = (mwIndex*)mxMalloc(rownum*sizeof(mwIndex));
    
    for (i=0; i<rownum; ++i) {
      row_ids[i] = (mwIndex)(rows[i]+0.1)-1;
      if (row_ids[i]<0 || row_ids[i]>=m) {
        mexErrMsgTxt("Row index in ROWS is out of range.");
      }
    }
    
    /* Do the actual computation */
    for (i=0; i<colnum; ++i) {
      col_id = (mwIndex)(cols[i]+0.1)-1;
      if (col_id<0 || col_id>=n) {
        mexErrMsgTxt("Column index in COLS is out of range.");
      }
      for (j=0; j<rownum; ++j) {
        y[j] += A[m*col_id+row_ids[j]]*x[i];
      }
    }
    
    mxFree(row_ids);
  }
  
  else {
    
    /* Do the actual computation */
    for (i=0; i<colnum; ++i) {
      col_id = (mwIndex)(cols[i]+0.1)-1;
      if (col_id<0 || col_id>=n) {
        mexErrMsgTxt("Column index in COLS is out of range.");
      }
      for (j=0; j<m; ++j) {
        y[j] += A[m*col_id+j]*x[i];
      }
    }
  }
  
  return;
}
