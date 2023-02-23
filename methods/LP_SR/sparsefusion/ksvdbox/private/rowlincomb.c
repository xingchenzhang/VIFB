/**************************************************************************
 *
 * File name: rowlincomb.c
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

#define	X_IN	   prhs[0]
#define A_IN     prhs[1]
#define ROWS_IN  prhs[2]
#define COLS_IN  prhs[3]


/* Output Arguments */

#define	Y_OUT	plhs[0]


void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray*prhs[])

{
  double *A, *x, *y, *rows, *cols;
  mwSize m,n,m1,n1,m2,n2,rownum,colnum;
  mwIndex *row_ids,*col_ids,i,j;
  int colnumspecified=0;
  
  
  /* Check for proper number of arguments */
  
  if (nrhs!=3 && nrhs!=4) {
    mexErrMsgTxt("Invalid number of input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  
  /* Check the input dimensions */
  
  m = mxGetM(A_IN);
  n = mxGetN(A_IN);
  if (!mxIsDouble(A_IN) || mxIsComplex(A_IN) || mxGetNumberOfDimensions(A_IN)>2) {
    mexErrMsgTxt("ROWLINCOMB requires that A be a double matrix.");
  }
  
  m1 = mxGetM(ROWS_IN);
  n1 = mxGetN(ROWS_IN);
  if (!mxIsDouble(ROWS_IN) || mxIsComplex(ROWS_IN) || (m1!=1 && n1!=1)) {
    mexErrMsgTxt("ROWLINCOMB requires that ROWS be an index vector of type double.");
  }
  rownum = (m1 > n1) ? m1 : n1;   /* the number of rows in the linear combination */
  
  m2 = mxGetM(X_IN);
  n2 = mxGetN(X_IN);
  if (!mxIsDouble(X_IN) || mxIsComplex(X_IN) || ((m2!=1) && (n2!=1))) {
    mexErrMsgTxt("ROWLINCOMB requires that X be a double vector.");
  }
  
  if (m2 != rownum && n2 != rownum) {
    mexErrMsgTxt("The length of X does not match the number of rows in ROWS.");
  }
  
  if (nrhs==4) {
    m1 = mxGetM(COLS_IN);
    n1 = mxGetN(COLS_IN);
    if (!mxIsDouble(COLS_IN) || mxIsComplex(COLS_IN) || (m1!=1 && n1!=1)) {
      mexErrMsgTxt("ROWLINCOMB requires that COLS be an index vector of type double.");
    }
    colnum = (m1 > n1) ? m1 : n1;   /* the number of columns */
    colnumspecified = 1;
    cols = mxGetPr(COLS_IN);
    
    Y_OUT = mxCreateDoubleMatrix(1, colnum, mxREAL);
  }
  else {
    cols = 0;
    Y_OUT = mxCreateDoubleMatrix(1, n, mxREAL);
  }
  
  
  /* Assign pointers to the various parameters */
  A = mxGetPr(A_IN);
  rows = mxGetPr(ROWS_IN);
  x = mxGetPr(X_IN);
  y = mxGetPr(Y_OUT);
  
  
  /* check row indices */
  
  row_ids = (mwIndex*)mxMalloc(rownum*sizeof(mwIndex));
  
  for (i=0; i<rownum; ++i) {
    row_ids[i] = (mwIndex)(rows[i]+0.1)-1;
    if (row_ids[i]<0 || row_ids[i]>=m) {
      mexErrMsgTxt("Row index in ROWS is out of range.");
    }
  }
  
  
  
  if (colnumspecified) {
    
    /* check column indices */
    col_ids = (mwIndex*)mxMalloc(colnum*sizeof(mwIndex));
    
    for (i=0; i<colnum; ++i) {
      col_ids[i] = (mwIndex)(cols[i]+0.1)-1;
      if (col_ids[i]<0 || col_ids[i]>=n) {
        mexErrMsgTxt("Column index in COLS is out of range.");
      }
    }
    
    /* Do the actual computation */
    for (j=0; j<colnum; ++j) {
      for (i=0; i<rownum; ++i) {
        y[j] += A[m*col_ids[j]+row_ids[i]]*x[i];
      }
    }
    
    mxFree(col_ids);
  }
  
  else {
    
    /* Do the actual computation */
    for (j=0; j<n; ++j) {
      for (i=0; i<rownum; ++i) {
        y[j] += A[m*j+row_ids[i]]*x[i];
      }
    }
  }
  
  
  mxFree(row_ids);
  
  return;
}
