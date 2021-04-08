/**************************************************************************
 *
 * File name: myblas.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 19.5.2009
 *
 *************************************************************************/


#include "myblas.h"
#include <ctype.h>


/* find maximum of absolute values */

mwIndex maxabs(double c[], mwSize m)
{
  mwIndex maxid=0, k;
  double absval, maxval = SQR(*c);   /* use square which is quicker than absolute value */

  for (k=1; k<m; ++k) {
    absval = SQR(c[k]);
    if (absval > maxval) {
      maxval = absval;
      maxid = k;
    }
  }
  return maxid;
}


/* compute y := alpha*x + y */

void vec_sum(double alpha, double x[], double y[], mwSize n)
{
  mwIndex i;

  for (i=0; i<n; ++i) {
    y[i] += alpha*x[i];
  }
}


/* compute y := alpha*A*x */

void mat_vec(double alpha, double A[], double x[], double y[], mwSize n, mwSize m)
{
  mwIndex i, j, i_n;
  double *Ax;

  Ax = mxCalloc(n,sizeof(double));

  for (i=0; i<m; ++i) {
    i_n = i*n;
    for (j=0; j<n; ++j) {
      Ax[j] += A[i_n+j] * x[i];
    }
  }

  for (j=0; j<n; ++j) {
    y[j] = alpha*Ax[j];
  }

  mxFree(Ax);

}


/* compute y := alpha*A'*x */

void matT_vec(double alpha, double A[], double x[], double y[], mwSize n, mwSize m)
{
  mwIndex i, j, n_i;
  double sum0, sum1, sum2, sum3;

  for (j=0; j<m; ++j) {
    y[j] = 0;
  }

  /* use loop unrolling to accelerate computation */

  for (i=0; i<m; ++i) {
    n_i = n*i;
    sum0 = sum1 = sum2 = sum3 = 0;
    for (j=0; j+4<n; j+=4) {
      sum0 += A[n_i+j]*x[j];
      sum1 += A[n_i+j+1]*x[j+1];
      sum2 += A[n_i+j+2]*x[j+2];
      sum3 += A[n_i+j+3]*x[j+3];
    }
    y[i] += alpha * ((sum0 + sum1) + (sum2 + sum3));
    while (j<n) {
      y[i] += alpha*A[n_i+j]*x[j];
      j++;
    }
  }
}


/* dot product */

double dotprod(double a[], double b[], mwSize n)
{
  double sum = 0;
  mwIndex i;
  for (i=0; i<n; ++i)
    sum += a[i]*b[i];
  return sum;
}


/* find maximum of vector */

mwIndex maxpos(double c[], mwSize m)
{
  mwIndex maxid=0, k;
  double val, maxval = *c;

  for (k=1; k<m; ++k) {
    val = c[k];
    if (val > maxval) {
      maxval = val;
      maxid = k;
    }
  }
  return maxid;
}


/* solve L*x = b */

void backsubst_L(double L[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= L[j*n+i]*x[j];
    }
    x[i] = rhs/L[i*n+i];
  }
}


/* solve L'*x = b */

void backsubst_Lt(double L[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= L[(i-1)*n+j]*x[j];
    }
    x[i-1] = rhs/L[(i-1)*n+i-1];
  }
}


/* solve U*x = b */

void backsubst_U(double U[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= U[j*n+i-1]*x[j];
    }
    x[i-1] = rhs/U[(i-1)*n+i-1];
  }
}


/* solve U'*x = b */

void backsubst_Ut(double U[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= U[i*n+j]*x[j];
    }
    x[i] = rhs/U[i*n+i];
  }
}


/* back substitution solver */

void backsubst(char ul, double A[], double b[], double x[], mwSize n, mwSize k)
{
  if (tolower(ul) == 'u') {
    backsubst_U(A, b, x, n, k);
  }
  else if (tolower(ul) == 'l') {
    backsubst_L(A, b, x, n, k);
  }
  else {
    mexErrMsgTxt("Invalid triangular matrix type: must be ''U'' or ''L''");
  }
}


/* solve equation set using cholesky decomposition */

void cholsolve(char ul, double A[], double b[], double x[], mwSize n, mwSize k)
{
  double *tmp;

  tmp = mxMalloc(k*sizeof(double));

  if (tolower(ul) == 'l') {
    backsubst_L(A, b, tmp, n, k);
    backsubst_Lt(A, tmp, x, n, k);
  }
  else if (tolower(ul) == 'u') {
    backsubst_Ut(A, b, tmp, n, k);
    backsubst_U(A, tmp, x, n, k);
  }
  else {
    mexErrMsgTxt("Invalid triangular matrix type: must be either ''U'' or ''L''");
  }

  mxFree(tmp);
}


/* perform a permutation assignment y := x(ind(1:k)) */

void vec_assign(double y[], double x[], mwIndex ind[], mwSize k)
{
  mwIndex i;

  for (i=0; i<k; ++i)
    y[i] = x[ind[i]];
}


/* print contents of matrix */

void printmat(double A[], int n, int m, char* matname)
{
  int i, j;
  mexPrintf("\n%s = \n\n", matname);

  if (n*m==0) {
    mexPrintf("   Empty matrix: %d-by-%d\n\n", n, m);
    return;
  }

  for (i=0; i<n; ++i) {
    for (j=0; j<m; ++j)
      mexPrintf("   %lf", A[j*n+i]);
    mexPrintf("\n");
  }
  mexPrintf("\n");
}


/* print contents of sparse matrix */

void printspmat(mxArray *a, char* matname)
{
  mwIndex *aJc = mxGetJc(a);
  mwIndex *aIr = mxGetIr(a);
  double *aPr = mxGetPr(a);

  int i;

  mexPrintf("\n%s = \n\n", matname);

  for (i=0; i<aJc[1]; ++i)
    printf("   (%d,1) = %lf\n", aIr[i]+1,aPr[i]);

  mexPrintf("\n");
}

