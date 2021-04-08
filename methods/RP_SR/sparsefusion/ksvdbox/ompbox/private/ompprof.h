/**************************************************************************
 *
 * File name: ompprof.h
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 11.4.2009
 *
 *
 * Collection of definitions and functions for profiling the OMP method.
 *
 *************************************************************************/


#ifndef __OMP_PROF_H__
#define __OMP_PROF_H__

#include "mex.h"
#include <time.h>



/**************************************************************************
 *
 * Constants and data types.
 *
 **************************************************************************/


/* constants denoting the various parts of the algorithm */

enum { DtX_TIME, XtX_TIME, DtR_TIME, MAXABS_TIME, DtD_TIME, LCHOL_TIME, COMPCOEF_TIME, 
       UPDATE_DtR_TIME, UPDATE_RESNORM_TIME, COMPRES_TIME, INDEXSORT_TIME };

       
       
/* profiling data container with counters for each part of the algorithm */
       
typedef struct profdata 
{
  clock_t prevtime;  /* the time when last initialization/call to addproftime() was performed */
  
  clock_t DtX_time;
  clock_t XtX_time;
  clock_t DtR_time;
  clock_t maxabs_time;
  clock_t DtD_time;
  clock_t Lchol_time;
  clock_t compcoef_time;
  clock_t update_DtR_time;
  clock_t update_resnorm_time;
  clock_t compres_time;
  clock_t indexsort_time;
  
  /* flags indicating whether profiling data was gathered */
  int DtX_time_counted;
  int XtX_time_counted;
  int DtR_time_counted;
  int DtD_time_counted;
  int update_DtR_time_counted;
  int resnorm_time_counted;
  int compres_time_counted;
  int indexsort_time_counted;
  
} profdata;



/**************************************************************************
 *
 * Initialize a profdata structure, zero all counters, and start its timer.
 *
 **************************************************************************/
void initprofdata(profdata *pd);


/**************************************************************************
 *
 * Add elapsed time from last call to addproftime(), or from initialization
 * of profdata, to the counter specified by comptype. comptype must be one
 * of the constants in the enumeration above.
 *
 **************************************************************************/
void addproftime(profdata *pd, int comptype);


/**************************************************************************
 *
 * Print the current contents of the counters in profdata.
 *
 * Parameters:
 *   pd - the profdata to print
 *   erroromp - indicates whether error-based (nonzero) or sparsity-based (zero)
 *              omp was performed.
 *   batchomp - indicates whether batch-omp (nonzero) or omp-cholesky (zero)
 *              omp was performed.
 *   signum   - number of signals processed by omp
 *
 **************************************************************************/
void printprofinfo(profdata *pd, int erroromp, int batchomp, int signum);


#endif

