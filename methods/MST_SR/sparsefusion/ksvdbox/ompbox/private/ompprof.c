/**************************************************************************
 *
 * File name: ompprof.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 11.4.2009
 *
 *************************************************************************/


#include "ompprof.h"


/* initialize profiling information */

void initprofdata(profdata *pd)
{
  pd->DtX_time = 0;
  pd->XtX_time = 0;
  pd->DtR_time = 0;
  pd->maxabs_time = 0;
  pd->DtD_time = 0;
  pd->Lchol_time = 0;
  pd->compcoef_time = 0;
  pd->update_DtR_time = 0;
  pd->update_resnorm_time = 0;
  pd->compres_time = 0;
  pd->indexsort_time = 0;
  
  pd->DtX_time_counted = 0;
  pd->XtX_time_counted = 0;
  pd->DtR_time_counted = 0;
  pd->DtD_time_counted = 0;
  pd->update_DtR_time_counted = 0;
  pd->resnorm_time_counted = 0;
  pd->compres_time_counted = 0;
  pd->indexsort_time_counted = 0;
  
  pd->prevtime = clock();
}


/* add elapsed time to profiling data according to specified computation */

void addproftime(profdata *pd, int comptype)
{
  switch(comptype) {
    case DtX_TIME:            pd->DtX_time            += clock()-pd->prevtime; pd->DtX_time_counted = 1; break;
    case XtX_TIME:            pd->XtX_time            += clock()-pd->prevtime; pd->XtX_time_counted = 1; break;
    case DtR_TIME:            pd->DtR_time            += clock()-pd->prevtime; pd->DtR_time_counted = 1; break;
    case DtD_TIME:            pd->DtD_time            += clock()-pd->prevtime; pd->DtD_time_counted = 1; break;
    case COMPRES_TIME:        pd->compres_time        += clock()-pd->prevtime; pd->compres_time_counted = 1; break;
    case UPDATE_DtR_TIME:     pd->update_DtR_time     += clock()-pd->prevtime; pd->update_DtR_time_counted = 1; break;
    case UPDATE_RESNORM_TIME: pd->update_resnorm_time += clock()-pd->prevtime; pd->resnorm_time_counted = 1; break;
    case INDEXSORT_TIME:      pd->indexsort_time      += clock()-pd->prevtime; pd->indexsort_time_counted = 1; break;
    case MAXABS_TIME:         pd->maxabs_time         += clock()-pd->prevtime; break;
    case LCHOL_TIME:          pd->Lchol_time          += clock()-pd->prevtime; break;
    case COMPCOEF_TIME:       pd->compcoef_time       += clock()-pd->prevtime; break;
  }
  pd->prevtime = clock();
}


/* print profiling info */

void printprofinfo(profdata *pd, int erroromp, int batchomp, int signum)
{
  clock_t tottime;
  
  tottime = pd->DtX_time + pd->XtX_time + pd->DtR_time + pd->DtD_time + pd->compres_time + pd->maxabs_time + 
            pd->Lchol_time + pd->compcoef_time + pd->update_DtR_time + pd->update_resnorm_time + pd->indexsort_time;
  
  mexPrintf("\n\n*****  Profiling information for %s  *****\n\n", erroromp? "OMP2" : "OMP");
  
  mexPrintf("OMP mode: %s\n\n", batchomp? "Batch-OMP" : "OMP-Cholesky");
  
  mexPrintf("Total signals processed: %d\n\n", signum);
  
  if (pd->DtX_time_counted) {
    mexPrintf("Compute DtX time:      %7.3lf seconds\n", pd->DtX_time/(double)CLOCKS_PER_SEC);
  }
  if (pd->XtX_time_counted) {
    mexPrintf("Compute XtX time:      %7.3lf seconds\n", pd->XtX_time/(double)CLOCKS_PER_SEC);
  }
  mexPrintf("Max abs time:          %7.3lf seconds\n", pd->maxabs_time/(double)CLOCKS_PER_SEC);
  if (pd->DtD_time_counted) {
    mexPrintf("Compute DtD time:      %7.3lf seconds\n", pd->DtD_time/(double)CLOCKS_PER_SEC);
  }
  mexPrintf("Lchol update time:     %7.3lf seconds\n", pd->Lchol_time/(double)CLOCKS_PER_SEC);
  mexPrintf("Compute coef time:     %7.3lf seconds\n", pd->compcoef_time/(double)CLOCKS_PER_SEC);
  if (pd->compres_time_counted) {
    mexPrintf("Compute R time:        %7.3lf seconds\n", pd->compres_time/(double)CLOCKS_PER_SEC);
  }
  if (pd->DtR_time_counted) {
    mexPrintf("Compute DtR time:      %7.3lf seconds\n", pd->DtR_time/(double)CLOCKS_PER_SEC);
  }
  if (pd->update_DtR_time_counted) {
    mexPrintf("Update DtR time:       %7.3lf seconds\n", pd->update_DtR_time/(double)CLOCKS_PER_SEC);
  }
  if (pd->resnorm_time_counted) {
    mexPrintf("Update resnorm time:   %7.3lf seconds\n", pd->update_resnorm_time/(double)CLOCKS_PER_SEC);
  }
  if (pd->indexsort_time_counted) {
    mexPrintf("Index sort time:       %7.3lf seconds\n", pd->indexsort_time/(double)CLOCKS_PER_SEC);
  }
  mexPrintf("---------------------------------------\n");
  mexPrintf("Total time:            %7.3lf seconds\n\n", tottime/(double)CLOCKS_PER_SEC);
}

