/**************************************************************************
 *
 * File name: omputils.h
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 11.4.2009
 *
 *
 * Various utility definitions and functions for the OMP library.
 *
 *************************************************************************/


#ifndef __OMP_UTILS_H__
#define __OMP_UTILS_H__

#include "mex.h"


/* constants for the representation mode of gamma */

extern const char FULL_GAMMA_STR[];      /* "full" */
extern const char SPARSE_GAMMA_STR[];    /* "sparse" */


#define FULL_GAMMA 1
#define SPARSE_GAMMA 2
#define INVALID_MODE 3



/**************************************************************************
 * Function checkmatrix:
 *
 * Verify that the specified mxArray is real, of type double, and has 
 * no more than two dimensions. If not, an error message is printed
 * and the mex file terminates.
 * 
 * Parameters:
 *   param - the mxArray to be checked
 *   name  - the name of the parameter (15 characters or less)
 *
 **************************************************************************/
void checkmatrix(const mxArray *param, char *name);


/**************************************************************************
 * Function checkscalar:
 *
 * Verify that the specified mxArray represents a real double scalar value. 
 * If not, an error message is printed and the mex file terminates.
 * 
 * Parameters:
 *   param - the mxArray to be checked
 *   name  - the name of the parameter (15 characters or less)
 *
 **************************************************************************/
void checkscalar(const mxArray *param, char *name);



/**************************************************************************
 * Convert number of seconds to hour, minute and second representation.
 *
 * Parameters:
 *   sectot - total number of seconds
 *   hrs, mins, secs - output hours (whole) and minutes (whole) and seconds
 *
 **************************************************************************/
void secs2hms(double sectot, int *hrs, int *mins, double *secs);


#endif

