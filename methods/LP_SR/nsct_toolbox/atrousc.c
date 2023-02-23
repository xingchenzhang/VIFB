/******************************************************************
* atrousc.c -  Written by Arthur Cunha. This routine builds up on 
*               zconv2D_OS.c written by Jason Laska
*
* Inputs:   x - A 2D signal
*           h - 2D filter
*           m - separable upsampling matrix
*         
* Outputs:  y - 2D result of convolution with filter 
*           upsampled by a m, only the 'valid' part is returned.
*           Similar to conv2(x,h,'valid'), where h is the upsampled
*           filter.
*  
*          
*
* Usage:    y = zconv2D_O(x,h,m);
*
* Notes:    This function does not actually upsample the filter, 
*           it computes the convolution as if the filter had been 
*           upsampled. This is the ultimate optimized version.
*           Further optimized for separable (diagonal) upsampling matrices.
*
* This is a MEX-FILE for matlab
*
/********************************************************/

#include "mex.h"
#include <math.h>

//Constants for matlab interfacing
#define OUT     plhs[0]
#define SIGNAL  prhs[0] //flip and shift
#define FILTER  prhs[1] //stationary
#define MMATRIX prhs[2]

//MACRO for converting positions to linear
#define LINPOS(row,col,collen) (row*collen)+col


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //Declarations
   double *FArray,*SArray,*outArray,*M;
/* FArray   - Filter coefficients
   SArray   - Signal coefficients
   outArray - Output coefficients
   M        - upsampling matrix 	*/
   int SColLength,SRowLength,FColLength,FRowLength,O_SColLength,O_SRowLength;
   int SFColLength,SFRowLength;
   int n1,n2,l1,l2,k1,k2,f1,f2, kk2, kk1;
   double sum;   
   int M0,M3,sM0,sM3;

    //Get the input sizes
    SColLength = mxGetM(SIGNAL); 
    SRowLength = mxGetN(SIGNAL);
    FColLength = mxGetM(FILTER); 
    FRowLength = mxGetN(FILTER);
    
    SFColLength = FColLength-1;
    SFRowLength = FRowLength-1;
    

	//Get The Data
    FArray = mxGetPr(FILTER);
    SArray = mxGetPr(SIGNAL);
    M = mxGetPr(MMATRIX);
    M0 = (int)M[0];    
    M3 = (int)M[3];   
    sM0 = M0-1;
    sM3 = M3-1;
    

	// Corrected Lengths

	O_SColLength = SColLength - M0*FColLength + 1;
	O_SRowLength = SRowLength - M3*FRowLength + 1;
	
	
    //Make output size and Allocate out vector
    
 
    OUT      = mxCreateDoubleMatrix(O_SColLength, O_SRowLength, mxREAL); 
    outArray = mxGetPr(OUT);	//outArray is new vector
 
	/* Convoluyion loop */

    for (n1=0;n1<O_SRowLength;n1++){
		for (n2=0;n2<O_SColLength;n2++){
			sum=0;		    
		    kk1 = n1 + sM0;;
			for (k1=0;k1<FRowLength;k1++){
  			    kk2 = n2 + sM3;
				for (k2=0;k2<FColLength;k2++){
					 f1 = SFRowLength - k1; /* flipped index */
					 f2 = SFColLength - k2;  		
					 sum+= FArray[LINPOS(f1,f2,FColLength)] * SArray[LINPOS(kk1,kk2,SColLength)];					
					 kk2+=M3;
				}
				kk1+=M0;
			} 
		    outArray[LINPOS(n1,n2,O_SColLength)] = sum;
		}
	}

    return;
}
