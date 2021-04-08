/*******************************************************
* zconv2D_OS.c - Written by Jason Laska
*
* Inputs:   x - A 2D signal
*           h - 2D filter
*           m - upsample matrix of type:
*               Mk^(l) = 2*D0^(l-2)*R3^Sl(k)
*               Sl(k) = 2*floor(k/2) - 2^(l-2) + 1
*               where D0 is a diagonal matrix and R3 
*               is a rotation matrix
*         
* Outputs:  y - 2D result of convolution with filter 
*           upsampled by a parallelogram matrix.
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
    /*Declarations*/
    double *FArray,*SArray,*outArray,*M;
    int SColLength,SRowLength,FColLength,FRowLength,NewFColLength,NewFRowLength;
    int n1,n2,l1,l2,SkipX,SkipY,SkipP;
    double sum;
    int indexx,indexy,outindexx,outindexy,mn1,mn2,mn2save;
    int M0,M1,M2,M3;
    int Start1,End1,Start2,End2;

    /*Get the input sizes*/
    SColLength = mxGetM(SIGNAL); 
    SRowLength = mxGetN(SIGNAL);
    FColLength = mxGetM(FILTER); 
    FRowLength = mxGetN(FILTER);
    
    /*Get The Data*/
    FArray = mxGetPr(FILTER);
    SArray = mxGetPr(SIGNAL);
    M = mxGetPr(MMATRIX);
    M0 = (int)M[0];
    M1 = (int)M[1];
    M2 = (int)M[2];
    M3 = (int)M[3];
    
    /*Calculate New Magical Lengths*/
    NewFRowLength = ((M[0]-1)*(FRowLength-1))+(M[2])*(FColLength-1) + FRowLength - 1;
    NewFColLength = ((M[3]-1)*(FColLength-1))+((M[1])*(FRowLength-1)) + FColLength - 1;

    /*Make output size and Allocate out vector*/
    OUT = mxCreateDoubleMatrix(SColLength, SRowLength, mxREAL); 
    outArray = mxGetPr(OUT);	/*outArray is new vector*/
 
    /*Initialize*/
    sum = 0;
    Start1 = NewFRowLength/2; End1 = SRowLength+NewFRowLength/2;
    Start2 = NewFColLength/2; End2 = SColLength+NewFColLength/2;
    mn1 = Start1%SRowLength;
    mn2 = mn2save = Start2%SColLength;
    
    
    /*Compute!*/
    for(n1=0; n1<SRowLength; n1++){
        for(n2=0; n2<SColLength; n2++)
        {
            outindexx = mn1;
            outindexy = mn2;
            for(l1=0; l1<FRowLength; l1++)
            {   
                indexx = outindexx;
                indexy = outindexy;
                for(l2=0; l2<FColLength; l2++)
                {
                    sum += SArray[LINPOS(indexx,indexy,SColLength)]*FArray[LINPOS(l1,l2,FColLength)];
                    /*indexx -= M2;
                    //if(indexx < 0)
                    //    indexx += SRowLength;
                    //if(indexx > (SRowLength-1))
                    //    indexx -= SRowLength;*/
                    indexy -= M3;    
                    if(indexy < 0)
                        indexy += SColLength;
                }
                outindexx -= M0;
                if(outindexx < 0)
                    outindexx += SRowLength;
                /*outindexy -= M1;    
                //if(outindexy < 0)
                //    outindexy += SColLength;
                //if(outindexy > (SColLength-1))
                //    outindexy -= SColLength;*/
            }
            outArray[LINPOS(n1,n2,SColLength)] = sum;
            sum = 0;
            mn2++;    
            if(mn2 > (SColLength-1))
                mn2 -= SColLength;
        }
        mn2 = mn2save;
        mn1++;
        if(mn1 > (SRowLength-1))
            mn1 -= SRowLength;
    }
        

    return;
}
