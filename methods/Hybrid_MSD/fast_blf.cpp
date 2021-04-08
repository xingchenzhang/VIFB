#include "mex.h" 
#include "include/array.h"
#include "include/fast_blf.h"


using namespace std; 


void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *in[]) 
{
	  double* matIm;
	  int width;
	  int height;
	  int dims;
	  int imSize;
	  int idx;
	  double* outSegInd;
	  int x;
	  int y;

    // Checking number of arguments
	  if(nlhs > 1){
        mexErrMsgTxt("Function has one return values");
        return;
    }

    if(nrhs != 3){
        mexErrMsgTxt("Usage: test(double im,double sigma_s,double sigma_r)");
        return;
    }


	 // Load in arrays and parameters
   matIm = (double*) mxGetData(in[0]);
	 width = mxGetDimensions(in[0])[1] ;
	 height = mxGetDimensions(in[0])[0] ;
	 dims=(int) mxGetNumberOfDimensions(in[0]);
	 
	//mexPrintf("dims=%d\n",dims );

	 if (!((dims == 3)||(dims == 2))){
		mexErrMsgTxt("The input image should be 3 dimension or 2 dimension");
	 }
	 
	imSize = width*height;
    
   // Convert to image. 
	image_type image(width,height);
	
	//for color  image
	if(dims==3){
	   for(y=0;y<height;y++){
		   for(x=0;x<width;x++){
			    idx = x * height + y;
			    image(x,y) = (20.0 * matIm[idx] + 40.0 * matIm[idx + imSize] + 1.0 * matIm[idx + 2 * imSize]) / 61.0 ;
		   }
	   }
  }
  
  //for grayscale  image
  if(dims==2){
	   for(y=0;y<height;y++){
		   for(x=0;x<width;x++){
			    idx = x * height + y;
			    image(x,y) =  matIm[idx] ;
		   }
	   }
  }
  
   // Fast binary filter.  
	double sigma_s = mxGetScalar(in[1]);
	double sigma_r = mxGetScalar(in[2]);
	
	image_type filtered_image(width,height);
	
	Image_filter::fast_LBF(image,image,sigma_s,sigma_r,false,	&filtered_image,&filtered_image);
		
		// Get the result image. 
	out[0] = mxCreateDoubleMatrix((mwSize)height, (mwSize)width, mxREAL);
    outSegInd = mxGetPr(out[0]);

	for ( x = 0; x < width; x++){
        for ( y = 0; y < height; y++){
            idx = x * height + y;
           outSegInd[idx] = filtered_image(x,y);
        }
  }

}
