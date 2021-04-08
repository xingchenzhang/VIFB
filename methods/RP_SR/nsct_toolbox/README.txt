Nonsubsampled Contourlet Toolbox (version 1.0.0, August 2004)

This toolbox contains Matlab files that implement the nonsubsample contourlet transform and its utility functions. Some Matlab files are the same as those in Contourlet Toolbox.

The main functions are the following: nsctdec, nsctrec, nsdfbdec, nsdfbrec.
+ nsctdec: Nonsubsampled contourlet decomposition.
+ nsctrec: Nonsubsampled contourlet reconstruction.
+ nsdfbdec: Nonsubsampled directional filter bank decomposition.
+ nsdfbrec: Nonsubsampled directional filter bank reconstruction.
 
In addition, there are several demos (decdemo) that provide examples on how to use these functions.  

Note:   There are three mex file (zconv2.c, zconv2S.c, atrousc.c) in the nonsubsampled Contourlet Toolbox that might need to be recompiled.
            This can be done by typing from the Matlab command window
            >> mex zconv2D_O.c  

            
References: 
1. Jianping Zhou, Arthur L. Cunha, and Minh N. Do, "Nonsubsampled contourlet transform: construction and application in enhancement", submitted, International Conference on Image Processing¡¯05, 2005.
2. Arthur L. Cunha, Jianping Zhou, and Minh N. Do, "Nonsubsampled contourlet transform: filter design and application in image denoising," submitted, International Conference on Image Processing¡¯05, 2005. 
3. Arthur L. Cunha, Jianping Zhou, and Minh N. Do, "Nonsubsampled contourlet transform: Theory, design, and applications", IEEE Transactions on Image Processing, to be submitted.



01/20/05



