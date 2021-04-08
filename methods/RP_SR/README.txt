%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This package contains the code which is associated with the following paper:

Yu Liu, Shuping Liu and Zengfu Wang, "A General Framework for Image Fusion Based on Multi-scale Transform and Sparse Representation", Information Fusion, vol. 24, no. 1, pp. 147-164, 2015.

version_1.1: Some "mexw64" files were added into the NSCT and KSVD toolboxes, and the code is effecitve for 64bit-MATLAB now! 

Edited by Yu Liu, 31-01-2015.   

Usage of this code is free for research purposes only. 

Please refer to the above publication if you use this code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Please notice that the package not only provides the implementation of the fusion methods described in the above paper, but also contains all the compared methods employed in the paper. Specifically, six multi-scale transform (MST)-based methods, the sparse representation (SR)-based method, and six MST-SR-based methods are involved in the package. The related MATLAB functions are listed as follows:

lp_fuse.m: Laplacian pyramid (LP)
rp_fuse.m: Ratio of low-pass pyramid (RP)
dwt_fuse.m: Discrete wavelet transform (DWT)
dtcwt_fuse.m: Dual-tree complex wavelet transform (DTCWT)
cvt_fuse.m: Curvelet transform(CVT)
nsct_fuse.m: Nonsubsampled contourlet transform (NSCT)
sparse_fusion.m: Sparse representation (SR) (under the ¡°sparsefusion¡± folder)
lp_sr_fuse.m: LP-SR
rp_sr_fused.m: RP-SR
dwt_sr_fuse.m: DWT-SR
dtcwt_sr_fuse.m: DTCWT-SR
cvt_fuse.m: CVT-SR
nsct_fuse.m: NSCT-SR

To run the code, please use the following three files:
MST_main.m
SR_main.m
MST_SR_main.m 

Therefore, the package is actually a new image fusion toolbox which contains many popular transform domain based fusion methods. It should be noted that many useful functions (or modified versions) contained in the widely used image fusion toolbox (provided by Dr. O. Rockinger, website: http://www.metapix.de/toolbox.htm) are employed in our toolbox, such as:
adb.m, dec.m, dec2.m, es.m, es2.m, selb.m, selc.m, undec.m, undec2.m, lp_fuse.m, rp_fuse.m

In addition, Professor Shutao Li from Hunan university also helps us a lot, the implementations of DTCWT, CVT and NSCT methods in this toolbox are provided by his research group. The related code is used in their following paper:
S. Li, B. Yang and J. Hu, ¡°Performance comparison of different multi-resolution transforms for image fusion¡±, Information Fusion, 2011. 
Moreover, several toolboxes about DTCWT, CVT and NSCT are required, and they have been contained in the package. These toolboxes can be downloaded from MATLAB CENTRAL.

The ksvd toolbox in the ¡°sparsefusion¡± folder is downloaded from Dr. R. Rubinstein¡¯s homepage: http://www.cs.technion.ac.il/~ronrubin/software.html. The function sparse_fusion.m is implemented by us based on the following paper:
B. Yang and S. Li, ¡°Multifocus image fusion and restoration with sparse representation¡±, IEEE Transactions on Instrumentation and Measurement, 2010.

We would like to thank Dr. Rockinger, Professor Li, and all the developers of the related MST and SR toolboxes for their contribution to our MST-SR image fusion toolbox.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Don¡¯t hesitate to contact me if you meet any problems when using this toolbox.
Author: Yu Liu                                                            
Email: liuyu1@mail.ustc.edu.cn; lyuxxz@163.com
Homepage: http://home.ustc.edu.cn/~liuyu1


