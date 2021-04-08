function ycb=cross_bilateral_filt2Df(x1,x2,sigmas,sigmar,ksize)

%%% 2D Cross/Joint Bilateral Filter.
%%% Allow a second image (x2) to shape the kernel’s weights (all dissimilarity comparisons are made 
%%% in second image)and filter weights are applied on first image (x1). 
%%% Both should have same size.

%%% Author : B. K. SHREYAMSHA KUMAR
%%% Created on 13-02-2012.
%%% Updated on 13-02-2012. 

half_ksize=floor(ksize/2);

%%% Gaussian Kernel Generation. 
gk=gauss_ker2D(sigmas,ksize);

%%% To take care of boundaries.
xm1=per_extn_im_fn(x1,ksize); 
xm2=per_extn_im_fn(x2,ksize); 

%%% Cross Bilateral Filter Implementation.
[MM,NN]=size(xm1);
for ii=half_ksize+1:MM-half_ksize
   for jj=half_ksize+1:NN-half_ksize
      xtemp1=xm1(ii-half_ksize:ii+half_ksize,jj-half_ksize:jj+half_ksize);      
      xtemp2=xm2(ii-half_ksize:ii+half_ksize,jj-half_ksize:jj+half_ksize);
      pix_diff=abs(xtemp2-xtemp2(half_ksize+1,half_ksize+1));
      rgk=exp(-(pix_diff.^2/(2*(sigmar)^2)));
      ycb(ii-half_ksize,jj-half_ksize)=sum(sum(xtemp1.*gk.*rgk))/sum(sum(gk.*rgk));
   end
end