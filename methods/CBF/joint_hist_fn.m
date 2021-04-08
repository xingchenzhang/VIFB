function [jhist_out]=joint_hist_fn(x1,x2)

%%% joint_hist_fn: Computes the joint histogram of Images x1 and x2.
%%% 
%%% jhist_out=joint_hist_fn(x) 
%%% Example:
%%%   Y = joint_hist_fn(X);  %% takes a pair of images x1 & x2 of equal size and returns the 2d joint histogram.
%%%
%%% Author : B. K. SHREYAMSHA KUMAR 
%%% Created on 21-10-2011.
%%% Updated on 21-10-2011.

[p,q]=size(x1);
jhist_out=zeros(256,256);
for ii=1:p
   for jj=1:q
      jhist_out(x1(ii,jj)+1,x2(ii,jj)+1)=jhist_out(x1(ii,jj)+1,x2(ii,jj)+1)+1;
   end
end