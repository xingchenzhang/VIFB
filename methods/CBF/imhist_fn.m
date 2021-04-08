function [hist_out]=imhist_fn(x)

%%% imhist_fn: Computes the histogram of Image x.
%%% 
%%% hist_out=imhist_fn(x) 
%%% Example:
%%%   Y = imhist_fn(X);    %% Computes the histogram of Image x.
%%%
%%% Author : B.K. SHREYAMSHA KUMAR 
%%% Created on 19-10-2011.
%%% Updated on 19-10-2011.

[p,q]=size(x);
hist_out=zeros(1,256);
for ii=1:p
   for jj=1:q
      hist_out(x(ii,jj)+1)=hist_out(x(ii,jj)+1)+1;
   end
end