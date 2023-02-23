function yT=symext(x,h,shift);

% FUNCTION Y = SYMEXT
% INPUT:  x, mxn image
%         h, 2-D filter coefficients
%         shift, optional shift
% OUTPUT: yT image symetrically extended (H/V symmetry)
%
% Performs symmetric extension for image x, filter h. 
% The filter h is assumed have odd dimensions.
% If the filter has horizontal and vertical symmetry, then 
% the nonsymmetric part of conv2(h,x) has the same size of x.
%
% Created by A. Cunha, Fall 2003;
% Modified 12/2005 by A. Cunha. Fixed a bug on wrongly 
% swapped indices (m and n). 

[m,n] = size(x);
[p,q] = size(h);
parp  = 1-mod(p,2) ;
parq  = 1-mod(q,2);

p2=floor(p/2);q2=floor(q/2);
s1=shift(1);s2=shift(2);

ss = p2 - s1 + 1;
rr = q2 - s2 + 1;

yT = [fliplr(x(:,1:ss)) x  x(:,n  :-1: n-p-s1+1)];
yT = [flipud(yT(1:rr,:)); yT ;  yT(m  :-1: m-q-s2+1,:)];
yT = yT(1:m+p-1 ,1:n+q-1);
 
