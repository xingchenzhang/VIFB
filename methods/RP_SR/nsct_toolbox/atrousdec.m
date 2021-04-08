function y = atrousdec(x,fname,Nlevels);

%  ATROUSDEC - computes the 2-D atrous decomposition using symmetric extension.
%  y = atrousdec(x,fname,L)
%  INPUT: x image
%         fname    - can be any filter available in the function atrousfilters
%         N levels - number of decomposition levels
%  OUTPUT: y vector cell. the first entry is the lowpass images, the next entries 
%            are the highpass images from coarser to finer scales
%  EXAMPLE: y = atrousdecd(x,'9-7',4)
%
%  History
%   Created on May, 2004 by Arthur Cunha
%   Modified on Aug 2004 by A. C.
%   Modified on Oct 2004 by A. C.
%   SEE ALSO: ATROUSREC, ATROUSFILTERS


[h0,h1,g0,g1] = atrousfilters(fname); % Obtain pyramid filters (the filters must be zero-phase!)
y = cell(1,Nlevels+1);

% First Level

shift = [1, 1]; % delay compensation
y0 = conv2(symext(x,h0,shift),h0,'valid');
y1 = conv2(symext(x,h1,shift),h1,'valid');


% Remaining levels

y{Nlevels+1} = y1;
x  = y0;
I2 = eye(2);
for i=1:Nlevels-1
  shift = -2^(i-1)*[1,1] + 2; L=2^i;
  y0 = atrousc(symext(x,upsample2df(h0,i),shift),h0,I2 * L);
  y1 = atrousc(symext(x,upsample2df(h1,i),shift),h1,I2 * L);
  y{Nlevels-i+1} = y1;
  x=y0;
end
y{1}=x;

