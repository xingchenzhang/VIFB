function x=satrousrec(y,fname);

%  SATROUSREC - computes the inverse of 2-D atrous decomposition computed with ATROUSDEC
%  y = satrousdrc(x,fname)
%  INPUT: x image
%         fname    - can be any filter available in the function atrousfilters
%
%  OUTPUT: reconstructed image
%
%  EXAMPLE: xr = satrousrec(y,'9-7');
%
%  History
%   Created on May, 2004 by Arthur Cunha
%   Modified on Aug 2004 by A. C.
%   Modified on Oct 2004 by A. C.
%   SEE ALSO: SATROUSDEC, ATROUSFILTERS

Nlevels=length(y)-1;
[h0,h1,g0,g1] = atrousfilters(fname);


% First Nlevels - 1 levels

x = y{1};
I2 = eye(2);
for i=Nlevels-1:-1:1
     y1=y{Nlevels-i+1};
     shift = -2^(i-1)*[1,1] + 2; % delay correction     
     L=2^i;
     x     = atrousc(symext(x,upsample2df(g0,i),shift),g0,L*I2)+ atrousc(symext(y1,upsample2df(g1,i),shift),g1,L*I2);
 end

% Reconstruct first level

shift=[1,1];
x     = conv2(symext(x,g0,shift),g0,'valid')+ conv2(symext(y{Nlevels+1},g1,shift),g1,'valid');


   
