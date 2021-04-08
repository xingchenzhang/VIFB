function x = nsfbrec( y0, y1, g0, g1, lev )
%  nsfbrec - computes the inverse of 2-D atrous decomposition at level lev
%  y =nsfb(x,fname)
%  INPUT: lowpass image y0 and highpass y1        
%  OUTPUT: reconstructed image at finer scale
%
%  EXAMPLE: xr = nsfbrec(y0,y1,g0,g1);
%
%  History
%   Created on May, 2004 by Arthur Cunha
%   Modified on Aug 2004 by A. C.
%   Modified on Oct 2004 by A. C.
%   SEE ALSO: NSFBDEC, ATROUSFILTERS

I2 = eye(2);
if lev ~= 0
     shift = -2^(lev-1)*[1,1] + 2; % delay correction     
     L=2^lev;
     x     = atrousc(symext(y0,upsample2df(g0,lev),shift),g0,L*I2) + ...
             atrousc(symext(y1,upsample2df(g1,lev),shift),g1,L*I2);
 else     
shift=[1,1];
x     = conv2(symext(y0,g0,shift),g0,'valid')+ conv2(symext(y1,g1,shift),g1,'valid');
end

   
