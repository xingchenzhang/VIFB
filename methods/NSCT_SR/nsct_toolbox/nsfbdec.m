function [y0, y1] = nsfbdec( x, h0, h1, lev )

%  nsfbdec - computes the ns pyramid decomposition 
%  at level lev with filters f1, f2
%  y = nsfbdec(x,h0,h1,L)
%  INPUT: x image at finer scale
%           h0, h1 atrous filters obtained from 'atrousfilters' 
%  OUTPUT: y1 - image at coarser scale, y2 - wavelet highpass output
% 
%  History
%   Adapted from atrousdec by Arthur Cunha
%   Modified on Aug 2004 by A. C.
%   Modified on Oct 2004 by A. C.
%   SEE ALSO: ATROUSREC, ATROUSFILTERS


 
if lev ~= 0 
    I2 = eye(2); % delay compensation
    shift = -2^(lev-1)*[1,1] + 2; L=2^lev;
    y0 = atrousc(symext(x,upsample2df(h0,lev),shift),h0,I2 * L);
    y1 = atrousc(symext(x,upsample2df(h1,lev),shift),h1,I2 * L); 
else
    % first Level
    shift = [1, 1]; % delay compensation
    y0 =  conv2(symext(x,h0,shift),h0,'valid');
    y1 =  conv2(symext(x,h1,shift),h1,'valid');
end
  

