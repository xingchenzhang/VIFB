function x = reErr(ref, sig)
%  
% psnr -- Compute relative error for images
%
% Usage:
%       x = reErr(ref, sig)
%
% Input:
%       ref         Reference image
%       sig         Modified image
%  
% Output:
%       x           reErr value
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov
%  

[Nrows Ncols depth] = size(ref);

num = mean((ref(:)-sig(:)).^2);
den = mean((sig(:)).^2);

x = sqrt(num/den);