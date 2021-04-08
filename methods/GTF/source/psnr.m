function x = psnr(ref, sig)
%  
% psnr -- Compute Peak Signal-to-Noise Ratio for images
%
% Usage:
%       x = snr(ref, sig)
%
% Input:
%       ref         Reference image
%       sig         Modified image
%  
% Output:
%       x           SNR value
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov
%  

[Nrows Ncols depth] = size(ref);

mse = sum((ref(:)-sig(:)).^2);

%  dv = Nrows*Ncols*depth*(max(ref(:))^2);
dv = Nrows*Ncols*depth;

x = 10*log10(dv/mse);
