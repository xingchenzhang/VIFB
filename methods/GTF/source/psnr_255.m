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


Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));

[Nrows Ncols depth] = size(ref);

mse = sum((255*Normalize(ref(:))-255*Normalize(sig(:))).^2);
dv = Nrows*Ncols*depth*(255*255);
x = 10*log10(dv/mse);
