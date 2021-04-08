function x = snr(ref, sig)
%  
% snr -- Compute Signal-to-Noise Ratio for images
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

mse = mean((ref(:)-sig(:)).^2);
dv = var(ref(:),1);
x = 10*log10(dv/mse);
