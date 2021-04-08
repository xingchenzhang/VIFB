function x = isnr(ref, sig, obs)
%  
% snr -- Compute Improvement Signal-to-Noise Ratio for images
%
% Usage:
%       x = isnr(ref, sig, obs)
%
% Input:
%       ref         Reference image
%       sig         Modified image
%       obs         Observed image
%  
% Output:
%       x           SNR value
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov
%  

mse1 = mean((ref(:)-sig(:)).^2);
mse2 = mean((obs(:)-sig(:)).^2);
x = 10*log10(mse2/mse1);
