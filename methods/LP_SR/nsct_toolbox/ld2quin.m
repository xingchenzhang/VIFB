function [h0, h1] = ld2quin(beta)
% LD2QUIN    Quincunx filters from the ladder network structure
%
% Construct the quincunx filters from an allpass filter (beta) using the
% ladder network structure
%
% Ref: Phong et al., IEEE Trans. on SP, March 1995

if all(size(beta) ~= 1)
    error('The input must be an 1-D fitler');
end

% Make sure beta is a row vector
beta = beta(:)';
    
lf = length(beta);
n = lf / 2;

if n ~= floor(n)
    error('The input allpass filter must be even length');
end

% beta(z1) * beta(z2)
sp = beta' * beta;

% beta(z1*z2^{-1}) * beta(z1*z2)
% Obtained by quincunx upsampling type 1 (with zero padded)
h = qupz(sp, 1);

% Lowpass quincunx filter
h0 = h;
h0(2*n, 2*n) = h0(2*n, 2*n) + 1;
h0 = h0 / 2;

% Highpass quincunx filter
h1 = -conv2(h, h0);
h1(4*n-1, 4*n-1) = h1(4*n-1, 4*n-1) + 1;
