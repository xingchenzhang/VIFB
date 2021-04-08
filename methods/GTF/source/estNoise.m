function[sigma2, cv, mu, sig2inv, cvinv, S, lmu ] = estNoise(I, w, N, zoom);
%
%  Input parameters
%   I     : Input (noisy) image
%   w     : radius of analysis window
%   N     : Number of ellements to be considered (histogram) when
%           estimating sigma2 and cv2
%   zoom  : Conditional flag to improve the estimation of
%           sigma2 and cv2
%  
%  
%  Output parameters
%  
%   sigma2  : estimation of the variance
%   cv2     : estimation of the coeficient of variation (cv2 = sigma^2 / mu^2)
%   mu      : estimation of the mean value
%   sig2inv : estimation of 1/sigma2
%   cvinv   : estimation of 1/cv2
%   S       : estimation of the local sample variance (same size as I)
%   lmu     : estimation of the local sample mean (same size as I)
%    

% Based on 
%   S. Aja-Fernandez, G. Vegas-Sanchez-Ferrero,
%   M. Martin-Fernandez, C. Alberola-Lopez
%   "Automatic noise estimation in images using local statistics.
%   Additive and multiplicative cases", 
%   Image and Vision Computing 27 (2009) 756-770
%
% Legal:
%   estNoise.m is part of the INMD software 
%   (http://sites.google.com/a/istec.net/prodrig/Home/sw). INMD is free 
%   software, you can redistribute it and/or modify it under the terms of 
%   the GNU General Public License (version 2).
%
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe


if(nargin < 4)
  zoom = 0;
  if(nargin < 3) 
    N = 100; 
  end
end




  L=10;
  b = 0.5*( 1 - cos(pi*(0:2*L)/L) );

  [Nrows, Ncols, depth] = size(I);


for d = 1:depth,

  % Local noise estimation
  [S(:,:,d) lmu(:,:,d)] = localnoise(I(:,:,d),w);


  sigma2(d) = estimateMode(S(:,:,d), N, zoom, b, L, w, 0, 1);
  sig2inv(d) = estimateMode(1./S, N, zoom, b, L, w, 1, 1);

  mu(d) = estimateMode(lmu(:,:,d), N, zoom, b, L, w, 0, 1);

  cv(d) = estimateMode(S(:,:,d)./(lmu(:,:,d).*lmu(:,:,d)), N, zoom, b, L, w, 0, 0);
  cvinv(d) = estimateMode((lmu(:,:,d).*lmu(:,:,d))./S(:,:,d), N, zoom, b, L, w, 1, 0);


end

%    sigma2 = estimateMode(S, N, zoom, b, L, w, 0, 1);
%    sig2inv = estimateMode(1./S, N, zoom, b, L, w, 1, 1);
%  
%    mu = estimateMode(lmu, N, zoom, b, L, w, 0, 1);
%  
%    cv = estimateMode(S./(lmu.*lmu), N, zoom, b, L, w, 0, 0);
%    cvinv = estimateMode((lmu.*lmu)./S, N, zoom, b, L, w, 1, 0);





%=================================================================

function[sigma2] = estimateMode(S, N, zoom, b, L, w, inv, factor)


v = S(:);

vmin = min(v);
vmax = max(v);
res = (vmax - vmin)/N;

[h bin] = hist(v, N);

if(inv == 1)
  tmp = conv(b, h) / (2*L+1);
  hs = tmp(L:N+L-1);
%    [dummy pos] = max(hs(end:-1:1));
  posvec = findmaxima( hs(end:-1:1) );
  if(length(posvec) >= 1) 
    pos = posvec(1);
  else
    [dummy pos] = max(h(end:-1:1));
  end
else
  [dummy pos] = max(h(end:-1:1));
end

if(zoom)
  
  ind = find( abs( (v - bin(N-pos+1))  ) < 2*res );

  [h bin] = hist(v(ind), N);

  tmp = conv(b, h) / (2*L+1);
  hs = tmp(L:N+L-1);

  maxpos = findmaxima(hs);


  if( length(maxpos) >= 2 )
    sval = bin( maxpos(1) );
  else
    [dummy pos] = max(h(end:-1:1));
    sval = bin(N-pos+1);
  end

else
  sval = bin(N-pos+1);
end

if factor == 1,
  w_size = (2*w + 1)*(2*w + 1);
  sigma2 = (w_size - 3)/(w_size - 1)*sval;
else
  sigma2 = sval;
end

%=================================================================

function minima = findminima(x)

minima = findmaxima(-x);


function maxima = findmaxima(x)


% Unwrap to vector
x = x(:);
% Identify whether signal is rising or falling
upordown = sign(diff(x));
% Find points where signal is rising before, falling after
maxflags = [upordown(1)<0; diff(upordown)<0; upordown(end)>0];
maxima   = find(maxflags);


