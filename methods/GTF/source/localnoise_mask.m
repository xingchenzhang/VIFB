function[sv, mu] = localnoise_mask(I, w, mask)
%  
%  Estimates the local sample variance and local mean value for
%  all the pixels of image I
%  
%  Input parameters
%   I : Input image
%   w :  radius of analysis window
%  
%  Output parameters
%   sv : local sample variance (same size as I)
%   mu : local sample mean value (same size as I)
%  

%
% Legal:
%   localnoise.m is part of the INMD software 
%   (http://sites.google.com/a/istec.net/prodrig/Home/sw). INMD is free 
%   software, you can redistribute it and/or modify it under the terms of 
%   the GNU General Public License (version 2).
%
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe


if(nargin < 2) 
  w = 1; 
end

[Nrows, Ncols] = size(I);

sv = zeros(Nrows, Ncols);
mu = zeros(Nrows, Ncols);


% -----------------

n_min = ((w+1:Ncols-w) - w); 
n_max = ((w+1:Ncols-w) + w);
ind_n = [];

w_size = 2*w + 1;
w_size2 = w_size*w_size;

for n = 1:length(n_min);

  ind_n = [ind_n, (n_min(n):n_max(n))];

end


% --------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   main block
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:Nrows,

  % --------------------
  for n = 1:w,

    [sv(k,n) mu(k,n)] = compute(I, mask, w, k, n, Nrows, Ncols);

  end
  % --------------------

    k_min = k - w; if(k_min <= 0) k_min = 1; end
    k_max = k + w; if(k_max > Nrows) k_max = Nrows; end

  ind_k = (k_min:k_max);

  tmpImg = I(ind_k(:), ind_n(:) );
  tmpMask = mask(ind_k(:), ind_n(:) ) == 0; % valid pixels

  L = w_size*length(ind_k);

  B = reshape(tmpImg, [L, length(n_min)]);

  C = reshape(tmpMask, [L, length(n_min)]);
  factor = sum(C,1); % number of valid pixels

  % mean
  mu_all = mean(C.*B,1);

  mu(k, w+1:Ncols-w) = mu_all*L./(factor + (factor==0));

  % variance
  %    sv(k, w+1:Ncols-w) = var(B);
  sv(k, w+1:Ncols-w) = var(C.*B);

  den = (factor - 1).*( (factor) >= 2 ) + ((factor)<=1);
  sv(k, w+1:Ncols-w) = ( sv(k, w+1:Ncols-w)*(L-1) - (L-factor).*(mu_all.^2) )./den;

  % --------------------
  for n = Ncols-w+1:Ncols,

    [sv(k,n) mu(k,n)] = compute(I, mask, w, k, n, Nrows, Ncols);

  end
  % --------------------

end



% ----------------------------------
% ----------------------------------
% ----------------------------------
function[p_var, p_mu] = compute(I, mask, w, k, n, Nrows, Ncols)

    k_min = k - w; if(k_min <= 0) k_min = 1; end
    k_max = k + w; if(k_max > Nrows) k_max = Nrows; end

    n_min = n - w; if(n_min <= 0) n_min = 1; end
    n_max = n + w; if(n_max > Ncols) n_max = Ncols; end


    p = I(k_min:k_max, n_min:n_max);


    tmp = mask(k_min:k_max, n_min:n_max);
    lmask = find(tmp(:) == 0);

    pvec = p(:);

    p_var = var(pvec(lmask));
    p_mu  = mean(pvec(lmask));
