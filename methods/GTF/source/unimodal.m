function[tau, Is] = unimodal(I, N)

if( nargin < 2)
  N = 100;
end

[I_hist, pos] = hist(I(:),N);

[vmax k_max] = max(I_hist);

pmax = pos(k_max);

vend = I_hist(end);
pend = pos(end);

m = (vend - vmax) / (pend - pmax);
alpha = pi/2 - atan(m);

k = pmax:(pend-pmax)/(N-k_max):pend;

dk = sqrt( (vmax-I_hist(k_max:end)).^2 + (pmax - k).^2 );
mk = (I_hist(k_max:end) - vmax) ./ (k-pmax);
alphak = pi/2 - atan(mk);

dpk = dk.*abs(sin(alpha-alphak));

%  figure; plot(dpk);

[d1 p1] = max(dpk);


tau = pmax + pos(p1);
Is = I >= tau;