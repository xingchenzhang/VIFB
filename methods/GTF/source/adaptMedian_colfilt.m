function[S] = adaptMedian_colfilt(I, w_max)
%  
%  

if nargin < 2
  w_max = 9;
end

[Nrows Ncols depth] = size(I);
S = uint8(zeros(Nrows, Ncols, depth));




for d = 1:depth,

% Get plane d
Ib = I(:,:,d);

% Initialization
w = 1;
p = 3;
cond = uint8(1);    % zero --> noiseless pixel


for l = w+2:2:w_max


  Imed = medfilt2(Ib, [p p]);
  Imin = colfilt(Ib, [p p], 'sliding', @(x) min(x));
  Imax = colfilt(Ib, [p p], 'sliding', @(x) max(x));


  m1 = uint8( ((Imin < Imed) .* (Imed < Imax)) );
  m2 = uint8( ((Imin < Ib) .* (Ib < Imax)) );

  if(l == w)
    S(:,:,d) =  m1.*(1-m2).*l;
  else
    S(:,:,d) = cond.*( m1.*(1-m2).*l ) + (1-cond).*S(:,:,d);
  end

  cond = cond.*(1-m1); % pixel that need further evaluation

  p = p+2;

end % _END_ FOR(l)

S(:,:,d) = cond.*(w_max+2) + (1-cond).*S(:,:,d);

end % _END_ FOR(d)

S = double(S);
