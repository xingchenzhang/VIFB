function Y = selb(M1, M2, mp)
%Y = selb(M1, M2, mp) coefficient selection for base image
%
%    M1  - coefficients A
%    M2  - coefficients B
%    mp  - switch for selection type
%          mp == 1: select A
%          mp == 2: select B
%          mp == 3: average A and B
% 
%    Y   - combined coefficients

%    (Oliver Rockinger 16.08.99)

switch (mp)
  case 1, Y = M1;
  case 2, Y = M2;
  case 3, Y = (M1 + M2) / 2;
  otherwise, error('unknown option');
end;
