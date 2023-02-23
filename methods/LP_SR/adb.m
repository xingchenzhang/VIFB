function Y = adb(X, bd)
%Y = ADB(X, bd) adds rows and/or columns by duplication on the 
%    lower resp. right side of the input matrix
%    
%    X     - input matrix
%    bd(1) - number of rows to add 
%    bd(2) - number of columns to add
%  
%    Y     - extended matrix

%    (Oliver Rockinger 16.08.99)

[z s] = size(X);

% copy interior
Y          = zeros(z+bd(1),s+bd(2));
Y(1:z,1:s) = X;

% add rows
if (bd(1) > 0)
  Y(z+1:z+bd(1),1:s) = X(z-1:-1:z-bd(1),1:s);
end;
% add columns
if (bd(2) > 0)
  Y(1:z,s+1:s+bd(2)) = X(1:z,s-1:-1:s-bd(2));
end;
% add corner
if (bd(1) > 0 & bd(2) > 0)
  Y(z+1:z+bd(1),s+1:s+bd(2)) = X(z-1:-1:z-bd(1),s-1:-1:s-bd(2));
end;  