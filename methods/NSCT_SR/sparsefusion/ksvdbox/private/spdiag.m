function Y = spdiag(V,K)
%SPDIAG Sparse diagonal matrices.
%   SPDIAG(V,K) when V is a vector with N components is a sparse square
%   matrix of order N+ABS(K) with the elements of V on the K-th diagonal. 
%   K = 0 is the main diagonal, K > 0 is above the main diagonal and K < 0
%   is below the main diagonal. 
%
%   SPDIAG(V) is the same as SPDIAG(V,0) and puts V on the main diagonal.
%
%   See also DIAG, SPDIAGS.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  June 2008


if (nargin<2)
  K = 0;
end

n = length(V) + abs(K);

if (K>0)
  i = 1:length(V);
  j = K+1:n;
elseif (K<0)
  i = -K+1:n;
  j = 1:length(V);
else
  i = 1:n;
  j = 1:n;
end

Y = sparse(i,j,V(:),n,n);