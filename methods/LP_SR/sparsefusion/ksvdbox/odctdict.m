function D = odctdict(n,L)
%ODCTDICT Overcomplete DCT dictionary.
%  D = ODCTDICT(N,L) returns the overcomplete DCT dictionary of size NxL
%  for signals of length N.
%
%  See also ODCT2DICT, ODCT3DICT, ODCTNDICT.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009


D = zeros(n,L);
D(:,1) = 1/sqrt(n);
for k = 2:L
  v = cos((0:n-1)*pi*(k-1)/L)';
  v = v-mean(v);
  D(:,k) = v/norm(v);
end
