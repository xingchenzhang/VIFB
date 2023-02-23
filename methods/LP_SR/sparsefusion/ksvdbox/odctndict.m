function D = odctndict(n,L,p)
%ODCTNDICT Multi-dimensional overcomplete DCT dictionary.
%  D = ODCTNDICT([N1 N2 ... Np],[L1 L2 ... Lp]) returns an overcomplete DCT
%  dictionary for p-dimensional signals of size N1xN2x...xNp. The number of
%  DCT atoms in the i-th dimension is Li, so the combined dictionary is of
%  size (N1*N2*...*Np) x (L1*L2*...*Lp).
%
%  D = ODCTNDICT([N1 N2 ... Np],L) specifies the total number of atoms in
%  the dictionary instead of each of the Li's individually. The Li's in
%  this case are selected so their relative sizes are roughly the same as
%  the relative sizes of the Ni's. Note that the actual number of atoms in
%  the dictionary may be larger than L, as rounding might be required for
%  the computation of the Li's.
%
%  D = ODCTNDICT(N,L,P) is shorthand for the call ODCTNDICT(N*ones(1,P),L),
%  and returns the overcomplete DCT dictionary for P-dimensional signals of
%  size NxNx...xN. L is the required size of the overcomplete dictionary,
%  and is rounded up to the nearest integer with a whole P-th root.
%
%  See also ODCTDICT, ODCT2DICT, ODCT3DICT.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009

if (nargin==2)
  p = length(n);
end

if (length(n)==1)
  n = n*ones(1,p);
end

if (length(L)==1 && p>1)
  N = prod(n);
  L = ceil((L*(n.^p/N).^(1/(p-1))).^(1/p));
end

D = odctdict(n(1),L(1));
for i = 2:p
  D = kron(D,odctdict(n(i),L(i)));
end
