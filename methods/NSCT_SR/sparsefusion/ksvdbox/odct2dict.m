function D = odct2dict(n,L)
%ODCT2DICT 2-D overcomplete DCT dictionary.
%  D = ODCT2DICT([N1 N2],[L1 L2]) returns an overcomplete DCT dictionary
%  for signals of size N1xN2. The number of dictionary atoms in the i-th
%  dimension is Li, so the combined dictionary is of size (N1*N2)x(L1*L2).
%
%  D = ODCT2DICT([N1 N2],L) specifies the total number of atoms in the
%  dictionary, instead of L1 and L2 individually. L1 and L2 in this case
%  are selected so that L1/L2 ~ N1/N2, by letting L1=ceil(sqrt(L*N1/N2))
%  and L2=ceil(sqrt(L*N2/N1)). Note that the actual number of atoms in the
%  dictionary may be larger than L.
%
%  D = ODCT2DICT(N,L) is shorthand for the call ODCT2DICT([N N],L), and
%  returns the overcomplete DCT dictionary for signals of size NxN. L is
%  the required size of the overcomplete dictionary, and is rounded up to
%  the nearest square integer.
%
%  See also ODCTDICT, ODCT3DICT, ODCTNDICT.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009

D = odctndict(n,L,2);
