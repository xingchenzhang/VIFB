function D = odct3dict(n,L)
%ODCT3DICT 3-D overcomplete DCT dictionary.
%  D = ODCT3DICT([N1 N2 N3],[L1 L2 L3]) returns an overcomplete DCT
%  dictionary for signals of size N1xN2xN3. The number of dictionary atoms
%  in the i-th dimension is Li, so the combined dictionary is of size
%  (N1*N2*N3) x (L1*L2*L3).
%
%  D = ODCT3DICT([N1 N2 N3],L) specifies the total number of atoms in the
%  dictionary instead of each of the Li's individually. The Li's in this
%  case are selected so their relative sizes are roughly the same as the
%  relative sizes of the Ni's. Note that the actual number of atoms in the
%  dictionary may be larger than L, as rounding might be required for the
%  computation of the Li's.
%
%  D = ODCT3DICT(N,L) is shorthand for the call ODCT2DICT([N N N],L), and
%  returns the overcomplete DCT dictionary for signals of size NxNxN. L is
%  the required size of the overcomplete dictionary, and is rounded up to
%  the nearest cubic integer.
%
%  See also ODCTDICT, ODCT2DICT, ODCTNDICT.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009

D = odctndict(n,L,3);
