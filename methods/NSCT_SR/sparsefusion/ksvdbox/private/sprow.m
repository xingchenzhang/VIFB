%SPROW Extract row of sparse matrix.
%  X = SPROW(A,J) where A is a sparse matrix, returns the nonzero values in
%  the row A(J,:).
%
%  [X,ID] = SPROW(A,J) also returns the column indices of the nonzeros.
%
%  Note that the call [X,ID] = SPROW(A,J) is equivalent (but more efficient
%  than) the Matlab code
%
%    IDS = find(A(J,:)); 
%    X = A(J,IDS);


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009
