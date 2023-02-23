%ROWLINCOMB Linear combination of matrix rows.
%  Y = ROWLINCOMB(X,A,ROWS) computes a linear combination of the rows of
%  the matrix A. The row indices are specified in the vector ROWS, and the
%  correspoinding coefficients are specified in the vector X. The vectors
%  ROWS and X must be of the same length. The call Y = ROWLINCOMB(X,A,ROWS)
%  is essentially equivalent to the command
%
%         Y = X'*A(ROWS,:) .
%
%  However, it is implemented much more efficiently.
%
%  Y = ROWLINCOMB(X,A,ROWS,COLS) only works on the columns of A specified
%  in COLS, returning a vector of length equal to COLS. This call is
%  essentially equivalent to the command
%
%         Y = X'*A(ROWS,COLS) .
%
%  See also COLLINCOMB.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009
