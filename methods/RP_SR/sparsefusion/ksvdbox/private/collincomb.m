%COLLINCOMB Linear combination of matrix columns.
%  Y = COLLINCOMB(A,COLS,X) computes a linear combination of the columns of
%  the matrix A. The column indices are specified in the vector COLS, and
%  the correspoinding coefficients are specified in the vector X. The
%  vectors COLS and X must be of the same length. 
%  
%  The call Y = COLLINCOMB(A,COLS,X) is essentially equivalent to
%
%         Y = A(:,COLS)*X .
%
%  However, it is implemented much more efficiently.
%
%  Y = COLLINCOMB(A,ROWS,COLS,X) only works on the rows of A specified
%  in ROWS, returning a vector of length equal to ROWS. This call is
%  essentially equivalent to the command
%
%         Y = A(ROWS,COLS)*X .
%
%  See also ROWLINCOMB.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009
