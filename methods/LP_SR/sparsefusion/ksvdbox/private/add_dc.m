function x = add_dc(y,dc,columns)
%ADD_DC Add DC channel to signals.
%   X = ADD_DC(Y,DC) adds the specified DC value to the (possibly
%   multi-dimensional) signal Y, returning the result as X. DC should be a
%   scalar value.
%
%   X = ADD_DC(Y,DC,'columns') where Y is a 2D matrix and DC is an array of
%   length size(Y,2), treats the columns of Y as individual 1D signals, 
%   adding to each one the corresponding DC value from the DC array. X is
%   the same size as Y and contains the resulting signals.
%
%   See also REMOVE_DC.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009


if (nargin==3 && strcmpi(columns,'columns')), columns = 1;
else columns = 0;
end

if (columns)
  x = addtocols(y,dc);
else
  x = y + dc;
end



