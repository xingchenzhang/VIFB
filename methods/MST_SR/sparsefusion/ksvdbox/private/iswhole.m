function z = iswhole(x,epsilon)
%ISWHOLE Determine whole numbers with finite precision.
%  Z = ISWHOLE(X,EPSILON) returns a matrix the same size as X, containing
%  1's where X is whole up to precision EPSILON, and 0's elsewhere. 
%
%  Z = ISWHOLE(X) uses the default value of EPSILON=1e-6.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  May 2005

if (nargin<2), epsilon = 1e-6; end
z = abs(round(x)-x) < epsilon;
