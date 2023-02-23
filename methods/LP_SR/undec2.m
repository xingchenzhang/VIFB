function Y = undec2(X)
%Y = undec2(X) upsampling of a matrix by 2
%
%    X - input matrix
%
%    Y - output matrix

%    (Oliver Rockinger 16.08.99)

[z s] = size(X);
Y     = zeros(2*z, 2*s); 

Y(1:2:2*z,1:2:2*s) = X;







