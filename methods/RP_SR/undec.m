function Y = undec(X,w);
%Y = undec(X) upsampling of a matrix
%
%    X - input matrix
%    w - choose rows or columns
%        w == 1: upsample rows by 2
%        w == 2: upsample columns by 2
%

%    (Oliver Rockinger 16.08.99)

[z s] = size(X);
if w == 1
  Y = zeros(2*z,s);
  Y(1:2:2*z,1:s) = X;
else
  Y = zeros(z,2*s);
  Y(1:z,1:2:2*s) = X;
end;







