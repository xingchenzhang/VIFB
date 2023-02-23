function Y = dec(X,w);
%Y = dec(X) downsampling of a matrix
%
%    X - input matrix
%    w - choose rows or columns
%        w == 1: downsample rows by 2
%        w == 2: downsample columns by 2
%
%    Y - output matrix

%    (Oliver Rockinger 16.08.99)

[a b] = size(X);
if w == 1
  Y = X(1:2:a, :);
else
  Y = X(:, 1:2:b);
end;
   
