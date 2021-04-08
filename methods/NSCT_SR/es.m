function Y = es(X, n, wo)
%Y = ES2(X, n) symmetric extension of a matrix on selected borders
%    
%    X  - input matrix
%    n  - number of rows/columns to extend
%    wo - where to extend
%       wo == 1: left and right
%				wo == 2: up and bottom
%				wo == 3: all borders
%
%    Y  - extended matrix

%    (Oliver Rockinger 16.08.99)

[z s] = size(X);

if wo == 1
  Y = zeros(z, s+2*n);
  Y(:,n:-1:1)  = X(:,2:1:n+1); 
  Y(:,n+1:1:n+s) = X(:,:);
  Y(:,n+s+1:1:s+2*n) = X(:,s-1:-1:s-n); 
end;

if wo == 2
  Y = zeros(z+2*n, s);
  Y(n:-1:1,:)  = X(2:1:n+1,:); 
  Y(n+1:n+z,:) = X(:,:);
  Y(n+z+1:1:z+2*n,:) = X(z-1:-1:z-n,:); 
end;

if wo == 3
  Y = zeros(z+2*n, s+2*n);
  Y(n+1:n+z,n:-1:1)  = X(:,2:1:n+1); 
  Y(n+1:n+z,n+1:1:n+s) = X;
  Y(n+1:n+z,n+s+1:1:s+2*n) = X(:,s-1:-1:s-n);
  Y(n:-1:1,n+1:s+n)  = X(2:1:n+1,:); 
  Y(n+1:n+z,n+1:s+n) = X;
  Y(n+z+1:1:z+2*n,n+1:s+n) = X(s-1:-1:s-n,:); 
end; 



