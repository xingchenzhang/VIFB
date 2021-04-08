function [Y,U] = MSVD(x)
% multiresolution SVD (MSVD)
% input-> x: image (spatial domain)
% outputs-> Y: one level MSVD decomposition of x
%           U: the unitary matrix (U in SVD)

[m,n] = size(x);
m = m/2; n = n/2;
A = zeros(4,m*n);
for j = 1:n
    for i = 1:m
        A(:,i + (j-1)*m) = reshape(x((i-1)*2+(1:2),(j-1)*2+(1:2)),4,1);
    end
end
[U,S] = svd(A);
T = U'*A;
Y.LL = reshape(T(1,:),m,n);
Y.LH = reshape(T(2,:),m,n);
Y.HL = reshape(T(3,:),m,n);
Y.HH = reshape(T(4,:),m,n);