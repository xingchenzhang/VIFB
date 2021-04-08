function[x] = IMSVD(Y,U)
% inverse MSVD (IMSVD)
% Inputs-> Y: MSVD coefficients & U: unitary matrix (U in SVD)
% output-> x: image (spaitial domain)

[m,n] = size(Y.LL);
mn = m*n;
T = zeros(4,mn);
T(1,:) = reshape(Y.LL,1,mn);
T(2,:) = reshape(Y.LH,1,mn);
T(3,:) = reshape(Y.HL,1,mn);
T(4,:) = reshape(Y.HH,1,mn);
A = U * T;  
x = zeros(m*2,n*2); 
for j = 1:n
    for i = 1:m
        x((i-1)*2+(1:2), (j-1)*2+(1:2)) = reshape(A(:,i+(j-1)*m),2,2);
    end
end