function A = extract_l1_feature(out)
% extract l1-norm features
[m, n, k] = size(out);
A_temp = zeros(m+2,n+2);

L1_norm = sum(abs(out),3);

A_temp(2:m+1, 2:n+1) = double(L1_norm);

A = A_temp;
end