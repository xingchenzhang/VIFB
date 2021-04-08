function out = Solve_Optimal(M,b1,b2,lambda)

[m, n]=size(b1);

smallNum = 0.000001;

a = abs(imfilter(b2,ones(7)/(sum(sum(ones(7)))),'replicate'));
A = 1./(a + smallNum);

A = A(:);
A = sparse(1:1:m*n,1:1:m*n,A(1:1:m*n)',m*n,m*n);

II=sparse(1:1:m*n,1:1:m*n,ones(1,m*n),m*n,m*n);

denominator=2*II+lambda*(A+A');

numerator=2*M(:)+lambda*(A+A')*b1(:);

out=denominator\numerator;
out=reshape(out,m,n);


end