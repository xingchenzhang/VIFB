function y1=curvelet_sr_fuse(I1,I2,n,D,overlap,epsilon)
%    CVT-SR
%    Input:
%    I1 - input image A
%    I2 - input image B
%    n - maximum decomposition level
%    D  - Dictionary for sparse representation
%    overlap - the overlapped pixels between two neighbor patches
%    epsilon - sparse reconstuction error
%    Output:
%    y1  - fused image   
%
%    The code is edited by Yu Liu, 01-09-2014.

I1=double(I1);
I2=double(I2);

% n=3;
is_real = 1;

finest = 1 ;
nbscales=n;
C1 = fdct_wrapping(I1, is_real, finest, nbscales);
C2 = fdct_wrapping(I2, is_real, finest, nbscales);
[M,N]=size(I1);
% C1 = fdct_usfft(I1,is_real,n);
% C2 = fdct_usfft(I2,is_real,n);

 C=C1;
%  se=[4 4 4 4 4;4 16 16 16 4;4 16 64 16 4;4 16 16 16 4;4 4 4 4 4]/256;
for scalelevel=2:n
    for orientationlevel=1:length(C1{scalelevel})
        E1=C1{scalelevel}{orientationlevel};
        E2=C2{scalelevel}{orientationlevel};
%         E1=conv2(E1,se,'same');
%         E2=conv2(E2,se,'same');
%         W=abs(E1)>abs(E2);
    um=3;
    A1 = ordfilt2(abs(es2(E1,floor(um/2))), um*um, ones(um));
  	A2 = ordfilt2(abs(es2(E2,floor(um/2))), um*um, ones(um));
    % second step
  	W = (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
        C{scalelevel}{orientationlevel}=C1{scalelevel}{orientationlevel}.*W+C2{scalelevel}{orientationlevel}.*~W;
    end
end
%C{1}{1}=(C1{1}{1}+C2{1}{1})/2;
C{1}{1}=sparse_fusion(C1{1}{1},C2{1}{1},D,overlap,epsilon);
% y1 = ifdct_usfft(C,is_real);


y1 = ifdct_wrapping(C, is_real, M, N);



