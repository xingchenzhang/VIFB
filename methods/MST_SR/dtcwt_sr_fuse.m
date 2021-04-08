function F=dtcwt_sr_fuse(I1,I2,N,D,overlap,epsilon)
%    DTCWT-SR
%    Input:
%    I1 - input image A
%    I2 - input image B
%    N - maximum decomposition level
%    D  - Dictionary for sparse representation
%    overlap - the overlapped pixels between two neighbor patches
%    epsilon - sparse reconstuction error
%    Output:
%    F  - fused image   
%
%    The code is edited by Yu Liu, 01-09-2014.

I1=double(I1);
I2=double(I2);
% N=2;
[Y1,h1] = dtwavexfm2(I1,N,'legall','qshift_06');
[Y2,h2] = dtwavexfm2(I2,N,'legall','qshift_06');

%Y=(Y1+Y2)/2;
Y=sparse_fusion(Y1,Y2,D,overlap,epsilon);

% [m,n]=size(h1{1,1}(:,:,1));
% [m1,n1]=size(h1{2,1}(:,:,1));
% E1=zeros(m,n);
% E_1=zeros(m1,n1);
% E2=E1;
% E_2=E_1;
% E1=cell(N,1);
% E2=E1;
% W=E1;
% for i=1:N
%     [m,n]=size(h1{i,1}(:,:,1));
%     E1{i}=zeros(m,n);
%     E2{i}=E1{i};
%   for j=1:6
%         E1{i}=E1{i}+abs(h1{i,1}(:,:,j).^2);
%         E2{i}=E2{i}+abs(h2{i,1}(:,:,j).^2);
% %         E_1=E_1+abs(h1{2,1}(:,:,j).^2);
% %         E_2=E_2+abs(h2{2,1}(:,:,j).^2);
%   end
% end
%   W=E1>E2;
%   W1=E_1>E_2;
%   W=medfilt2(W);
%   W1=medfilt2(W1);

h=h1;
for i=1:N
%     W{i}=E1{i}>E2{i};
%     W{i}=medfilt2(W{i});
    for j=1:6
    %   W=abs(h1{i,1}(:,:,j))>abs(h2{i,1}(:,:,j));
        um=3;
        A1 = ordfilt2(abs(es2(h1{i,1}(:,:,j),floor(um/2))), um*um, ones(um));
        A2 = ordfilt2(abs(es2(h2{i,1}(:,:,j),floor(um/2))), um*um, ones(um));
        %second step
        W = (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
        h{i,1}(:,:,j)=h1{i,1}(:,:,j).*W+h2{i,1}(:,:,j).*~W;
    %   h{2,1}(:,:,j)=h1{2,1}(:,:,j).*W1+h2{2,1}(:,:,j).*~W1;
    end
end
F= dtwaveifm2(Y,h,'legall','qshift_06');
  
    


