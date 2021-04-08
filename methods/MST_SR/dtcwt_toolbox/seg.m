%function I=seg(I1)


 clc;
 clear all;
%  I1=imread('E:\sun\Í¼Æ¬\Defocused2.jpg');
% % %I2=imread('E:\sun\Í¼Æ¬\fig2.jpg');
%  I1=colorspace('yuv<-rgb',I1);
% % %I2=colorspace('yuv<-rgb',I2);
%  I1=I1(:,:,1);
% %Y2=I1(:,:,2);
% 
I1=imread('D:\MATLAB701\work\Lena.jpg');
% % I1=imread('E:\sun\Í¼Æ¬\clockA.bmp');
% % I2=imread('E:\sun\Í¼Æ¬\clockB.bmp');
I1=double(I1);
[R,C]=size(I1);
[ax,ay]=gfilter(I1,2);
IG=((ax.^2)+(ay.^2)).^0.5;
n=1;
[Yl,Yh] = dtwavexfm2(I1,n,'near_sym_b','qshift_b');
%[r,c]=size(Yh{1,1}(:,:,1));
for i=1:n
   for j=1:6
A{i,1}(:,:,j)=abs(Yh{i,1}(:,:,j));
   end
end
%ÖÐÖµÂË²¨
for i=1:n
    X=7+2*i;
S1{i,1}(:,:,1)=I_medfilt(A{i,1}(:,:,1),X,105);
S{i,1}(:,:,1)=I_medfilt(S1{i,1}(:,:,1),X,15);

S1{i,1}(:,:,2)=I_medfilt(A{i,1}(:,:,2),X,75);
S{i,1}(:,:,2)=I_medfilt(S1{i,1}(:,:,2),X,-15);

S1{i,1}(:,:,3)=I_medfilt(A{i,1}(:,:,3),X,165);
S{i,1}(:,:,3)=I_medfilt(S1{i,1}(:,:,3),X,75);

S1{i,1}(:,:,4)=I_medfilt(A{i,1}(:,:,4),X,15);
S{i,1}(:,:,4)=I_medfilt(S1{i,1}(:,:,4),X,-75);

S1{i,1}(:,:,5)=I_medfilt(A{i,1}(:,:,5),X,135);
S{i,1}(:,:,5)=I_medfilt(S1{i,1}(:,:,5),X,45);

S1{i,1}(:,:,6)=I_medfilt(A{i,1}(:,:,6),X,45);
S{i,1}(:,:,6)=I_medfilt(S1{i,1}(:,:,6),X,-45);
end
se1 = strel('square',3);
for k=1:n
   for i=1:6
[ax,ay]=gfilter(S{k,1}(:,:,i),2);
TG{k,1}(:,:,i)=((ax.^2)+(ay.^2)).^0.5;                                     
M{k,1}(:,:,i)=imerode(S{k,1}(:,:,i)./(2^n),se1);% ÐÎÌ¬Ñ§¸¯Ê´
   end
end
for k=1:n
[r,c]=size(TG{k,1}(:,:,1));
 N=r*c*6;
   for i=1:6
    Max=max(max(TG{k,1}(:,:,i)));
        T_G{k,1}(:,:,i)=TG{k,1}(:,:,i)/Max;
W{k,1}(:,:,i)=N/sum(sum(T_G{k,1}(:,:,i).^2)); 
   end
end

TG1=zeros(R,C);
E=zeros(R,C);
for k=1:n
  for i=1:6
    TGG{k,1}(:,:,i)=W{k,1}(:,:,i).*T_G{k}(:,:,i);
    TGg{k,1}(:,:,i)=imresize(TGG{k,1}(:,:,i),[R,C],'bilinear');%%²åÖµ
    TG1=TG1+TGg{k,1}(:,:,i);
    
     ME{k,1}(:,:,i)=imresize(M{k,1}(:,:,i),[R,C],'bilinear');
      E=E+ME{k,1}(:,:,i);
  end
end
    
WT=median(TG1(:));
WI=3*median(IG(:));
   
% TG1=sort(TG2(:));
% WT=4*median(TG1);
% IG1=sort(IG(:));
% WI=median(IG1);
a=0.8;
b=7;
for i=1:R
    for j=1:C
        h=E(i,j)/a-b;
        if h<0
            k=0;
        else
            k=h;
        end
        Activity(i,j)=exp(k);
        I_G(i,j)=IG(i,j)/Activity(i,j);
        GS(i,j)=I_G(i,j)/WI+TG1(i,j)/WT;
        %GS(i,j)=IG(i,j)/(Activity(i,j)*WI)+TG1(i,j)/WT;   
    end
end

TH=0.45*median(GS(:));

I3 = imhmin(GS,TH); 
L = watershed(I3);
I=L~=0;
figure;imshow(mat2gray(I))
%I_seg=spectral(L,I1)







