 clc;
 clear all;
I1=imread('E:\sun\ͼƬ\clockA.bmp');
I2=imread('E:\sun\ͼƬ\clockB.bmp');
%  I1=colorspace('yuv<-rgb',I1);
% % %I2=colorspace('yuv<-rgb',I2);
%  I1=I1(:,:,1);
% %Y2=I1(:,:,2);
% 
% % %I1=imread('D:\MATLAB701\work\Lena.jpg');
% % I1=imread('E:\sun\ͼƬ\clockA.bmp');
% % I2=imread('E:\sun\ͼƬ\clockB.bmp');
I1=double(I1);
I2=double(I2);
I=(I1+I2)/2;        %??
n=1;
[Y1,h1] = dtwavexfm2(I1,n,'near_sym_b','qshift_b');
[Y2,h2] = dtwavexfm2(I2,n,'near_sym_b','qshift_b');
[Y,h] = dtwavexfm2(I,n,'near_sym_b','qshift_b');


[r,c]=size(h1{:,:,1});
for i=1:6
   h_seg(:,:,i)=seg(h{1,1}(:,:,i));%??
   F=zeros(r,c);
  L2=bwlabel(h_seg(:,:,i));
  ma=max(max(L2));
      for k=1:ma
    I=L2==k;
     index1=find(h_seg(:,:,i)==1);
     F(index1)=W(i,k)*I1(index)+(1-W(i,k))*I2(index1);
      end
  Fu(:,:,i)=F;
end
Fuse= dtwaveifm2(Yl,Fu,'near_sym_b','qshift_b');
h=Entrophy(Fuse);











