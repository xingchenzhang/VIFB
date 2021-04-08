% This code is the implementation of "Multi-scale Guided Image and Video Fusion: A Fast and Efficient Approach" 
%Cite this article as:
% Bavirisetti, D.P., Xiao, G., Zhao, J. et al. Circuits Syst Signal Process (2019).
%https://doi.org/10.1007/s00034-019-01131-z
 

%%
 clc;
clear all;
close all;

% Guided image filter parameters
r=9;eps=10^3;

         

%% load source images
 I1=double(imread('source24_1.tif'));
 I2=double(imread('source24_2.tif'));
%  % refer my website
% https://sites.google.com/view/durgaprasadbavirisetti/datasets?authuser=0
% for all the datasets used in the paper.


%% apply multi-scale guided image fusion on source images
tic
F = fuse_MGF(I1, I2, r, eps);
toc
%% display source images and the fused image
figure, imshow((I1), []);
figure, imshow((I2),[]);
figure, imshow((F),[]);
%% apply image contrast enhancement techniques to further enhance the fused image quality
% we either implemented them or used built-in matlab functions.
%They are
%1) Histogram Equalization (HE),
%2) (Bi-histogram equalization) BHE, 
%3) Recursive mean-separate histogram equalization (RMSHE), 
%4) Bi-bi-histogram equalization with variable enhancement degree (BBHEwVED), 
%5) Recursively separated and weighted histogram equalization (RSWHE) 
%6) Contrast Limited Adaptive Histogram Equalization (CLAHE).
%References:
%[1]  R. Gonzalez and R.Woods, “The book”, Digital Image Processing, 2002. 
%[2] Y. Kim, “Contrast enhancement using brightness preserving bi-histogram equalization”, IEEE Transactions on Consumer Electronics, vol. 43, no 1, pp 1-8, 2002.
%[3] S. Chen and A. Ramli, “Contrast enhancement using recursive mean-separate histogram equalization for scalable brightness preservation,” IEEE Transactions on Consumer Electronics, vol. 49, no. 4, pp. 1301-1309, 2003.
%[4] K. Murahira, T.Kawakami and A. Taguchi, “Modified histogram equalization for image contrast enhancement”, in 4th International Symposium on Communications, Control and Signal Processing (ISCCSP), 2010, IEEE, pp. 1-5.
%[5] M. Kim and M. Chung, “Recursively separated and weighted histogram equalization for brightness preservation and contrast enhancement”, IEEE Transactions on Consumer Electronics, vol. 54, no.3, pp. 1389-1397, 2008.
%[6] Zuiderveld, Karel. "Contrast Limited Adaptive Histogram Equalization." Graphic Gems IV. San Diego: Academic Press Professional, 1994. 474–485.


 % Contrast Limited Adaptive Histogram Equalization (CLAHE). (MATLAB built-in function)
  F1 = adapthisteq(uint8(F));
 
% Histogram Equalization (HE) (MATLAB built-in function)
 F2 = histeq(uint8(F));
 % Bi-histogram equalization) BHE (Implemented). 
 F3 = BHEg_gray(uint8(F));
% Recursive mean-separate histogram equalization (RMSHE) (Implemented)
 F4= RMSHEg_gray(uint8(F));
 %Bi-bi-histogram equalization with variable enhancement degree (BBHEwVED) (Implemented)
 F5 = BBHEwVEDg_gray(uint8(F));
 % Recursively separated and weighted histogram equalization (RSWHE) (Implemented)
 F6 = RSWHEg_gray(uint8(F));

 figure(4), imshow((F1),[]);
 figure(5), imshow((F2),[]);
 figure(6), imshow((F4),[]);
 figure(7), imshow((F5),[]);
 figure(8), imshow((F6),[]);
 figure(9), imshow((F3),[]);

%This is the basic implementation of the MGFF algorithm. 
 % I have also implemented codes to fuse RGB and sequence of images. I can
 % send them upon your request(bdps1989@gmail.com)
  
