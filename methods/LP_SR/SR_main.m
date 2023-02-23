clear all;
close all;
clc;

addpath(genpath('sparsefusion'));
load('Dictionary/D_100000_256_8.mat');

[imagename1 imagepath1]=uigetfile('source_images\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the first input image');
image_input1=imread(strcat(imagepath1,imagename1));    
[imagename2 imagepath2]=uigetfile('source_images\*.jpg;*.bmp;*.png;*.tif;*.tiff;*.pgm;*.gif','Please choose the second input image');
image_input2=imread(strcat(imagepath2,imagename2));    

figure;imshow(image_input1);
figure;imshow(image_input2);

if size(image_input1)~=size(image_input2)
    error('two images are not the same size.');
end

img1=double(image_input1);
img2=double(image_input2);

overlap = 6;                    
epsilon=0.1;

tic;
imgf=sparse_fusion(img1,img2,D,overlap,epsilon);
toc

image_fusion=uint8(imgf);
figure;imshow(image_fusion);
imwrite(image_fusion,'Results/fused_sr.tif');
