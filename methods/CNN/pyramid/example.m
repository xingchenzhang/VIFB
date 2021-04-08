function example;

clear;
clc;
close all;

%dyramic,mode=1
%I=load_images('..\sourceimages\ArchSequence\',1);
%I=load_images('..\sourceimages\ForrestSequence\',1);
I=load_images('..\sourceimages\campus\',1);

%static,mode=0

%I=load_images('..\sourceimages\memorial\',1);
%I=load_images('..\sourceimages\BelgiumHouse\',1);
%I=load_images('..\sourceimages\garage\',1);


%I=load_images('..\sourceimages\flash\',1);
tic;
R = exposure_fusion(I,[1 1 1]);
toc;
figure('Name','Result'); 
imshow(R); 

imwrite(R,'campus_ef.jpg')
