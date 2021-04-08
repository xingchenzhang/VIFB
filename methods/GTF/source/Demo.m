clc;
clear;
close all;
warning off;
addpath(genpath(cd));

I=double(imread('E:\lichang\1.Image fusion total variation\img\tank\tank_IR.png'))/255;
%V=double(imread('E:\lichang\1.Image fusion total variation\img\tank\tank_VIS.png'))/255;
% I=histeq(I);
% V=histeq(V);
% figure; imshow(I);
% figure; imshow(V);
% imwrite(I,'new_4917_IR.png','png');
% imwrite(V,'new_4917_VIS.png','png');
% I=rgb2gray(I);
% V=rgb2gray(V);
%proposed
nmpdef;
pars_irn = irntvInputPars('l2tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 5;
pars_irn.U0         = I;

pars_irn.variant       = NMP_TV_SUBSTITUTION;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;
Final_Metric=zeros(6,1); 
tic;
for i=[0.1:0.1:0.9]%[0.01:0.01:0.09 0.1:0.1:0.9 1:10 20:10:100]
   
    X = irntv(I, {}, i, pars_irn);
    X=im2gray(X);
    imwrite(X,strcat('Proposed_tank_',strcat(num2str(i),'.png')),'png');
    %figure; 
    %imshow(X);
%     Result = Metric(uint8(abs(I)*255),uint8(abs(V)*255),uint8(abs(X*255)));
%     temp=Result.Total;
%     Final_Metric=max(Final_Metric,temp);
end
