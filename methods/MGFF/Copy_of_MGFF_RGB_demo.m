% This code is the implementation of "Multi-scale Guided Image and Video Fusion: A Fast and Efficient Approach" 
%Cite this article as:
% Bavirisetti, D.P., Xiao, G., Zhao, J. et al. Circuits Syst Signal Process (2019).
%https://doi.org/10.1007/s00034-019-01131-z
 
%% MGFF algorithm to fuse color images
function img=run_MGFF(imgVI, imgIR, outputPath, method, visualization)


    %%
%     clc;
%     clear all;
%     close all;

    % Guided image filter parameters
    r=9;eps=10^3;

    %% load source images
     I1=double(imread('chairs1.jpg'));
     I2=double(imread('chairs5.jpg'));
     
     
    %  % refer my website
    % https://sites.google.com/view/durgaprasadbavirisetti/datasets?authuser=0
    % for all the datasets used in the paper.


    %% apply multi-scale guided image fusion on source images
    tic
    F = fuse_MGF_RGB(I1, I2, r, eps);
    toc
    %% display source images and the fused image
    figure, imshow(uint8(I1), []);
    figure, imshow(uint8(I2),[]);
    figure, imshow((F),[]);
 
end
