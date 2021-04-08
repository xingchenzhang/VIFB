% This is the main function of the paper "Infrared and Visual Image Fusion 
% through Infrared Feature Extraction and Visual Information, Infrared Physics & Technology, 2017."
% Implemented by Zhang Yu (uzeful@163.com).
%
% The interface is created by the authors of VIFB.

function img = run_IFEVIP(imgVI, imgIR, visualization)
    
    % para settings
    QuadNormDim = 512;
    QuadMinDim = 32;
    GaussScale = 9;
    MaxRatio = 0.001;
    StdRatio = 0.8;
    
    tic;
    
    imgVis = imread(imgVI.img);
    imgIR = imread(imgIR.img);
    
    if size(imgIR,3)==3
        imgIR=rgb2gray(imgIR);
    end

    % image fusion
    result = BGR_Fuse(imgVis, imgIR, QuadNormDim, QuadMinDim, GaussScale, MaxRatio, StdRatio);
    
    % show image
    if visualization == 1
        figure, imshow(result)
    end
    toc;

    img = result;
    
end