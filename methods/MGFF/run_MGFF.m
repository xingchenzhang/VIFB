% This code is the implementation of "Multi-scale Guided Image and Video Fusion: A Fast and Efficient Approach" 
% Cite this article as:
% Bavirisetti, D.P., Xiao, G., Zhao, J. et al. Circuits Syst Signal Process (2019).
%https://doi.org/10.1007/s00034-019-01131-z
% 
% The interface is created by the authors of VIFB.
 
function img = run_MGFF(imgVI, imgIR, visualization)

    % Guided image filter parameters
    r=9;eps=10^3;
    
    I1 = double(imread(imgIR.img));
    if size(I1,3)==1
        I1 = repmat(I1,[1,1,3]);
    end
    
    I2 = double(imread(imgVI.img));     

    %% apply multi-scale guided image fusion on source images
    tic
    F = fuse_MGF_RGB(I1, I2, r, eps);
    toc
    
    fuseimage = im2uint8(F);
    img = fuseimage;
    
    %% display source images and the fused image
    if visualization == 1
        figure, imshow(uint8(I1), []);
        figure, imshow(uint8(I2),[]);
        figure, imshow((F),[]);
    end 
end
