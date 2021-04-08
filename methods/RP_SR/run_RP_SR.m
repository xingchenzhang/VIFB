% This is the code for PR_SR image fusion algorithm.
% Y. Liu, S. Liu, and Z. Wang, ¡°A general framework for image fusion based on multi-scale transform and sparse representation,¡± Information 
% Fusion, vol. 24, pp. 147¨C164, 2015.
% 
% The interface is created by the authors of VIFB.

function img = run_RP_SR(imgVI, imgIR, visualization)

    addpath(genpath('sparsefusion'));
    addpath(genpath('dtcwt_toolbox'));
    addpath(genpath('fdct_wrapping_matlab'));
    addpath(genpath('nsct_toolbox'));

    load('sparsefusion/Dictionary/D_100000_256_8.mat');
    
    % IR image
    image_input1 = imread(imgIR.img);

    % VI image
    image_input2 = imread(imgVI.img);
     
    if visualization == 1
        figure;imshow(image_input1);
        figure;imshow(image_input2);
    end

    img1=double(image_input1);
    img2=double(image_input2);

    overlap = 6;                    
    epsilon=0.1;
    level=4;

    tic;    
    if size(img2,3)==1   
        imgf = rp_sr_fuse(img1,img2,level,3,3,D,overlap,epsilon);    
    elseif size(img1,3) == 1    
        imgf = zeros(size(img2)); 
        for i=1:3
            imgf(:,:,i) = rp_sr_fuse(img1,img2(:,:,i),level,3,3,D,overlap,epsilon);    
        end    
    else
        imgf = zeros(size(img2)); 
        for i=1:3
            imgf(:,:,i) = rp_sr_fuse(img1(:,:,i),img2(:,:,i),level,3,3,D,overlap,epsilon);    
        end    
    end
    toc;

    if visualization == 1
        figure;imshow(uint8(imgf));
    end
    img = uint8(imgf);
end

