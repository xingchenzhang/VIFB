% This code is in association with the following paper
% "Ma J, Zhou Z, Wang B, et al. Infrared and visible image fusion based on visual saliency map and weighted least square optimization[J].
% Infrared Physics & Technology, 2017, 82:8-17."
% Authors: Jinlei Ma, Zhiqiang Zhou, Bo Wang, Hua Zong
% Code edited by Jinlei Ma, email: majinlei121@163.com
% 
% Interface is created by the author of VIFB.

function img = run_VSMWLS(imgVI, imgIR, visualization)

    % IR image
    I1 = double(imread(imgIR.img))/255;
        
    % VI image
    I2 = double(imread(imgVI.img))/255;
    
    if visualization == 1
        figure;imshow(I1);
        figure;imshow(I2);
    end

    tic;
    if size(I2, 3) == 1
        fused = WLS_Fusion(I1, I2);
    elseif size(I1,3) == 1
        fused = zeros(size(I2));
        for i=1:3
            fused(:,:,i) = WLS_Fusion(I1,I2(:,:,i));    
        end       
    else
        fused = zeros(size(I2));
        for i=1:3
           fused(:,:,i) = WLS_Fusion(I1(:,:,i),I2(:,:,i));    
        end    
    end     
    toc;
    
    if visualization == 1
        figure;imshow(fused);
    end
    
    img = im2uint8(fused);
end