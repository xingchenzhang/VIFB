% This is the code of Hybrid_MSD algorithm:
% 
% The code of Hybrid_MSD is provided by the authors of Hybrid_MSD.
% The interface is created by the authors of VIFB.

function img = run_Hybrid_MSD(imgVI, imgIR, visualization)

    %    The following is an implementation of an infrared and visible image
    %    fusion algorithm proposed in the paper:
    %
    %      "Perceptual fusion of infrared and visible images through a hybrid 
    %      multi-scale decomposition with Gaussian and bilateral filters"
    %      Information Fusion, 2016.
    %    
    %    This code is for testing purpose only.
    %    Some of the test images were obtained at
    %      http://www.imagefusion.org
    %      http://www.ece.lehigh.edu/SPCRL/IF/image_fusion.htm
    %
    %    Zhiqiang Zhou, Beijing Institute of Technology
    %    May, 2015
    
    path_Vis = imgVI.img;
    path_IR  = imgIR.img;
    
    % IR image
    img1 = imread(path_IR);
    
    % VI image
    img2 = imread(path_Vis);    
    
    img1 = double(img1);
    img2 = double(img2);

    if visualization == 1
        paraShow1.fig = 'Visible image';
        paraShow2.fig = 'Infrared image';
        ShowImageGrad(img2, paraShow2);
        ShowImageGrad(img1, paraShow1);
    end
    
    tic;   
    if size(img2, 3) == 1
        fuseimage = Hybrid_MSD(img2, img1);
    elseif size(img1,3) == 1
        fuseimage = zeros(size(img2));
        for i=1:3
            fuseimage(:,:,i) = Hybrid_MSD(img2(:,:,i),img1);    
        end       
    else
        fuseimage = zeros(size(img2));
        for i=1:3
           fuseimage(:,:,i) = Hybrid_MSD(img2(:,:,i),img1(:,:,i));    
        end    
    end       
    toc;
    if visualization == 1
        figure,imshow(fuseimage, []);
    end 
    img = uint8(fuseimage);   
    
end
    
function res = Hybrid_MSD(img1, img2)
    
    nLevel = 4;
    lambda = 30;
    % lambda = 3000;
    
    %% ---------- Hybrid Multi-scale Decomposition --------------
    sigma = 2.0;
    sigma_r = 0.05;
    k = 2;

    M1 = cell(1, nLevel+1);
    M1L = cell(1, nLevel+1);
    M1{1} = img1;
    M1L{1} = img1;
    M1D = cell(1, nLevel+1);
    M1E = cell(1, nLevel+1);
    sigma0 = sigma;
    for j = 2:nLevel+1,
        w = floor(3*sigma0);
        h = fspecial('gaussian', [2*w+1, 2*w+1], sigma0);   
        M1{j} = imfilter(M1{j-1}, h, 'symmetric');
        %M1L{j} = 255*bfilter2(M1L{j-1}/255,w,[sigma0, sigma_r/(k^(j-2))]);
        M1L{j} = 255*fast_bfilter2(M1L{j-1}/255,[sigma0, sigma_r/(k^(j-2))]);

        M1D{j} = M1{j-1} - M1L{j};
        M1E{j} = M1L{j} - M1{j};

        sigma0 = k*sigma0;
    end

    M2 = cell(1, nLevel+1);
    M2L = cell(1, nLevel+1);
    M2{1} = img2;
    M2L{1} = img2;
    M2D = cell(1, nLevel+1);
    M2E = cell(1, nLevel+1);
    sigma0 = sigma;
    for j = 2:nLevel+1,
        w = floor(3*sigma0);
        h = fspecial('gaussian', [2*w+1, 2*w+1], sigma0);   
        M2{j} = imfilter(M2{j-1}, h, 'symmetric');
        %M2L{j} = 255*bfilter2(M2L{j-1}/255,w,[sigma0, sigma_r/(k^(j-2))]);
        M2L{j} = 255*fast_bfilter2(M2L{j-1}/255,[sigma0, sigma_r/(k^(j-2))]);

        M2D{j} = M2{j-1} - M2L{j};
        M2E{j} = M2L{j} - M2{j};

        sigma0 = k*sigma0;
    end

    %% ---------- Multi-scale Combination --------------
    for j = nLevel+1:-1:3
    b2 = abs(M2E{j});
    b1 = abs(M1E{j});
    R_j = max(b2-b1, 0);
    Emax = max(R_j(:));
    P_j = R_j/Emax;

    C_j = atan(lambda*P_j)/atan(lambda);

    % Base level combination
    sigma0 = 2*sigma0;
    if j == nLevel+1
        w = floor(3*sigma0);
        h = fspecial('gaussian', [2*w+1, 2*w+1], sigma0);
        lambda_Base = lambda;
        %lambda_Base = 30;
        C_N = atan(lambda_Base*P_j)/atan(lambda_Base);
        C_N = imfilter(C_N, h, 'symmetric');
        MF = C_N.*M2{j} + (1-C_N).*M1{j};
    end

    % Large-scale combination
    sigma0 = 1.0;
    w = floor(3*sigma0);
    h = fspecial('gaussian', [2*w+1, 2*w+1], sigma0);   
    C_j = imfilter(C_j, h, 'symmetric');

    D_F = C_j.*M2E{j}+ (1-C_j).*M1E{j};
    MF = MF + D_F;
    D_F = C_j.*M2D{j}+ (1-C_j).*M1D{j};
    MF = MF + D_F;
    end 

    % Small-scale combination
    sigma0 = 0.2;
    w = floor(3*sigma0);
    h = fspecial('gaussian', [2*w+1, 2*w+1], sigma0);   
    C_0 = double(abs(M1E{2}) < abs(M2E{2}));
    C_0 = imfilter(C_0, h, 'symmetric');
    D_F = C_0.*M2E{2} + (1-C_0).*M1E{2};
    MF = MF + D_F;  
    C_0 = abs(M1D{2}) < abs(M2D{2});
    C_0 = imfilter(C_0, h, 'symmetric');
    D_F = C_0.*M2D{2} + (1-C_0).*M1D{2};
    MF = MF + D_F;

    %% ---------- Fusion Result --------------
    % FI = ImRegular(MF);   % The intensities are regulated into [0, 255]
    FI = max(min(MF,255), 0);
    res = FI;
end
