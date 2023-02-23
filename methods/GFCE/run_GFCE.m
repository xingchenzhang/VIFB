% This is the code for GFCE image fusion algorithms:
% Z. Zhou, M. Dong, X. Xie, and Z. Gao, ¡°Fusion of infrared and visible images for night-vision context enhancement,¡±
% Applied optics, vol. 55, no. 23, pp. 6480¨C6490, 2016
%
% The code of GFCE is provided by the authors of GFCE.
% The interface is created by the authors of VIFB. Necessary modifications
% are made to integrated GFCE to VIFB.

function img = run_GFCE(imgVI, imgIR, visualization)

    %    The following is an implementation of the guided filter (GF) based
    %    context enhancement (GFCE) through fusion of infrared and visible
    %    images.
    %    
    %    Some of the test images were obtained at
    %      http://www.imagefusion.org
    %      http://www.ece.lehigh.edu/SPCRL/IF/image_fusion.htm
    %
    %    Zhiqiang Zhou, Beijing Institute of Technology
    %    Apr. 2016
    
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
        fuseimage = GFCE(img2, img1);
    elseif size(img1,3) == 1
        fuseimage = zeros(size(img2));
        for i=1:3
            fuseimage(:,:,i) = GFCE(img2(:,:,i),img1);    
        end       
    else
        fuseimage = zeros(size(img2));
        for i=1:3
           fuseimage(:,:,i) = GFCE(img2(:,:,i),img1(:,:,i));    
        end    
    end   
    toc;
    if visualization == 1
        figure,imshow(fuseimage, []);
    end 
    img = uint8(fuseimage);        
end   
    
function res = GFCE(img1, img2)

    nLevel = 4;
    
    %% ---------- Visibility enhancement for visible image--------------
    img1E = Ehn_GF(img1);
    img1 = img1E;
    %% ---------- Infrared image normalization--------------
    mi = min(img2(:));
    ma = max(img2(:));
    img2 = (img2-mi)/(ma-mi)*255;
    %% ---------- Automatic parameter selection --------------
    Rs = Relative_PS(img2, img1);
    if Rs<0.8
        lambda = 100;
    else
        if Rs>1.6
            lambda = 2000;
        else
            lambda = 2500*Rs - 1900;
        end
    end    

    %% ---------- Hybrid multiscale decomposition based on guided filter--------------
    sigma = 2;  k = 2;
    r0 = 2;     eps0 = 0.1;  
    l = 2;

    M1 = cell(1, nLevel+1);
    M1L = cell(1, nLevel+1);
    M1{1} = img1/255;
    M1L{1} = M1{1};
    M1D = cell(1, nLevel+1);
    M1E = cell(1, nLevel+1);
    sigma0 = sigma;
    r = r0;
    eps = eps0;
    for ii = 2:nLevel+1,

    %     % ***using fast guided filter, which has the potential to achieve real-time performance when codes are fully optimized
    %     % ***NOTE: large subsampling ratio may cause problem for fusion of some source images.
    %     s = max(1, r/2); % subsampling ratio
    %     M1{ii} = fastguidedfilter_md(M1{ii-1}, M1{ii-1}, r, 100^2, s);  
    %     M1L{ii} = fastguidedfilter_md(M1L{ii-1}, M1L{ii-1}, r, eps^2, s);

        M1{ii} = guidedfilter(M1{ii-1}, M1{ii-1}, r, 100^2);  
        M1L{ii} = guidedfilter(M1L{ii-1}, M1L{ii-1}, r, eps^2);    

        M1D{ii} = M1{ii-1} - M1L{ii};
        M1E{ii} = M1L{ii} - M1{ii};

        sigma0 = k*sigma0;
        r = k*r;
        eps = eps/l;
    end

    M2 = cell(1, nLevel+1);
    M2L = cell(1, nLevel+1);
    M2{1} = img2/255;
    M2L{1} = M2{1};
    M2D = cell(1, nLevel+1);
    M2E = cell(1, nLevel+1);
    sigma0 = sigma;
    r = r0;
    eps = eps0;
    for ii = 2:nLevel+1,
    %     s = max(1, r/2);
    %     M2{ii} = fastguidedfilter_md(M2{ii-1}, M2{ii-1}, r, 100^2, s);
    %     M2L{ii} = fastguidedfilter_md(M2L{ii-1}, M2L{ii-1}, r, eps^2, s);
        M2{ii} = guidedfilter(M2{ii-1}, M2{ii-1}, r, 100^2);
        M2L{ii} = guidedfilter(M2L{ii-1}, M2L{ii-1}, r, eps^2);    

        M2D{ii} = M2{ii-1} - M2L{ii};
        M2E{ii} = M2L{ii} - M2{ii};

        sigma0 = k*sigma0;
        r = k*r;
        eps = eps/l;
    end

    %% ---------- Fusion --------------

    for j = nLevel+1:-1:3
    D2 = abs(M2E{j});
    D1 = abs(M1E{j});
    R = max(D2-D1, 0);
    Rmax = max(R(:));
    P = R/Rmax;

    Cj = atan(lambda*P)/atan(lambda);

    sigma_b = 2*sigma0;
    if j == nLevel+1
        w = floor(3*sigma_b);
        h = fspecial('gaussian', [2*w+1, 2*w+1], sigma_b);
        lambda0 = lambda;
        Cb = atan(lambda0*P)/atan(lambda0);
        Cb = imfilter(Cb, h, 'symmetric');
        MB = Cb.*M2{nLevel+1} + (1-Cb).*M1{nLevel+1};
    end

    sigma_c = 1;
    w = floor(3*sigma_c);
    h = fspecial('gaussian', [2*w+1, 2*w+1], sigma_c);   
    Cj = imfilter(Cj, h, 'symmetric');

    md = Cj.*M2E{j}+ (1-Cj).*M1E{j};
    MB = MB + md;
    md = Cj.*M2D{j}+ (1-Cj).*M1D{j};
    MB = MB + md;
    end 

    sigma_t = 1;
    w = floor(3*sigma_t);
    h = fspecial('gaussian', [2*w+1, 2*w+1], sigma_t);   
    C11 = double(abs(M1E{2}) < abs(M2E{2}));
    C11 = imfilter(C11, h, 'symmetric');
    md = C11.*M2E{2}+ (1-C11).*M1E{2};
    MB = MB + md;  
    C10 = double(abs(M1D{2}) < abs(M2D{2}));
    md = C10.*M2D{2}+ (1-C10).*M1D{2};
    MB = MB + md;
    FI = min(round(MB*275), 255);
    FI = max(FI, 0);
    res = FI;
end

