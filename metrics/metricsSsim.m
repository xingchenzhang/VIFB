% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
%
% Reference for the metric:
% Z. Wang, A. C. Bovik, H. R. Sheikh, E. P. Simoncelli et al.,“Image quality assessment: from error visibility 
% to structural similarity,” IEEE transactions on image processing, vol. 13, no. 4, pp. 600–612, 2004.

function res = metricsSsim(img1,img2,fused)  

    fused = double(fused); 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);

    if b == 1
        g = Ssim(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = Ssim(img1(:,:,k), img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = Ssim(img1(:,:,k), img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end

end


function output = Ssim(img1,img2,fused)  

    ssimVI = ssim(img1, fused);
    ssimIR = ssim(img2, fused);
    output = ssimVI + ssimIR;
    
end

function mssim = ssim(img1,fused)  
    % ========================================================================  
    %ssim的算法主要参考如下论文：  
    %Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image  
    % quality assessment: From error visibility to structural similarity,"  
    % IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612,  
    % Apr. 2004.  
    %  首先对图像加窗处理，w=fspecial('gaussian', 11, 1.5);  
    %                 (2*ua*ub+C1)*(2*sigmaa*sigmab+C2)  
    %   SSIM(A,B)=————————————————————————  
    %              (ua*ua+ub*ub+C1)(sigmaa*sigmaa+sigmab*sigmab+C2)  
    %     C1=（K1*L）;  
    %     C2=(K2*L);   K1=0.01，K2=0.03  
    %     L为灰度级数，L=255  
    %-------------------------------------------------------------------  
    %     ima - 比较图像A  
    %     imb - 比较图像B  
    %  
    % ssim_map - 各加窗后得到的SSIM（A,B|w）组成的映射矩阵  
    %    mssim - 对加窗得到的SSIM（A,B|w）求平均，即最终的SSIM（A,B）  
    %  siga_sq - 图像A各窗口内灰度值的方差  
    %  sigb_sq - 图像B各窗口内灰度值的方差  
    %-------------------------------------------------------------------  
    %  Cool_ben  
    %========================================================================  

    w = fspecial('gaussian', 11, 1.5);  %window 加窗  
    K(1) = 0.01;                      
    K(2) = 0.03;                      
    L = 255;       

    s=size(size(img1));
    if s(2)==3 %判断是灰度图还是RGB
        img1=rgb2gray(img1);
    end 

    s1=size(size(fused));
    if s1(2)==3 %判断是灰度图还是RGB
        fused=rgb2gray(fused);
    end 

    img1 = double(img1);  
    fused = double(fused);  

    C1 = (K(1)*L)^2;  
    C2 = (K(2)*L)^2;  
    w = w/sum(sum(w));  

    ua   = filter2(w, img1, 'valid');%对窗口内并没有进行平均处理，而是与高斯卷积，  
    ub   = filter2(w, fused, 'valid'); % 类似加权平均  
    ua_sq = ua.*ua;  
    ub_sq = ub.*ub;  
    ua_ub = ua.*ub;  
    siga_sq = filter2(w, img1.*img1, 'valid') - ua_sq;  
    sigb_sq = filter2(w, fused.*fused, 'valid') - ub_sq;  
    sigab = filter2(w, img1.*fused, 'valid') - ua_ub;  

    ssim_map = ((2*ua_ub + C1).*(2*sigab + C2))./((ua_sq + ub_sq + C1).*(siga_sq + sigb_sq + C2));     
    mssim = mean2(ssim_map);  

    %return mssim 
end