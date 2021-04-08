% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
%
% Reference for the metric:
% P. Jagalingam and A. V. Hegde, ¡°A review of quality metrics for fused image,¡± Aquatic Procedia, vol. 4, 
% no. Icwrcoe, pp. 133¨C142, 2015.

function res = metricsPsnr(img1,img2,fused)
   
    fused = double(fused); 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);
    img1 = double(img1);
    img2 = double(img2);

    if b == 1
        g = Psnr(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = Psnr(img1(:,:,k), img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = Psnr(img1(:,:,k), img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end

end

function PSNR = Psnr(img1,img2,fused)

    B=8;               
    MAX=2^B-1;        

    MES = (mse(img1, fused) + mse(img2, fused))./2.0;
    PSNR=20*log10(MAX/sqrt(MES));
    
end

function res0 = mse(a, b)
    if size(a,3) > 1
        a = rgb2gray(a);  
    end

    if size(b,3) > 1
        b = rgb2gray(b); 
    end

    [m, n]=size(a);
    temp=sqrt(sum(sum((a-b).^2)));
    res0=temp/(m*n);
end
