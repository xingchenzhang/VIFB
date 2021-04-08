% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB.
%
% Reference for the metric:
% P. Jagalingam and A. V. Hegde, ¡°A review of quality metrics for fused image,¡± Aquatic 
% Procedia, vol. 4, no. Icwrcoe, pp. 133¨C142, 2015

function res = metricsRmse(img1,img2,fused) 
  
    fused = double(fused); 
    % Get the size of img 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);
    img1 = double(img1);
    img2 = double(img2);

    if b == 1
        g = Rmse(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = Rmse(img1(:,:,k), img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = Rmse(img1(:,:,k), img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end

end


function output = Rmse(img1,img2,fused) 

   rmseVF = mse(img1, fused);
   rmseIF = mse(img2, fused);
   
   rmse = rmseVF + rmseIF;   
   output = rmse./2.0;
   
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

