% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
%
% Reference for the metric:
% Y.-J. Rao, ¡°In-fibre bragg grating sensors,¡± Measurement science and technology, vol. 8, no. 4, p. 355, 1997.

function res  = metricsVariance(img1,img2,fused)

    fused = double(fused); 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);

    if b == 1
        [a,g] = Variance(img1,img2,fused);
        img_var = g;
    elseif b1 == 1
        for k = 1 : b 
           [a,g(k)] = Variance(img1(:,:,k), img2,fused(:,:,k)); 
        end 
        img_var = mean(g); 
    else
        for k = 1 : b 
           [a,g(k)] = Variance(img1(:,:,k), img2(:,:,k),fused(:,:,k)); 
        end 
        img_var = mean(g); 
    end
    res = img_var;
end

function [img_mean,img_var] = Variance(img1,img2,fused)
    if size(fused,3) > 1 
        fused=rgb2gray(fused);  
    end
    fused = double(fused); 
    [r, c] = size(fused); 

    % Mean value 
    img_mean = mean(mean(fused)); 

    % Variance 
    img_var = sqrt(sum(sum((fused - img_mean).^2)) / (r * c ));
end