% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB.
%
% Reference for the metric:
% G. Cui, H. Feng, Z. Xu, Q. Li, and Y. Chen, ¡°Detail preserved fusion of visible and infrared images using regional
% saliency extraction and multi-scale image decomposition,¡± Optics Communications, vol. 341, pp. 199 ¨C 209, 2015

function res = metricsAvg_gradient(img1,img2,fused)
    if nargin == 3 
        fused = double(fused); 
        [r,c,b] = size(fused); 

        dx = 1; 
        dy = 1; 
        for k = 1 : b 
            band = fused(:,:,k); 
            [dzdx,dzdy] = gradient(band,dx,dy); 
            s = sqrt((dzdx .^ 2 + dzdy .^2) ./ 2); 
            g(k) = sum(sum(s)) / ((r - 1) * (c - 1)); 
        end 
        res = mean(g); 
    else 
        error('Wrong number of input!'); 
    end

 