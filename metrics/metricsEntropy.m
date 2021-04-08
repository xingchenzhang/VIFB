% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
%
% V. Aardt and Jan, ¡°Assessment of image fusion procedures using entropy, image quality, and multispectral 
% classification,¡± Journal of Applied Remote Sensing, vol. 2, no. 1, p.023522, 2008.

function res = metricsEntropy(img1,img2,fused) 
    
    fused = double(fused); 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);

    if b == 1
        g = Entropy(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = Entropy(img1(:,:,k),img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = Entropy(img1(:,:,k),img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end

end

function output = Entropy(img1,img2,fused) 
    h = fused;
    s=size(size(h));
    if s(2)==3 
        h1=rgb2gray(h);
    else
        h1=h;
    end    
    h1=double(h1);
    [m,n]=size(h1);
    X=zeros(1,256);
    result=0;
    for i=1:m
        for j=1:n
            X(floor(h1(i,j))+1)=X(floor(h1(i,j))+1)+1;
        end
    end
    for k=1:256
        P(k)=X(k)/(m*n);
        if (P(k)~=0)
            result=P(k)*log2(P(k))+result;
        end
    end
    result=-result;
    EN1=result;
    output = EN1;
end

