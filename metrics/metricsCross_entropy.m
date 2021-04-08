% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
% 
% Reference for the metric:
% D. M. Bulanon, T. Burks, and V. Alchanatis, ¡°Image fusion of visible and thermal images for fruit detection,¡±
% Biosystems Engineering, vol. 103, no. 1, pp. 12¨C22, 2009.

function res = metricsCross_entropy(img1,img2,fused)

    fused = double(fused); 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);

    if b == 1
        g = cross_entropy_single(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = cross_entropy_single(img1(:,:,k), img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = cross_entropy_single(img1(:,:,k), img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end

end

function output = cross_entropy_single(img1,img2,fused)

    cross_entropyVI = cross_entropy(img1,fused);
    cross_entropyIR = cross_entropy(img2,fused);
    output = (cross_entropyVI + cross_entropyIR)./2.0;

end

function res0 = cross_entropy(img1,fused)
    s=size(size(img1));
    if s(2)==3 
        f1=rgb2gray(img1);
    else
        f1=img1;
    end 

    s1=size(size(fused));
    if s1(2)==3 
        f2=rgb2gray(fused);
    else
        f2=fused;
    end 

    G1=double(f1);
    G2=double(f2);
    [m1,n1]=size(G1);
    [m2,n2]=size(G2);
    m2=m1;
    n2=n1;
    X1=zeros(1,256);
    X2=zeros(1,256);
    result=0;

    for i=1:m1
        for j=1:n1
            X1(G1(i,j)+1)=X1(G1(i,j)+1)+1;
            X2(G2(i,j)+1)=X2(G2(i,j)+1)+1;
        end
    end

    for k=1:256
        P1(k)=X1(k)/(m1*n1);
        P2(k)=X2(k)/(m1*n1);
        if((P1(k)~=0)&(P2(k)~=0))
            result=P1(k)*log2(P1(k)/P2(k))+result;
        end
    end
    res0=result;
end