% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
%
% Reference for the metric:
% G. Qu, D. Zhang, and P. Yan, “Information measure for performance of image fusion,” Electronics letters, vol. 38, no. 7,
% pp. 313–315, 2002.

function res = metricsMutinf(img1,img2,fused)

    fused = double(fused); 
    % Get the size of img 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);

    if b == 1
        g = Mutinf_single(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = Mutinf_single(img1(:,:,k), img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = Mutinf_single(img1(:,:,k), img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end

end


function output = Mutinf_single(img1,img2,fused)

  miVI = mutinf(img1, fused);
  miIR = mutinf(img2, fused);
  mi = miVI + miIR;
  output = mi;
  
end


function mi = mutinf(a, b)

    if size(a,3) > 1
        a=rgb2gray(a);
    end
    
    if size(b,3) > 1
        b=rgb2gray(b);
    end

    a=double(a);
    b=double(b);
    [Ma,Na] = size(a);
    [Mb,Nb] = size(b);
    M=min(Ma,Mb);
    N=min(Na,Nb);

    %初始化直方图数组
    hab = zeros(256,256);
    ha = zeros(1,256);
    hb = zeros(1,256);

    %归一化
    if max(max(a))~=min(min(a))
        a = (a-min(min(a)))/(max(max(a))-min(min(a)));
    else
        a = zeros(M,N);
    end

    if max(max(b))-min(min(b))
        b = (b-min(min(b)))/(max(max(b))-min(min(b)));
    else
        b = zeros(M,N);
    end

    a = double(int16(a*255))+1;
    b = double(int16(b*255))+1;

    %统计直方图
    for i=1:M
        for j=1:N
           indexx =  a(i,j);
           indexy = b(i,j) ;
           hab(indexx,indexy) = hab(indexx,indexy)+1;%联合直方图
           ha(indexx) = ha(indexx)+1;%a图直方图
           hb(indexy) = hb(indexy)+1;%b图直方图
       end
    end

    %计算联合信息熵
    hsum = sum(sum(hab));
    index = find(hab~=0);
    p = hab/hsum;
    Hab = sum(sum(-p(index).*log(p(index))));

    %计算a图信息熵
    hsum = sum(sum(ha));
    index = find(ha~=0);
    p = ha/hsum;
    Ha = sum(sum(-p(index).*log(p(index))));

    %计算b图信息熵
    hsum = sum(sum(hb));
    index = find(hb~=0);
    p = hb/hsum;
    Hb = sum(sum(-p(index).*log(p(index))));

    %计算a和b的互信息
    mi = Ha+Hb-Hab;

    %计算a和b的归一化互信息
    mi1 = hab/(Ha+Hb); 

end