% The source code is from https://github.com/zhengliu6699/imageFusionMetrics/blob/master/metricChen.m
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
% 
% Reference for the metric
% H. Chen and P. K. Varshney, “A human perception inspired quality metric for image fusion based on regional 
% information,” Information fusion, vol. 8, no. 2, pp. 193–207, 2007.

function res = metricsQcv(img1,img2,fused)

    fused = double(fused); 
    % Get the size of img 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);

    if b == 1
        g = Qcv(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = Qcv(img1(:,:,k),img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = Qcv(img1(:,:,k),img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end
end



function output = Qcv(im1, im2, fused)

    % function res=metricChen(im1,im2,fused)
    %
    % This function implements Chen's algorithm for fusion metric.
    % im1, im2 -- input images;
    % fused      -- fused image;
    %
    % IMPORTANT: The size of the images need to be 2X.
    % See also: evalu_fusion.m
    %
    % Z. Liu [July 2009]
    %

    % Ref: A human perception inspired qualtiy metric for image fusion based on
    % regional information, Information Fusion, Vol.8, pp193-207, 2007
    % by Hao Chen et al.
    % 

    %% set up the constant values

    alpha_c=1;
    alpha_s=0.685;
    f_c=97.3227;
    f_s=12.1653;

    % local window size 16 x 16
    windowSize=16;

    % alpha = 1, 2, 3, 4, 5, 10, 15. This value is adjustable.;-)
    alpha=5;


    %% pre-processing
    % im1=double(im1);
    % im2=double(im2);
    % fused=double(fused);

    % im1=normalize1(im1);
    % im2=normalize1(im2);
    % fused=normalize1(fused);
% 
%         s=size(size(im1));
%         if s(2)==3 %判断是灰度图还是RGB
%             im1=rgb2gray(im1);
%         else
%             im1=im1;
%         end 
% 
%         s1=size(size(im2));
%         if s1(2)==3 %判断是灰度图还是RGB
%             im2=rgb2gray(im2);
%         else
%             im2=im2;
%         end 
% 
%         s2=size(size(fused));
%         if s2(2)==3 %判断是灰度图还是RGB
%             fused=rgb2gray(fused);
%         else
%             fused=fused;
%         end 

        im1 = double(im1);
        im2 = double(im2);
        fused = double(fused);

        im1=normalize1(im1);
        im2=normalize1(im2);
        fused=normalize1(fused);

    %% Step 1: extract edge information

    flt1=[-1 0 1 ; -2 0 2 ; -1 0 1];
    flt2=[-1 -2 -1; 0 0 0; 1 2 1];

    % 1) get the map

    fuseX=filter2(flt1,fused,'same');
    fuseY=filter2(flt2,fused,'same');
    fuseG=sqrt(fuseX.*fuseX+fuseY.*fuseY);

    buffer=(fuseX==0);
    buffer=buffer*0.00001;
    fuseX=fuseX+buffer;
    fuseA=atan(fuseY./fuseX);


    img1X=filter2(flt1,im1,'same');
    img1Y=filter2(flt2,im1,'same');
    im1G=sqrt(img1X.*img1X+img1Y.*img1Y);

    buffer=(img1X==0);
    buffer=buffer*0.00001;
    img1X=img1X+buffer;
    im1A=atan(img1Y./img1X);

    img2X=filter2(flt1,im2,'same');
    img2Y=filter2(flt2,im2,'same');
    im2G=sqrt(img2X.*img2X+img2Y.*img2Y);

    buffer=(img2X==0);
    buffer=buffer*0.00001;
    img2X=img2X+buffer;
    im2A=atan(img2Y./img2X);

    %% calculate the local region saliency

    % seperate the image into local regions
    [hang,lie]=size(im1);
    H=floor(hang/windowSize);
    L=floor(lie/windowSize);

    fun=@(x) sum(sum(x.^alpha)); 

    ramda1=blkproc(im1G, [windowSize windowSize],[0 0],fun);
    ramda2=blkproc(im2G, [windowSize windowSize],[0 0],fun);

    %% similarity measurement

    f1=im1-fused;
    f2=im2-fused;

    % filtering with CSF filters
    % first build the three filters in frequency domain:

    [u,v]=freqspace([hang,lie],'meshgrid');

    %---------------
    % length 1
    %u=lie/2*u; v=hang/2*v;

    % length 2
    %u=lie/4*u; v=hang/4*v;

    % length 3
    u=lie/8*u; v=hang/8*v;


    r=sqrt(u.^2+v.^2);

    % Mannos-Skarison's filter
    theta_m=2.6*(0.0192+0.144*r).*exp(-(0.144*r).^1.1);

    % Daly's filter
    % avoid being divided by zero
    index=find(r==0);
    r(index)=1;

    buff=0.008./(r.^3)+1;
    %buff(index)=0;
    buff=buff.^(-0.2);
    buff1=-0.3*r.*sqrt(1+0.06*exp(0.3*r));

    theta_d=((buff).^(-0.2)).*(1.42*r.*exp(buff1));
    theta_d(index)=0;
    clear buff;
    clear buff1;

    % Ahumada filter
    theta_a=alpha_c*exp(-(r/f_c).^2)-alpha_s*exp(-(r/f_s).^2);

    % second, filter the image in frequency domain
    ff1=fft2(f1); 
    ff2=fft2(f2);

    % here filter 1 is used.
    Df1=ifft2(ifftshift(fftshift(ff1).*theta_m));
    Df2=ifft2(ifftshift(fftshift(ff2).*theta_m));

    fun2=@(x) mean2(x.^2);
    D1=blkproc(Df1, [windowSize windowSize],[0 0],fun2);
    D2=blkproc(Df2, [windowSize windowSize],[0 0],fun2);

    %% global quality

    % pay attention to this equation, the author might be WRONG!!!
    %Q=sum(sum(ramda1.*D1+ramda2.*D2/(ramda1+ramda2)));

    Q=sum(sum(ramda1.*D1+ramda2.*D2))/sum(sum(ramda1+ramda2));

    output=Q;

end

function RES=normalize1(data)

    % function RES=normalize1(data)
    %
    % This function is to NORMALIZE the data. 
    % The data will be in the interval 0-255 (gray level) and pixel value has
    % been rounded to an integer.
    % 
    % See also: normalize.m 
    %
    % Z. Liu @NRCC (Aug 24, 2009)

    data=double(data);
    da=max(data(:));
    xiao=min(data(:));
    if (da==0 & xiao==0)
        RES=data;
    else
        newdata=(data-xiao)/(da-xiao);
        RES=round(newdata*255);
    end
end


