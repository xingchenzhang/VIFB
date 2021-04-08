% This is the program of TIF:
% D. P. Bavirisetti and R. Dhuli, ¡°Two-scale image fusion of visible and infrared images using saliency detection,¡±
% Infrared Physics & Technology, vol. 76, pp. 52¨C64, 2016 
%
% Codes are provided by the authors of TIF. 
% The interface is created by the authors of VIFB.

function img = run_TIF(imgVI, imgIR, visualization)

    % IR image
    img1 = imread(imgIR.img);

    % VI image
    img2 = imread(imgVI.img);
    
    tic;   
    if size(img2, 3) == 1
        fuseimage = TIF(img2, img1);
    elseif size(img1,3) == 1
        fuseimage = zeros(size(img2));
        for i=1:3
            fuseimage(:,:,i) = TIF(img2(:,:,i),img1);    
        end       
    else
        fuseimage = zeros(size(img2));
        for i=1:3
           fuseimage(:,:,i) = TIF(img2(:,:,i),img1(:,:,i));    
        end    
    end       
    toc;
    if visualization == 1
        figure,imshow(fuseimage, []);
    end 
    img = uint8(fuseimage);   
end
    
function res = TIF(I1, I2)
    
    N=3;
    Med1= medfilt2(I1, [N N]);
    M=35;
    h=1/(M*M)*ones(M);
    b1=imfilter(double(I1),double(h),'circular');
    d1=double(I1)-b1; 
    S1=(b1-double(Med1)).^2;
    Med2= medfilt2(I2, [N N]);
    b2=imfilter(double(I2),double(h),'circular');
    d2=double(I2)-b2;
    S2=(b2-double(Med2)).^2;
    w1=S1./(S1+S2);
    w2=S2./(S1+S2);
    F1=double(w1).*double(d1)+double(w2).*double(d2);
    F2=0.5*b1+0.5*b2;
    FF=double(F1)+F2;

    res = FF;   
end
