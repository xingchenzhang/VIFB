% This is the program of FPDE: 
% D. P. Bavirisetti, G. Xiao, and G. Liu, ¡°Multi-sensor image
% fusion based on fourth order partial differential equations,¡± in
% 2017 20th International Conference on Information Fusion
% (Fusion). IEEE, 2017, pp. 1¨C9.
% Codes are provided by the authors of FPDE.
%
% The interface is created by the authors of VIFB. Necessary modifications
% are made to be integrated into VIFB. 

function img = run_FPDE(imgVI, imgIR, visualization)

    % IR image
    I1 = imread(imgIR.img);

    % VI image
    I2 = imread(imgVI.img);

    if visualization == 1
        figure, imshow((uint8(I1)));
    end

    if visualization == 1
        figure, imshow(uint8(I2));
    end

    tic;
    if size(I2, 3) == 1
        fuseimage = FPDE(I1, I2);
    elseif size(I1,3) == 1
        fuseimage = zeros(size(I2));
        for i=1:3
            fuseimage(:,:,i) = FPDE(I1,I2(:,:,i));    
        end       
    else
        fuseimage = zeros(size(I2));
        for i=1:3
           fuseimage(:,:,i) = FPDE(I1(:,:,i),I2(:,:,i));    
        end    
    end   
    toc;
    if visualization == 1
        figure,imshow(fuseimage, []);
    end 
    img = uint8(fuseimage);
end
    
function res = FPDE (I1,I2)

    % Assigning values to the parameters
    n=15; 
    dt=0.9;
    k=4;

    % Decomposing input images as base and detail layers

    [A1]=fpdepyou(I1,n);
    [A2]=fpdepyou(I2,n);
    D1=double(I1)-double(A1);
    D2=double(I2)-double(A2);

    A(:,:,1)=A1;
    A(:,:,2)=A2;

    D(:,:,1)=D1;
    D(:,:,2)=D2;

    % Detail layer fusion 

    C1 = cov([D1(:) D2(:)]);
    [V11, D11] = eig(C1);
    if D11(1,1) >= D11(2,2)
        pca1 = V11(:,1)./sum(V11(:,1));
    else  
        pca1 = V11(:,2)./sum(V11(:,2));
    end
    imf1 = pca1(1)*D1 + pca1(2)*D2;

    % Base layer fusion 
    imf2=(0.5*A1+0.5*A2);

    % Final fused image
    res=(double(imf1)+double(imf2));
end
