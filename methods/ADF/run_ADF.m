% This is the program of ADF from the paper:
%
% D. P. Bavirisetti and R. Dhuli, ¡°Fusion of infrared and visible
% sensor images based on anisotropic diffusion and karhunenloeve transform,¡± IEEE Sensors Journal, vol. 16, no. 1, pp.
% 203¨C209, 2016.
%
% The function ADF() is provided by the authors of ADF.
% The interface is created by the authors of VIFB.

function img = run_ADF(imgVI, imgIR, visualization)
     % IR image
     I1 = imread(imgIR.img);

     % VI image
     I2 = imread(imgVI.img);

     if visualization == 1
        figure, imshow((uint8(I1)));
        figure, imshow(uint8(I2));
     end
     
     tic;
     if size(I2, 3) == 1
         fuseimage = ADF(I1, I2);
     elseif size(I1,3) == 1
         fuseimage = zeros(size(I2));
         for i=1:3
            fuseimage(:,:,i) = ADF(I1,I2(:,:,i));    
         end       
     else
         fuseimage = zeros(size(I2));
         for i=1:3
            fuseimage(:,:,i) = ADF(I1(:,:,i),I2(:,:,i));    
         end    
     end
     toc;   

     if visualization == 1
         figure, imshow((fuseimage), [])
     end

     img = uint8(fuseimage);
end

function res = ADF(I1, I2)
     %ANISOTROPIC DIFFUSION
     num_iter = 10;
     delta_t = 0.15;
     kappa = 30;
     option = 1;

     A1 = anisodiff2D(I1,num_iter,delta_t,kappa,option);
     A2= anisodiff2D(I2,num_iter,delta_t,kappa,option);

     D1=double(I1)-A1;
     D2=double(I2)-A2;

     C1 = cov([D1(:) D2(:)]);
     [V11, D11] = eig(C1);
     if D11(1,1) >= D11(2,2)
        pca1 = V11(:,1)./sum(V11(:,1));
     else  
        pca1 = V11(:,2)./sum(V11(:,2));
     end

     imf1 = pca1(1)*D1 + pca1(2)*D2;
     imf2=(0.5*A1+0.5*A2);

     res=(double(imf1)+double(imf2));
end
 