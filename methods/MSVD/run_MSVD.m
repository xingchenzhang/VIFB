% This is the code for MSVD image fusion algorithm.
% V. Naidu, ¡°Image fusion technique using multi-resolution singular value decomposition,¡± Defence Science Journal,
% vol. 61, no. 5, pp. 479¨C484, 2011.
% 
% The interface is created by the authors of VIFB.

function img = run_MSVD(imgA, imgB, visualization)

    addpath(genpath(cd));

    % IR image
    I = double(imread(imgB.img))/255;

    % VI image
    V = double(imread(imgA.img))/255;

    I = double(I);
    V = double(V);

    calc_metric = 0; % Calculate the metrices is time consuming, it is used for quantitative evaluation. Set it to 0 if you do not want to do it.

    tic;
    if size(V,3)==1   %for gray images       
        % %Pixel-and region-based image fusion with complex wavelets(2007)
        [M,N]=size(I);
        I4=imresize(I,[M+mod(M,2) N+mod(N,2)]);
        V4=imresize(V,[M+mod(M,2) N+mod(N,2)]);

        %Image Fusion technique using Multi-resolution singular Value decomposition(2011)
        %apply MSVD
        tic;
        [Y1, U1] = MSVD(I4);
        [Y2, U2] = MSVD(V4);

        %fusion starts
        X6.LL = 0.5*(Y1.LL+Y2.LL);

        D  = (abs(Y1.LH)-abs(Y2.LH)) >= 0; 
        X6.LH = D.*Y1.LH + (~D).*Y2.LH;
        D  = (abs(Y1.HL)-abs(Y2.HL)) >= 0; 
        X6.HL = D.*Y1.HL + (~D).*Y2.HL;
        D  = (abs(Y1.HH)-abs(Y2.HH)) >= 0; 
        X6.HH = D.*Y1.HH + (~D).*Y2.HH;

        %XX = [X.LL, X.LH; X.HL, X.HH];
        U = 0.5*(U1+U2);

        %apply IMSVD
        X6 = IMSVD(X6,U);
        toc;

        X6 = mat2gray(X6);
        imgf = X6;

    elseif  size(I,3)==1 
        imgf = zeros(size(V)); 
        M=size(I,1);
        N=size(I,2);
        for i=1:3
            I4=imresize(I,[M+mod(M,2) N+mod(N,2)]);
            V4=imresize(V(:,:,i),[M+mod(M,2) N+mod(N,2)]);

            [Y1, U1] = MSVD(I4);
            [Y2, U2] = MSVD(V4);

            %fusion starts
            X6{i}.LL = 0.5*(Y1.LL+Y2.LL);

            D  = (abs(Y1.LH)-abs(Y2.LH)) >= 0; 
            X6{i}.LH = D.*Y1.LH + (~D).*Y2.LH;
            D  = (abs(Y1.HL)-abs(Y2.HL)) >= 0; 
            X6{i}.HL = D.*Y1.HL + (~D).*Y2.HL;
            D  = (abs(Y1.HH)-abs(Y2.HH)) >= 0; 
            X6{i}.HH = D.*Y1.HH + (~D).*Y2.HH;

            %XX = [X.LL, X.LH; X.HL, X.HH];
            U = 0.5*(U1+U2);

            %apply IMSVD
            X6{i} = IMSVD(X6{i},U);

            imgf(:,:,i) = X6{i};
        end 

    else
        imgf = zeros(size(I)); 
        M=size(I,1);
        N=size(I,2);
        for i=1:3
            I4=imresize(I(:,:,i),[M+mod(M,2) N+mod(N,2)]);
            V4=imresize(V(:,:,i),[M+mod(M,2) N+mod(N,2)]);

            [Y1, U1] = MSVD(I4);
            [Y2, U2] = MSVD(V4);

            %fusion starts
            X6{i}.LL = 0.5*(Y1.LL+Y2.LL);

            D  = (abs(Y1.LH)-abs(Y2.LH)) >= 0; 
            X6{i}.LH = D.*Y1.LH + (~D).*Y2.LH;
            D  = (abs(Y1.HL)-abs(Y2.HL)) >= 0; 
            X6{i}.HL = D.*Y1.HL + (~D).*Y2.HL;
            D  = (abs(Y1.HH)-abs(Y2.HH)) >= 0; 
            X6{i}.HH = D.*Y1.HH + (~D).*Y2.HH;

            %XX = [X.LL, X.LH; X.HL, X.HH];
            U = 0.5*(U1+U2);

            %apply IMSVD
            X6{i} = IMSVD(X6{i},U);
            imgf(:,:,i) = X6{i};
        end
    end
    toc;

    img = im2uint8(imgf);
end
