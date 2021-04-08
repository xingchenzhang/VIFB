% The source code is from online.
% The interface is modified by the authors of VIFB to integrate it into VIFB.
%
% Reference for the metric:
% C. S. Xydeas and P. V. V., ¡°Objective image fusion performance measure,¡± Military Technical Courier, vol. 36, no. 4,
% pp. 308¨C309, 2000

function res = metricsQabf(img1, img2, fused)

    fused = double(fused); 
    % Get the size of img 
    [m,n,b] = size(fused); 
    [m1,n1,b1] = size(img2);

    if b == 1
        g = Qabf(img1,img2,fused);
        res = g;
    elseif b1 == 1
        for k = 1 : b 
           g(k) = Qabf(img1(:,:,k), img2,fused(:,:,k)); 
        end 
        res = mean(g); 
    else
        for k = 1 : b 
            g(k) = Qabf(img1(:,:,k), img2(:,:,k),fused(:,:,k)); 
        end 
        res = mean(g); 
    end

end

function Qabf1= Qabf(strA, strB, strF)

    % model parameters
        L=1; Tg=0.9994;kg=-15;Dg=0.5;Ta=0.9879;ka=-22;Da=0.8;    

    % Sobel Operator
    h1=[1 2 1;0 0 0;-1 -2 -1]; h2=[0 1 2;-1 0 1;-2 -1 0]; h3=[-1 0 1;-2 0 2;-1 0 1];
    % if y is the response to h1 and x is the response to h3;
    % then the intensity is sqrt(x^2+y^2) and orientation is arctan(y/x);
     pA = double(strA);
     pB = double(strB); 

     if size(pA,3) > 1
         pA=rgb2gray(pA);
     end

     if size(pB,3) > 1
        pB=rgb2gray(pB);
     end

     if size(strF,3) > 1
        strF=rgb2gray(strF);
     end

     pA = 255*im2double(pA);
     pB = 255*im2double(pB);
     pF = 255*im2double(strF); 

    SAx = conv2(pA,h3,'same'); SAy = conv2(pA,h1,'same');
    gA = sqrt(SAx.^2 + SAy.^2); 
    [M,N] = size(SAx); aA = zeros(M,N);
    for i=1:M
        for j=1:N
            if ( SAx(i,j) == 0 ) aA(i,j) = pi/2;
            else
                aA(i,j) = atan(SAy(i,j)/SAx(i,j));
            end
        end
    end

    SBx = conv2(pB,h3,'same'); SBy = conv2(pB,h1,'same');
    gB = sqrt(SBx.^2 + SBy.^2); 
    [M,N] = size(SBx); aB = zeros(M,N);
    for i=1:M
        for j=1:N
            if ( SBx(i,j) == 0 ) aB(i,j) = pi/2;
            else
                aB(i,j) = atan(SBy(i,j)/SBx(i,j));
            end
        end
    end

    SFx = conv2(pF,h3,'same'); SFy = conv2(pF,h1,'same');
    gF = sqrt(SFx.^2 + SFy.^2); 
    [M,N] = size(SAx); aF = zeros(M,N);
    for i=1:M
        for j=1:N
            if ( SFx(i,j) == 0 ) aF(i,j) = pi/2;
            else
                aF(i,j) = atan(SFy(i,j)/SFx(i,j));
            end
        end
    end

    % the relative strength and orientation value of GAF,GBF and AAF,ABF;
    GAF = zeros(M,N); AAF = zeros(M,N); QgAF = zeros(M,N); QaAF = zeros(M,N); QAF = zeros(M,N);
    for i=1:M
        for j=1:N

            if ( gA(i,j) > gF(i,j))  GAF(i,j) = gF(i,j)/gA(i,j);
            else
                if ( gA(i,j) == gF(i,j) )  GAF(i,j) = gF(i,j);
                else
                    GAF(i,j) = gA(i,j) / gF(i,j);
                end
            end 
            AAF(i,j) = 1 - abs(aA(i,j)-aF(i,j))/(pi/2);

            QgAF(i,j) = Tg / (1 + exp(kg*( GAF(i,j) - Dg )));
            QaAF(i,j) = Ta / (1 + exp(ka*( AAF(i,j) - Da )));

            QAF(i,j) = QgAF(i,j) * QaAF(i,j);
        end
    end

    GBF = zeros(M,N); ABF = zeros(M,N); QgBF = zeros(M,N); QaBF = zeros(M,N); QBF = zeros(M,N);
    for i=1:M
        for j=1:N

            if ( gB(i,j) > gF(i,j))  GBF(i,j) = gF(i,j)/gB(i,j);
            else
                if ( gB(i,j) == gF(i,j) )  GBF(i,j) = gF(i,j);
                else
                    GBF(i,j) = gB(i,j) / gF(i,j);
                end
            end 
            ABF(i,j) = 1 - abs(aB(i,j)-aF(i,j))/(pi/2);

            QgBF(i,j) = Tg / (1 + exp(kg*( GBF(i,j) - Dg )));
            QaBF(i,j) = Ta / (1 + exp(ka*( ABF(i,j) - Da )));

            QBF(i,j) = QgBF(i,j) * QaBF(i,j);
        end
    end

    deno = sum(sum( gA + gB ));
    nume = sum(sum( QAF.*gA + QBF.*gB ));
    Qabf1 = nume / deno;
    Qabf=Qabf1;
    Qabf;

end