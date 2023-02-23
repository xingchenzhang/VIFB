disp(' ');
disp('fdct_wrapping_demo_denoise_enhanced.m - Image denoising using Curvelets');
disp(' ');
disp('This is an extension of the file fdct_wrapping_demo_denoise.m')
disp('We try different tricks to obtain better image denoising: (1) curvelets');
disp('at the finest scale, (2) complex block thresholding, (3) cycle spinning,');
disp('(4) ell-1 iterations. The user can fine-tune some parameters, see the');
disp('beginning of the code.')
disp(' ');
disp('By Laurent Demanet, 2005');
disp(' ');

% fdct_wrapping_demo_denoise_enhanced.m - Image denoising using Curvelets

N = 512;
sigma = 0.1*255;            % noise stdev is 10%

% Tuning parameters
finest = 1;                 % 1: curvelets at finest scale, 2: wavelets at finest scale
nbscales = log2(N) - 3;
nbangles_coarse = 8;
nsigmas_coarse = 2.2;       % threshold proportional to nsigmas_coarse*sigma at all scales except finest
nsigmas_fine = 2.5;         % threshold proportional to nsigmas_fine*sigma at finest scale
nshifts = 1;                % number of translations (per dimension) considered in the cycle-spinning
nell1 = 0;                  % number of ell-1 iterations
neighb_weight = 0.5;        % for group thresholding, weight assigned to neighboring curvelets
tuning_neighb = 0.6;        % for group thresholding, parameter used to renormalize the weighted sum of coefficients squared
% Don't forget that you can choose between diagonal real, diagonal complex
% and block complex thersholding by uncommenting portions of the code

img = double(imread('barbara512.jpg'));
noise = randn(N) * sigma;
noisy_img = img + noise;

% get L2norm, put into E 
F = ones(N,N);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
disp('Computing L^2 norms ...');
tic;
C = fdct_wrapping(X,0,finest,nbscales,nbangles_coarse);

E = cell(size(C));
for s=1:length(C)
  E{s} = cell(size(C{s}));
  for w=1:length(C{s})
    A = C{s}{w};
    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
  end
end

disp('Computing parameters ...');
[X_rows, X_cols, F_rows, F_cols, N_rows, N_cols] = fdct_wrapping_param(C,N,N);

% Cycle spinning
n = 0;
restored_img = 0;
for xshift = 1:nshifts
    for yshift = 1:nshifts

        shift_img = circshift(noisy_img,[xshift yshift]);
        n = n + 1;

        disp(['Direct transform, shift nr. ',num2str(n),' ...']);
        C = fdct_wrapping(shift_img,0,finest,nbscales,nbangles_coarse);

        % Thresholding
        disp('Thresholding ...')
        thresh = nsigmas_coarse * sigma;
        for j = 1:length(C)
            if j == length(C), thresh = nsigmas_fine * sigma; end;
            for l = 1:length(C{j})
                
                thresh_jl = thresh*E{j}{l};
                
% Uncomment for 'diagonal real thresholding'
%                recjl = real(C{j}{l});     
%                imcjl = imag(C{j}{l});
%                recjl = recjl .* (abs(recjl) > thresh_jl);
%                imcjl = imcjl .* (abs(imcjl) > thresh_jl);
%                C{j}{l} = recjl + sqrt(-1)*imcjl;                

% Uncomment for 'diagonal complex thresholding'
%                modcjl = abs(C{j}{l});     
%                argcjl = C{j}{l} ./ modcjl;
%                modcjl = modcjl .* (modcjl > thresh_jl);
%                C{j}{l} = argcjl .* modcjl;

% Uncomment for 'block complex thresholding'
                modcjl = abs(C{j}{l});          
                argcjl = C{j}{l} ./ modcjl;
                rowstep = N/N_rows{j}{l};
                colstep = N/N_cols{j}{l};
                evenquad = ~mod(ceil(l*4/length(C{j})),2);
                if evenquad,
                    if (j == 1)|(finest==2 & j==nbscales), fcolsjl = 1; else fcolsjl = F_cols{j}{l}; end;
                    rowshift = - round(F_rows{j}{l}/fcolsjl * rowstep);
                    testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[1 0]).^2 + circshift(modcjl,[-1 0]).^2 + ...
                        circshift(modcjl,[rowshift 1]).^2 + circshift(modcjl, [-rowshift -1]).^2));
                else
                    if (j == 1)|(finest==2 & j==nbscales), frowsjl = 1; else frowsjl = F_rows{j}{l}; end;
                    colshift = - round(F_cols{j}{l}/frowsjl * colstep);
                    testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[0 1]).^2 + circshift(modcjl,[0 -1]).^2 + ...
                        circshift(modcjl,[1 colshift]).^2 + circshift(modcjl, [-1 -colshift]).^2));
                end;
                testcjl = testcjl ./ sqrt(1+4*neighb_weight*tuning_neighb);
                modcjl = modcjl .* (testcjl > thresh_jl);
                C{j}{l} = argcjl .* modcjl;

            end
        end

        disp('Inverse transform ...');
        temp_restored = real(ifdct_wrapping(C,0,N,N));
        
        % L1 iterations

        for nupdate = 1:nell1
            nupdate
            disp(['Direct transform, within ell-1 iteration, shift nr.',num2str(n),' ...']);
            D = fdct_wrapping(temp_restored,0,finest,nbscales,nbangles_coarse);
            disp('Thresholding ...');
            thresh = nsigmas_coarse * sigma;
            for j = 1:length(C)
                if j == length(C), thresh = nsigmas_fine * sigma; end;
                for l = 1:length(C{j})
                    thresh_jl = thresh*E{j}{l};
                    
% Uncomment for 'diagonal real thresholding'
%                    redjl = real(D{j}{l});   
%                    imdjl = imag(D{j}{l});
%                    recjl = real(C{j}{l});
%                    imcjl = imag(C{j}{l});
%                    redjl = (recjl - redjl) .* (abs(recjl) > thresh_jl);
%                    imdjl = (imcjl - imdjl) .* (abs(imcjl) > thresh_jl);
%                    D{j}{l} = redjl + sqrt(-1)*imdjl;


% Uncomment for 'diagonal complex thresholding'
%                    modcjl = abs(C{j}{l});       
%                    moddjl = abs(D{j}{l});          
%                    argcjl = C{j}{l} ./ (modcjl+1e-16);
%                    argdjl = D{j}{l} ./ (moddjl+1e-16);
%                    moddjl = (modcjl - moddjl) .* (modcjl > thresh_jl);
%                    D{j}{l} = argdjl .* moddjl;

% Uncomment for 'block complex thresholding'
                    modcjl = abs(C{j}{l});          
                    moddjl = abs(D{j}{l});
                    argcjl = C{j}{l} ./ (modcjl+1e-16);
                    argdjl = D{j}{l} ./ (moddjl+1e-16);
                    rowstep = N/N_rows{j}{l};
                    colstep = N/N_cols{j}{l};
                    evenquad = ~mod(ceil(l*4/length(C{j})),2);
                    if evenquad,
                        if (j == 1)|(finest==2 & j==nbscales), fcolsjl = 1; else fcolsjl = F_cols{j}{l}; end;
                        rowshift = - round(F_rows{j}{l}/fcolsjl * rowstep);
                        testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[1 0]).^2 + circshift(modcjl,[-1 0]).^2 + ...
                            circshift(modcjl,[rowshift 1]).^2 + circshift(modcjl, [-rowshift -1]).^2));
                    else
                        if (j == 1)|(finest==2 & j==nbscales), frowsjl = 1; else frowsjl = F_rows{j}{l}; end;
                        colshift = - round(F_cols{j}{l}/frowsjl * colstep);
                        testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[0 1]).^2 + circshift(modcjl,[0 -1]).^2 + ...
                            circshift(modcjl,[1 colshift]).^2 + circshift(modcjl, [-1 -colshift]).^2));
                    end;
                    testcjl = testcjl ./ sqrt(1+4*neighb_weight*tuning_neighb);
                    moddjl = (modcjl - moddjl) .* (testcjl > thresh_jl);
                    D{j}{l} = argdjl .* moddjl;
                    
                end
            end
            disp('Inverse transform ...')
            temp_update = real(ifdct_wrapping(D,0,N,N));
            max(max(abs(temp_update)))

            temp_restored = temp_restored + temp_update;
        end;
        
        temp_restored = circshift(temp_restored,[-xshift, -yshift]);
        restored_img = (n-1)/n*restored_img + 1/n*temp_restored;
        
    end
end

MSE = sum(sum((img-restored_img).^2))/N^2;
PSNR = 20*log10(255/sqrt(MSE));
disp(['PSNR = ',num2str(PSNR)]);
disp(['Time elapsed = ',num2str(toc)]);

figure(1); clf; imagesc(img); colormap gray; axis('image');
title('Original Image')
figure(2); clf; imagesc(noisy_img); colormap gray; axis('image');
title('Noisy Image')
figure(3); clf; imagesc(restored_img); colormap gray; axis('image');
title('Restored image');
