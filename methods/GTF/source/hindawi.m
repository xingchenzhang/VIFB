function[] = hindawi()
%  
% Simulation code for [1], 
%
% [1] Paul Rodriguez "Review: Total Variation Regularization Algorithms for 
%  		      Images Corrupted With Different Noise Models"
%  
% 
%
% Legal:
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%  




l2TV         = 1;
l2TV_nqp     = 2;
l1TV         = 3;
l1TV_adapt   = 4;
poiTV_denoi  = 5;
poiTV_deconv = 6;
gammaTV      = 7;
mixTV        = 8;
nTV          = mixTV;


% Test Images (for each case)

l2Imgs = {'goldhl', 'Clena'};           % To be corrupted with Gaussian noise
l1Imgs = {'Cbarb', 'boats'};            % To be corrupted with Salt & Pepper noise
poiImgs = {'cman512', 'cpeppers'};     % To be corrupted with Poisson noise
gammaImgs = {'tank', 'Clena'};          % To be corrupted with Gamma noise
mixImgs = {'cman512'};                  % To be corrupted with mixed Gaussian + Impulse noise


Tests = cell(nTV-1, 1);

Tests{l2TV} = l2Imgs;
Tests{l2TV_nqp} = l2Imgs;
Tests{l1TV} = l1Imgs;
Tests{l1TV_adapt} = l1Imgs;
Tests{poiTV_denoi} = poiImgs;
Tests{poiTV_deconv} = poiImgs;

%  Tests{gammaTV} = gammaImgs;
%  Tests{mixTV} = mixImgs;


% kernel

  kernel = fspecial('disk',3.2);

  K = @(x) imfilter(x, kernel, 'symmetric','conv');

  KT = @(x) K(x);
  KC = {K, KT};


% Noise level (define the noise level to corrupt images)

l2Noise = [0.1 0.2];
l1Noise = [0.3 0.5 0.8];
randomNoise = [0.1 0.2 0.3];
poiNoise = [5 30 100 255];

Noise = cell(nTV-1, 1);

Noise{l2TV} = l2Noise;
Noise{l2TV_nqp} = l2Noise;
Noise{l1TV} = l1Noise;
Noise{l1TV_adapt} = l1Noise;
Noise{poiTV_denoi} = poiNoise;
Noise{poiTV_deconv} = poiNoise;



% TV Noise Model (to corrupt imput images)

NoiseModel = cell(nTV-1, 1);

NoiseModel{l2TV} = @(Img, sigma) imnoise(K(Img),'gaussian', 0, (sigma^2) );
NoiseModel{l2TV_nqp} = @(Img, sigma) imnoise(K(Img),'gaussian', 0, (sigma^2) );
NoiseModel{l1TV} = @(Img, spnoise) imnoise(Img, 'salt & pepper', spnoise);
NoiseModel{l1TV_adapt} = @(Img, spnoise) imnoise(Img, 'salt & pepper', spnoise);

NoiseModel{poiTV_denoi} = @(Img, M) poissrnd( M*NormalizeMax(Img) );
NoiseModel{poiTV_deconv} = @(Img, M) poissrnd( M*NormalizeMax(K(Img)) );

% ----------------------
% Define TV parameters
% ----------------------

% Loops

TVLoops = zeros(nTV-1, 1);

TVLoops(l2TV)           = 5;
TVLoops(l2TV_nqp)       = 5;
TVLoops(l1TV)           = 4;
TVLoops(l1TV_adapt)     = 4;

TVLoops(poiTV_denoi)    = 6;
TVLoops(poiTV_deconv)   = 6;

% regularization parameters

TVLambda = zeros(nTV-1, 5);

TVLambda(l2TV,1) = 0.05; TVLambda(l2TV,2) = 0.1; 
TVLambda(l2TV_nqp,1) = 0.05; TVLambda(l2TV_nqp,2) = 0.1; 

TVLambda(l1TV,1) = 1.2; TVLambda(l1TV,2) = 1.4; TVLambda(l1TV,3) = 1.5;
TVLambda(l1TV_adapt,1) = 1.0; TVLambda(l1TV_adapt,2) = 1.0; TVLambda(l1TV_adapt,3) = 1.0;

TVLambda(poiTV_denoi,1:4) = 2./[0.35, 0.15, 0.075, 0.025 ];




%  Wrappers to TV functions)

RestoreTV = cell(nTV-1, 1);

RestoreTV{l2TV} = @(b, lambda, loops, KC) runl2TV(b, lambda, loops, KC);
RestoreTV{l2TV_nqp} = @(b, lambda, loops, KC) runl2TV_nqp(b, lambda, loops, KC);
RestoreTV{l1TV} = @(b, lambda, loops, dummy) runl1TV(b, lambda, loops, []);
RestoreTV{l1TV_adapt} = @(b, lambda, loops, dummy) runl1TV_Adapt(b, lambda, loops);

RestoreTV{poiTV_denoi} = @(b, lambda, loops, dummy) run_poiTV(b, lambda, loops, {});
RestoreTV{poiTV_deconv} = @(b, lambda, loops, KC) run_poiTV(b, lambda, loops, KC);

% additional vars

t1 = zeros(nTV-1, 2, 2);
t2 = zeros(nTV-1, 2, 2);

nmpdef;


%%%%%%%%%%%%%%%%%%%%%%%%
%        Setup         %
%%%%%%%%%%%%%%%%%%%%%%%%

SSIM_CODE = exist('ssim_index');
if( SSIM_CODE == 0 )
  disp('NOTE:');
  disp('  The function ssim_index (code for [SSIM]) is not in your path...');
  disp('  You may download it from:');
  disp('  http://www.ece.uwaterloo.ca/~z70wang/research/ssim/');
  disp('  ');
  disp('  Disabling SSIM reports');
  disp('  ');
  disp('  [SSIM]  Z. Wang, A. Bovik, H. Sheikh and E. Simoncelli, ');
  disp('         "Image quality assessment: From error visibility to');
  disp('         structural similarity " ');
  disp('         IEEE Transactions on Image Processing, 2004, 13:4(600-612).');

end



t=5;

str_all = sprintf(' \n');



while( iscell(Tests{t}) )

  str = sprintf(' \n'); str_all = [ str_all str ];
  str = sprintf(' Img \t\t Noise \t\t SNR (db) \t\t Time (s) \t\t SSIM \n');
  str_all = [ str_all str ];


  NImgs = length( Tests{t} );

  for k = 1: NImgs,

  switch lower( Tests{t}{k} )

    case{'lena'}
      I = double( imread('gray_imgs/lena_gray_512.png') ) / 255;
      SSIM_COLOR = 1;

    case{'peppers'}
      I = double( imread('gray_imgs/peppers_gray.png') ) / 255;
      SSIM_COLOR = 1;

    case{'boats'}
      I = double( imread('gray_imgs/boats_gray.png') ) / 255;
      SSIM_COLOR = 1;

    case{'goldhl'}
      I = double( imread('gray_imgs/goldhill_gray.png') ) / 255;
      SSIM_COLOR = 1;

    case{'cman'}
      I = double( imread('gray_imgs/cameraman_256x256.tiff') ) / 255;
      SSIM_COLOR = 1;

    case{'cman512'}
      I = double( imread('gray_imgs/cameraman_512x512.tiff') ) / 255;
      SSIM_COLOR = 1;

    case{'clena'}
      I = double( imread('color_imgs/lena_color_512.png') ) / 255;
      SSIM_COLOR = 0;

    case{'cbarb'}
      I = double( imread('color_imgs/barbara_color.png') ) / 255;
      SSIM_COLOR = 0;

    case{'cpeppers'}
      I = double( imread('color_imgs/peppers_color.png') ) / 255;
      SSIM_COLOR = 0;

    case{'cmandrill'}
      I = double( imread('color_imgs/mandrill_color.png') ) / 255;
      SSIM_COLOR = 0;

  end % --- END (switch) ---


  for n = 1:length(Noise{t})

    % add noise
    b = NoiseModel{t}( I, Noise{t}(n) );
    figure; imagesc( Normalize(b) );  colormap gray; axis image; axis off;


    % Restore via TV
    [u t1(t,k,n) t2(t,k,n)] = RestoreTV{t}(b, TVLambda(t,n), TVLoops(t), KC );

    % Restore performance
    [snrRec(t,k,n), ssimRec(t,k,n)] = computePerf(I, u, SSIM_CODE*SSIM_COLOR, 0);

    if(n==1)
      str = sprintf('%s \t\t %1.2f \t\t %2.2f \t\t %2.2f (%2.2f) \t\t %1.3f \n', ...
                    lower( Tests{t}{k} ), Noise{t}(n), snrRec(t,k,n), t1(t,k,n), t2(t,k,n), ssimRec(t,k,n) );
    else
      str = sprintf('     \t\t %1.2f \t\t %2.2f \t\t %2.2f (%2.2f) \t\t %1.3f \n', ...
                     Noise{t}(n), snrRec(t,k,n), t1(t,k,n), t2(t,k,n), ssimRec(t,k,n) );
    end
    str_all = [ str_all str ];


  end % --- END (FOR(n)) ---

  end % --- END (FOR(k)) ---

  disp(str_all);
  

 t = t + 1;

 if(t==6) break; end

end


end

% ============================================

%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

function[u] = Normalize(x)

u = (x - min(x(:)))/(max(x(:)) - min(x(:)));


end

% ============================================

function[u] = NormalizeMax(x)

u = x/max(abs(x(:)));


end

% ============================================

function[u, timePerf, t2] = runl2TV(b, lambda, loops, KC)


    pars = irntvInputPars('l2tv');

    pars.pcgtol_ini = 1e-4;
    pars.loops      = loops;
%      pars.U0         = b;

    pars.adapt_epsR   = 1;
    pars.epsR_cutoff  = 0.01;
    pars.adapt_epsF   = 1;
    pars.epsF_cutoff  = 0.05; 

    tic;
    u = irntv(b, KC, lambda, pars);
    timePerf = toc;

    t2 = 0;

end

% ============================================

function[u, timePerf, t2] = runl2TV_nqp(b, lambda, loops, KC)

    pars = irntvInputPars('l2tv_nqp');

    pars.pcgtol_ini = 1e-4;
    pars.loops      = loops;
%      pars.U0         = b;

    pars.adapt_epsR   = 1;
    pars.epsR_cutoff  = 0.01;
    pars.adapt_epsF   = 1;
    pars.epsF_cutoff  = 0.1; 

    pars.loops_NQP    = 15;
    pars.gamma_NQP    = 0.5e-3;
    pars.alpha_NQP    = 5e-1;
    pars.vmax_NQP     = 1+1;            % maximum value of output image


    tic;
    u = irntv(b+1, KC, lambda, pars)-1;
    timePerf = toc;

    t2 = 0;

end

% ============================================

function[u, timePerf, t2] = runl1TV(b, lambda, loops, KC)


    if(nargin < 6)
      flagTitle = 0;
    end

    pars = irntvInputPars('l1tv');

    pars.pcgtol_ini = 1e-4;
    pars.epsF       = 1e-2;    
    pars.epsR       = 1e-4;
    pars.loops      = loops;


    tic;
    u = irntv(b, KC, lambda, pars);
    timePerf = toc;

    t2 = 0;

end


% ============================================

function[u, timePerf, t2] = runl1TV_Adapt(b, lambda, loops)


tic;

S0 = adaptMedian(b);

t2 = toc;

InputDims = size(S0);


if(length(InputDims) == 2)

    lambda_adapt  = 1e-6*(S0 == 0) + lambda*5*(S0 ~= 0);

else

    for d=1:3,
      lambda_adapt(:,:,d)  = 1e-6*(S0(:,:,d) == 0) + lambda*5.0*(S0(:,:,d) ~= 0);
    end

end

% ---------- Setup for IRN_adapt  ------------------

pars_irn_adapt = irntvInputPars('l1tv_adapt');

pars_irn_adapt.adapt_epsR   = 1;
pars_irn_adapt.epsR_cutoff  = 0.01;
pars_irn_adapt.adapt_epsF   = 1;
pars_irn_adapt.epsF_cutoff  = 0.1;

pars_irn_adapt.pcgtol_ini = 1e-4;

pars_irn_adapt.loops      = 1;

pars_irn_adapt.U0      = b;   % necessary the for adapt case


  tic;
  u = irntv(b, [], lambda_adapt, pars_irn_adapt);

  for m = 2:loops

    pars_irn_adapt.U0 = u;
    lambda_adapt = lambda_adapt*1.2;

    u = irntv(b, [], lambda_adapt, pars_irn_adapt);

  end

  timePerf = toc;



end


% ============================================

% ============================================

function[u, timePerf, t2] = run_poiTV(b, lambda, loops, KC)

    pars = irntvInputPars('tv_poisson_adapt');

    pars.pcgtol_ini = 1e-4;
    pars.loops      = 1;
%      pars.U0         = b;

    pars.adapt_epsR   = 1;
    pars.epsR_cutoff  = 0.1;
    pars.adapt_epsF   = 1;
    pars.epsF_cutoff  = 0.15; 

    pars.loops_NQP    = 35;
    pars.gamma_NQP    = 5e-3;
    pars.alpha_NQP    = 5e-1;
    pars.vmax_NQP     = 1+1;            % maximum value of output image

    T = lambda;

    tic;
    u = irntv(b+0.01, KC, lambda, pars)-0.01;

    

    for ln=2:loops,

      [sig1, cv1, mu1, siginv1, cvinv1, S1, lmu1] = estNoise(u, 2*3+1, 100, 1);
      [sig2, cv2, mu2, siginv2, cvinv2, S2, lmu2] = estNoise(lmu1, 2*3+1, 100, 1);

      pars.U0 = u;
      vmin = min( pars.U0(:) );
      if(vmin < 0) pars.U0 = pars.U0 - vmin + 0.01; end

      alpha = 1.0;
      [Nrows Ncols Ndims] = size(u);

      mask = zeros(size(S2)); 

      for d = 1:Ndims,
        maskTmp = zeros(Nrows*Ncols,1);

        vS1 = S2(:,:,d); vS1 = vS1(:);

        tau = unimodal(S1(:,:,d));
        p2 = find( vS1 < alpha*tau );
        maskTmp(p2) = (1 - Normalize( vS1(p2) ) )*(1-0.85) + 0.85;

        p3 = find(vS1 >= alpha*tau);
        maskTmp(p3) = 1.0;

        mask(:,:,d) = reshape(maskTmp, [Nrows Ncols]);
      
      end % _END_ FOR(d)

      T = T.*mask;

      u = irntv(b+0.1, KC, T, pars)-0.1;

    end

    % ====================================
    timePerf = toc;

    t2 = 0;

end

% ============================================

function[snrRec, ssimRec] = computePerf(Img, u, flagSSIM, flagTitle)

    snrRec = snr(Img, u); 
    if(flagSSIM)
      ssimRec = ssim_index(255*Normalize(Img), 255*Normalize(u));
    else
      ssimRec = NaN;
    end

    figure; imagesc( Normalize(u) );  colormap gray; axis image; axis off;
    if(flagTitle)
      title(sprintf('Restored Image\n SNR: %4.1fdB. SSIM: %1.2f ', ...
               snrRec, ssimRec));
    end

end