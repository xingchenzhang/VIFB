function[S0 IRN_Ig Adapt_Ig, Ig_L1, lambda_adapt] = eusipco11_spatially(Img)

if nargin < 1
    Img='clena';
end


SHOW_IMGS = true;
SHOW_HIST = false;

%  HEURISTICS = true;
HEURISTICS = false;


example = 'l1denoise';
%  example = 'l1deconv';

loops = 8;

%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  switch lower(Img)

    case{'lena'}
      Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

    case{'peppers'}
      Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;

    case{'goldhl.'}
      Ig = double( imread('gray_imgs/goldhill_gray.png') ) / 255;

    case{'clena'}
      Ig = double( imread('color_imgs/lena_color_512.png') ) / 255;

    case{'none'}
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting test_adaptL1 code...');
      return;

    otherwise
      error('Not a valid image\n');

  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %       Blurred & noisy images        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      kernel = fspecial('disk',3.2);
      K = @(x) imfilter(x, kernel, 'symmetric','conv');
      KT = @(x) K(x);
      KC = {K, KT};

      IgBlur = K(Ig);

      Ig_01L1 = imnoise(Ig, 'salt & pepper', 0.1);
      Ig_03L1 = imnoise(Ig, 'salt & pepper', 0.3);
      Ig_05L1 = imnoise(Ig, 'salt & pepper', 0.5);
      Ig_07L1 = imnoise(Ig, 'salt & pepper', 0.7);
      Ig_09L1 = imnoise(Ig, 'salt & pepper', 0.9);

      Ig_zoom = zeros(size(Ig));
      Ig_zoom(1:2:end,1:2:end) = Ig(1:2:end,1:2:end);


      IgBlur_01L1 = imnoise(IgBlur, 'salt & pepper', 0.1);
      IgBlur_03L1 = imnoise(IgBlur, 'salt & pepper', 0.3);
      IgBlur_05L1 = imnoise(IgBlur, 'salt & pepper', 0.5);
      IgBlur_07L1 = imnoise(IgBlur, 'salt & pepper', 0.7);
      IgBlur_09L1 = imnoise(IgBlur, 'salt & pepper', 0.9);


if( strcmp(example,'l1denoise') || strcmp(example,'all') )
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Ig_L1 = Ig_09L1;
Ig_L1 = Ig_07L1;
%  Ig_L1 = Ig_zoom;

lambda_irn    = 1.4;

% pre-processing
t =tic;
S0 = adaptMedian_colfilt(Ig_L1);
toc(t)

InputDims = size(S0);


if(length(InputDims) == 2)

  if(HEURISTICS)
    lambda_adapt  = 1e-6*(S0 == 0) + 0.85*(S0 == 3) + ...
                    1.1*(S0 == 5) + 1.2*(S0 == 7) + 1.4*(S0 == 9);
  else
    lambda_adapt  = 1e-6*(S0 == 0) + 1*(S0 > 0);
  end

else

  if(HEURISTICS)

    for d=1:3,
    lambda_adapt(:,:,d)  = 1e-6*(S0(:,:,d) == 0) + 1.05*(S0(:,:,d) == 3) + ...
                    1.3*(S0(:,:,d) == 5) + 1.4*(S0(:,:,d) == 7) + 1.6*(S0(:,:,d) == 9);
    end
  else

    for d=1:3,
      lambda_adapt(:,:,d)  = 1e-6*(S0(:,:,d) == 0) + 1.0*(S0(:,:,d) > 0);
    end
  end

end


% ---------- Setup for IRN (original)  ------------------

pars_irn = irntvInputPars('l1tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.1;

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 8;


% ---------- Setup for IRN_adapt  ------------------

pars_irn_adapt = irntvInputPars('l1tv_adapt');

pars_irn_adapt.adapt_epsR   = 1;
pars_irn_adapt.epsR_cutoff  = 0.01;
pars_irn_adapt.adapt_epsF   = 1;
pars_irn_adapt.epsF_cutoff  = 0.1;

pars_irn_adapt.pcgtol_ini = 1e-4;

pars_irn_adapt.loops      = 2;

pars_irn_adapt.U0      = Ig_L1;   % necessary the for adapt case

% ---------- 
% ---------- 

  tic;
  IRN_Ig = irntv(Ig_L1, {}, lambda_irn, pars_irn);
  tirn_denoise = toc;

  irn_snr_denoise = snr(Ig, IRN_Ig);



  tic;
  Adapt_Ig = irntv(Ig_L1, {}, lambda_adapt, pars_irn_adapt);
  tirn_denoise = toc;

  adapt_snr_denoise = snr(Ig, Adapt_Ig);


  if(SHOW_IMGS)

    if(SHOW_HIST)
    figure; hist( lambda_adapt(:), 100)
    title('Lambda 1st pass');
    end

    figure; imagesc( Normalize(IRN_Ig) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN. SNR: %4.1fdB.\n ', ...
                  irn_snr_denoise));

    figure; imagesc( Normalize(Adapt_Ig) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN Adapt. SNR: %4.1fdB.\n ', ...
                  adapt_snr_denoise));


  end


% ---------- New pass IRN_adapt  ------------------



for m = 2:loops

  Ipass = Adapt_Ig;
  dI = Ipass - Ig_L1;
  lambda_adapt = adaptLambda(dI, lambda_adapt, 0.5, S0, 0.8);

  pars_irn_adapt.U0      = Ipass;   % necessary the for adapt case

  tic;
  Adapt_Ig = irntv(Ipass, {}, lambda_adapt, pars_irn_adapt);
  tirn_denoise = toc;

  adapt_snr_denoise = snr(Ig, Adapt_Ig);


  if(SHOW_IMGS)

    if(SHOW_HIST)
    figure; hist( lambda_adapt(:), 100)
    title(sprintf('Lambda %dth pass',m));
    end

    figure; imagesc( Normalize(Adapt_Ig) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN Adapt. SNR: %4.1fdB.\n ', ...
                  adapt_snr_denoise));

  end


end

end % _END_ if( strcmp(example,'l1denoise') || strcmp(example,'all') )


% --------------------------------------------------------------------
% --------------------------------------------------------------------


if( strcmp(example,'l1deconv') || strcmp(example,'all') )
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    L1 Deconvolve        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


IgBlur_L1 = IgBlur_03L1;
Ig_L1 = IgBlur_L1;


lambda_irn    = 0.9;

% pre-processing

W0 = adaptMedian_colfilt(Ig_L1);
InputDims = size(W0);

if(length(InputDims) == 2)

  if(HEURISTICS)
    lambda_adapt  = 1e6*(W0 == 0) + (1/0.035)*(W0 == 3) + ...
                  (1/0.07)*(W0 == 5) + (1/0.1)*(W0 == 7) + (1/0.2)*(W0 == 9);
  else
    lambda_adapt  = 1e6*(W0 == 0) + 1.0*(W0 ~= 0);
  end

else

  if(HEURISTICS)
    for d=1:3,
      lambda_adapt(:,:,d)  = 1e6*(W0(:,:,d) == 0) + (1/0.15)*(W0(:,:,d) == 3) + ...
                (1/0.3)*(W0(:,:,d) == 5) + (1/0.45)*(W0(:,:,d) == 7) + (1/0.6)*(W0(:,:,d) == 9);
    end
  else
    for d=1:3, lambda_adapt(:,:,d)  = 1e6*(W0(:,:,d) == 0) + 1.0*(W0(:,:,d) ~= 0); end
  end

end


% ---------- Setup for IRN (original)  ------------------

pars_irn = irntvInputPars('l1tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 6;


% ---------- Setup for IRN_adapt  ------------------

pars_irn_adapt = irntvInputPars('l1tv_adapt');

pars_irn_adapt.adapt_epsR   = 1;
pars_irn_adapt.epsR_cutoff  = 0.01;
pars_irn_adapt.adapt_epsF   = 1;
pars_irn_adapt.epsF_cutoff  = 0.05;

pars_irn_adapt.pcgtol_ini = 1e-4;

pars_irn_adapt.loops      = 2;

pars_irn_adapt.U0      = Ig_L1;   % necessary the for adapt case

% ---------- 
% ---------- 

  tic;
  IRN_Ig = irntv(Ig_L1, KC, lambda_irn, pars_irn);
  tirn_denoise = toc;

  irn_snr_denoise = snr(Ig, IRN_Ig);



  tic;
%    Adapt_Ig = irntv(Ig_L1+0.1, KC, lambda_adapt, pars_irn_adapt) - 0.1;
  Adapt_Ig = irntv(Ig_L1, KC, lambda_adapt, pars_irn_adapt);
  tirn_denoise = toc;

  adapt_snr_denoise = snr(Ig, Adapt_Ig);


  if(SHOW_IMGS)

    if(SHOW_HIST)
    figure; hist( lambda_adapt(:), 100)
    title('Lambda 1st pass');
    end

    figure; imagesc( Normalize(IRN_Ig) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN. SNR: %4.1fdB.\n ', ...
                  irn_snr_denoise));

    figure; imagesc( Normalize(Adapt_Ig) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN Adapt. SNR: %4.1fdB.\n ', ...
                  adapt_snr_denoise));


  end


% ---------- New pass IRN_adapt  ------------------



for m = 2:loops

  Ipass = Adapt_Ig;
  dI = K(Ipass) - IgBlur_L1;
%    lambda_adapt = adaptLambda(dI, lambda_adapt, 0.5, W0, length(kernel), 1.5);

  lambda_adapt = adaptLambda(dI, lambda_adapt, 0.5, W0, 0.8);


  pars_irn_adapt.U0      = Ipass;   % necessary the for adapt case

  tic;
%    Adapt_Ig = irntv(K(Ipass)+0.1, KC, lambda_adapt, pars_irn_adapt)-0.1;
  Adapt_Ig = irntv(K(Ipass), KC, lambda_adapt, pars_irn_adapt);
  tirn_denoise = toc;

  adapt_snr_denoise = snr(Ig, Adapt_Ig);


  if(SHOW_IMGS)

    if(SHOW_HIST)
    figure; hist( lambda_adapt(:), 100)
    title(sprintf('Lambda %dth pass',m));
    end

    figure; imagesc( Normalize(Adapt_Ig) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN Adapt. SNR: %4.1fdB.\n ', ...
                  adapt_snr_denoise));

  end

end



end % _END_ if( strcmp(example,'l1deconv') || strcmp(example,'all') )

return

%==========================================================
%==========================================================
%==========================================================


%==========================================================
