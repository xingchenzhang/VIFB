
% Simulation code for [1], presented in ISPA'11 (http://www.isispa.org/)
%
% [1] Paul Rodriguez "Total Variation Regularization for Poisson Vector-Valued 
%                     Image Restoration with a Spatially Adaptive Regularization 
%                     Parameter Selection"
%                     accepted to 7th International Symposium on Image and Signal 
%                     Processing and Analysis (ISPA 2011)
%  
% 
%
% Legal:
%   ispa11.m is based on NUMIPAD (http://numipad.sf.net). NUMIPAD is free 
%   software, you can redistribute it and/or modify it under the terms of 
%   the GNU General Public License (version 2).
%
%   The NUMIPAD library is being developed under U.S. Government contract
%   W-7405-ENG-36 for Los Alamos National Laboratory.
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov



%  SHOW_IMGS = false;
SHOW_IMGS = true;

NATURAL = false;

%  example = 'ispa11'
example = 'poiss_deconv';
%  example = 'poiss_denoise';


%  Img = 'lena';
Img = 'peppers';
%  Img = 'cman512';
%  Img = 'goldhl.';
%  Img = 'satell.';
%  Img = 'all';
%  Img = 'dubrov1';

%  Color = 0;  % grayscale images
Color = 1;  % color images

Nloops = 1; % This value will be change to 10 if example == ispa11

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if( strcmp(example, 'ispa11') )

  disp('  Running the "ispa11" example...');
  disp('  This script reproduces the Exprimental Results of ');
  disp('  "Non-negative Quadratic Programming Total Variation');
  disp('   Regularization for Poisson Grayscale and Vector-Valued');
  disp('   Image Restoration"');
  disp(' ')
  disp('  This script takes sometime (ten-trial average)');
  disp(' ')
  Nloops = 10;

  if( ~strcmp(Img, 'all') ) 

    disp('NOTE:');
    disp('  If you run the "ispa11" example, you should set Img = "all" ...');
    disp('  setting Img = "all" ');
    Img = 'all';
  end

end


str_all = [];

if( strcmp(lower(Img), 'all') )

  all = 1;
  if( strcmp(example, 'ispa11') )
    Img = {'cman', 'peppers'};
  else
    Img = {'lena', 'peppers', 'goldhl.', 'cman', 'cman512', 'Clena', 'Cpeppers', 'Cmandrill'};
  end

  NImgs = length(Img);


else

  all = 0;
  NImgs = 1;

  if(Color == 0)

  switch lower(Img)

    case{'lena'}
      Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

    case{'peppers'}
      Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;

    case{'goldhl.'}
      Ig = double( imread('gray_imgs/goldhill_gray.png') ) / 255;

    case{'cman'}
      Ig = double( imread('gray_imgs/cameraman_256x256.tiff') ) / 255;

    case{'cman512'}
      Ig = double( imread('gray_imgs/cameraman_512x512.tiff') ) / 255;

    case{'none'}
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting ispa11 code...');
      return;

    otherwise
      error('Not a valid image\n');
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting ispa11 code...');
      
  end % _END_ SWITCH(Img)

  else

  SSIM_CODE = 0;
  switch lower(Img)

    case{'lena'}
      Ig = double( imread('color_imgs/lena_color_512.png') ) / 255;

    case{'peppers'}
      Ig = double( imread('color_imgs/peppers_color.png') ) / 255;

    case{'mandrill'}
      Ig = double( imread('color_imgs/mandrill_color.png') ) / 255;

    case{'dubrov1'}
%        Itmp = double( imread('color_imgs/dubrovnik1.png') ) / 255;
%        Ig = Itmp(1400:2000, 400:1300,: );
      Ig = double( imread('color_imgs/dubrovnik1-full.png') ) / 255;
      NATURAL = true;

    case{'none'}
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting ispa11 code...');
      return;

    otherwise
      error('Not a valid image\n');
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting ispa11 code...');

  end % _END_ SWITCH(Img)

  end % _END_ IF(color)

end

%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));


vst_poiss = @(x) ( sqrt( abs(x) + 1 ) + sqrt( abs(x) ) );

%%%%%%%%%%%%%%%%%%%%%%%%
%       kernels        %
%%%%%%%%%%%%%%%%%%%%%%%%


for k = 1:NImgs,

  if( all == 1 )
    kImg = Img{k};

    switch lower(example)

    case{'ispa11'}

      switch lower(kImg)
  
        case{'cman'}
          Ig = double( imread('gray_imgs/cameraman_256x256.tiff') ) / 255;
          Color = 0;

        case{'peppers'}
          Ig = double( imread('color_imgs/peppers_color.png') ) / 255;
          Color = 1;
          SSIM_CODE = 0;

      end % _END_ switch lower(Img)

    otherwise

      switch lower(kImg)
  
        case{'cman'}
          Ig = double( imread('gray_imgs/cameraman_256x256.tiff') ) / 255;
          Color = 0;

        case{'cman512'}
          Ig = double( imread('gray_imgs/cameraman_512x512.tiff') ) / 255;
          Color = 0;

        case{'lena'}
          Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

        case{'peppers'}
          Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;

        case{'goldhl.'}
          Ig = double( imread('gray_imgs/goldhill_gray.png') ) / 255;
  
        case{'Clena'}
          Ig = double( imread('color_imgs/lena_color_512.png') ) / 255;
          Color = 1;
          SSIM_CODE = 0;

        case{'Cpeppers'}
          Ig = double( imread('color_imgs/peppers_color.png') ) / 255;
          Color = 1;
          SSIM_CODE = 0;

        case{'Cmandrill'}
          Ig = double( imread('color_imgs/mandrill_color.png') ) / 255;
          Color = 1;
          SSIM_CODE = 0;



      end % _END_ switch lower(Img)

    end % _END_ switch lower(example)

  else

    kImg = Img;

  end %_END_ IF(all)


  disp(' ');
  disp(kImg);
  disp(' ');

  if( strcmp(lower(kImg), 'cman') )
    kernel = ones(7,7);
    kernel = kernel / sum(kernel(:));
  else
    kernel = fspecial('disk',3.2);
  end

  K = @(x) imfilter(x, kernel, 'symmetric','conv');

  KT = @(x) K(x);
  KC = {K, KT};


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %       Blurred & noisy images        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nrows Ncols Ndims] = size(Ig);
Ndata = Nrows*Ncols;

for loops=1:Nloops

  switch lower(example)

    case{'poiss_deconv', 'ispa11'}

      M = [ 5, 30, 100, 255 ];

      for m=1:length(M)
        IgBlur = K(Ig*M(m));
        vmax(m) = max(abs(IgBlur(:)));
        Ig_p{m} = poissrnd( IgBlur );
      end

    case{'poiss_denoise'}

      M = [ 5, 30, 100, 255 ];

if NATURAL
      for m=1:length(M)
        Ig_p{m} = M(m)*Ig;
      end

else
      for m=1:length(M)
        vmax(m) = max(abs(Ig(:)));
        Ig_p{m} = poissrnd( (M(m)/vmax(m))*Ig );
      end
end

  end % _END_ switch lower(example)

if(loops == 1)
  irn_snr = zeros(length(M),Nloops);
  irn_isnr = zeros(length(M),Nloops);
  irn_ssim = zeros(length(M),Nloops);
  irn_mae = zeros(length(M),Nloops);
  tirn_decon = zeros(length(M),Nloops);
end

      if(SHOW_IMGS)

      for m=1:length(M)

        p_snr(m) = snr(M(m)*Ig, Ig_p{m}); 
        p_mae(m) = norm(M(m)*Ig(:) - Ig_p{m}(:), 1) / Ndata;
        if(SSIM_CODE)
          p_ssim(m) = ssim_index(255*Normalize(Ig), 255*Normalize(Ig_p{m}));
        else
          p_ssim(m) = NaN;
        end

        figure; imagesc( Normalize(Ig_p{m}) ); 
        colormap gray; axis image; axis off;
        title(sprintf('Noisy/Blurred Image.\n SNR: %2.1fdB.  SSIM: %1.3fdB, MAE: %3.3f', ...
              p_snr(m), p_ssim(m), p_mae(m)));

      end

      end % _END_ IF(SHOW_IMGS)



%-----------------------------------------------------------------------------


  if( strcmp(example,'poiss_deconv') || strcmp(example,'all') || ( strcmp(example, 'ispa11') ) )

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %     Poisson Deconvolved     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noise level: max pixel value=30

pars = irntvInputPars('tv_poisson_adapt');

if(Color == 0) 
  SSIM_CODE=1;
  lambda  = 1./[0.25, 0.075, 0.03, 0.012 ]; 
%    lambda  = 1.5./[0.25, 0.075, 0.03, 0.012 ]; 
    update_lSigma  = [0.85 0.85 0.85 0.85];
    update_hSigma  = [0.025 0.025 0.025 0.025];

end

if(Color == 1) 
%    lambda  = 1./[0.25, 0.075, 0.03, 0.02 ]; 
  lambda  = 1./[0.1, 0.075, 0.03, 0.02 ]; 
    update_lSigma  = [0.85 0.85 0.9 0.9];
    update_hSigma  = [0.025 0.025 0.025 0.025];
end


  pars.loops        = 1;

  pars.adapt_epsF   = 0;
  pars.epsF_cutoff  = 1e-4;
  pars.adapt_epsR   = 1;
  pars.epsR_cutoff  = 0.05;

  pars.loops_NQP    = 45;
  pars.gamma_NQP    = 5e-3;
  pars.alpha_NQP    = 5e-1;

for m=1:length(lambda)

if(Color == 0)
  pars.U0 = conv2( Ig_p{m}, ones(3,3)/9, 'same');
else
  pars.U0 = zeros(size(Ig_p{m}));
  pars.U0(:,:,1) = conv2( Ig_p{m}(:,:,1), ones(3,3)/9, 'same');
  pars.U0(:,:,2) = conv2( Ig_p{m}(:,:,2), ones(3,3)/9, 'same');
  pars.U0(:,:,3) = conv2( Ig_p{m}(:,:,3), ones(3,3)/9, 'same');
end

  tic;

  T = lambda(m);

  IRNNQP_Ig = irntv(Ig_p{m}+0.1, KC, T, pars)-0.1;

  for ln=2:6,

    pars.U0 = IRNNQP_Ig;
    vmin = min( pars.U0(:) );
    if(vmin < 0) pars.U0 = pars.U0 - vmin + 0.01; end

    [sig1, cv1, mu1, siginv1, cvinv1, S1, lmu1] = estNoise(IRNNQP_Ig, 2*2+1, 100, 1);
%      [sig2, cv2, mu2, siginv2, cvinv2, S2, lmu2] = estNoise(lmu1, 2*2+1, 100, 1);


    mask = zeros(size(S1));

    for d = 1:Ndims,
      maskTmp = zeros(Nrows*Ncols,1);

%        vS1 = S2(:,:,d); vS1 = vS1(:);
      vS1 = S1(:,:,d); vS1 = vS1(:);

if 0
      alpha = 1.0;

      tau = unimodal(S1(:,:,d));
      p2 = find( vS1 < alpha*tau );
      maskTmp(p2) = (1 - Normalize( vS1(p2) ) )*(1-update_lSigma(m)) + update_lSigma(m);

      p3 = find(vS1 >= alpha*tau);
      maskTmp(p3) = 1.0;


%        p1 = find(vS1 < sig1(d));
%        maskTmp(p1) = Normalize( vS1(p1) )*(1-update_lSigma(m))*0.5 + update_lSigma(m);
%  
%        tau = unimodal(S1(:,:,d));
%        p2 = find( (vS1 >= sig1(d)).*(vS1 < alpha*tau) );
%        maskTmp(p2) = (1 - Normalize( vS1(p2) ) )*(1-update_lSigma(m))*0.5 + update_lSigma(m);
%  
%        p3 = find(vS1 >= alpha*tau);
%        maskTmp(p3) = Normalize( vS1(p3) )*update_hSigma(m) + 1.0;

else

%        alpha = 2.0;
      tau = unimodal(S1(:,:,d));
      alpha = tau/sig1(d);

      p1 = find(vS1 < sig1(d)*alpha);
      maskTmp(p1) = (Normalize( vS1(p1)  )*0.15 + 0.85);

      p2 = find(vS1 >= sig1(d)*alpha);
%        maskTmp(p2) = Normalize( vS1(p2) )*0.025 + 1.0;
      maskTmp(p2) = Normalize( vS1(p2) )*.25 + 1.0;
end

%        p1 = find(vS1 < sig1(d)*alpha);
%        maskTmp(p1) = Normalize( vS1(p1)  )*0.15 + 0.85;
%  
%        p2 = find(vS1 >= sig1(d)*alpha);
%        maskTmp(p2) = Normalize( vS1(p2) )*0.025 + 1.0;

      mask(:,:,d) = reshape(maskTmp, [Nrows Ncols]);
    end

%      mask = zeros(size(S2)); mask = mask(:);
%      vS1 = S2(:);
%      p1 = find(vS1 < sig1*alpha);
%      mask(p1) = Normalize( vS1(p1)  )*0.2 + 0.8;
%  
%      p2 = find(vS1 >= sig1*alpha);
%      mask(p2) = Normalize( vS1(p2) )*0.05 + 1.0;
%  
%      mask = reshape(mask, size(S2) );

    T = T.*mask;

    IRNNQP_Ig = irntv(Ig_p{m}+0.01, KC, T, pars)-0.01;

  end

%    figure; imagesc( 1./T ); colormap gray; axis image; axis off
  tirn_deconv(m, loops) = toc;


  irn_snr(m, loops) = snr(M(m)*Ig, IRNNQP_Ig); 
  irn_isnr(m, loops) = isnr(M(m)*Ig, IRNNQP_Ig, Ig_p{m}); 

%    irn_mae(m, loops) = norm(Ig(:)*M(m) - IRNNQP_Ig(:), 1) / Ndata;

        if(SSIM_CODE)
          irn_mae(m, loops) = ssim_index(255*Normalize(Ig*M(m)), 255*Normalize(IRNNQP_Ig));
        else
          irn_mae(m, loops) = NaN;
        end


  if(SSIM_CODE)
    irn_ssim(m, loops) = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig));
  else
    irn_ssim(m,loops) = NaN;
  end

  if(SHOW_IMGS)

    figure; imagesc( Normalize(IRNNQP_Ig) ); 
    colormap gray; axis image; axis off;
%      title(sprintf('Denoised Image - IRN-NQP.\n SNR: %2.1fdB. ISNR : %2.1fdB. SSIM: %1.3fdB. MAE: %3.3f', ...
%            irn_snr(m, loops),  irn_isnr(m, loops), irn_ssim(m, loops),irn_mae(m, loops)));

  end % _END_ IF(SHOW_IMGS)

end % _END_ FOR(m)

if(loops == Nloops)

  mae_mean  = mean(irn_mae,2);
  time_mean = mean(tirn_deconv, 2);
  snr_mean = mean(irn_snr, 2);
  isnr_mean = mean(irn_isnr, 2);

  if(k==1) str_all = sprintf(' \n');
  else str = sprintf(' \n'); str_all = [ str_all str ];
  end

  if(strcmp( lower(kImg), 'cman' ) )
    str_deconv = sprintf('(blurred by a 7x7 uniform filter)');
  else
    str_deconv = sprintf('(blurred by a 7x7 out-of-focus kernel)');
  end

  if(Color == 0) str = sprintf('Img: grayscale %s %s\n', lower(kImg), str_deconv);
  else str = sprintf('Img: color %s %s\n', lower(kImg), str_deconv);
  end

  str_all = [ str_all str ];

  str = sprintf(' M \t MAE \t SNR \t ISNR    Time\n');
  str_all = [ str_all str ];

  for m=1:length(M)
    str = sprintf('%d \t%1.2f \t%1.2f \t%1.2f    %1.2f \n', M(m), mae_mean(m), snr_mean(m), isnr_mean(m), time_mean(m) );
    str_all = [ str_all str ];
  end

  str = sprintf(' \n');
  str_all = [ str_all str ];

  if( (all==1) && (k==NImgs) ) disp(str_all); end
  if( all==0 ) disp(str_all); end


end % _END_ IF(Nloops)

%---------------------------------------------------



  end % _END_ if( strcmp(example,'l1deconv') || strcmp(example,'all') )


%-----------------------------------------------------------------------------


  if( strcmp(example,'poiss_denoise') ) 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %      Poisson Denoising      %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%noise level: max pixel value=30


if(Color == 0) 
  SSIM_CODE = 1;
  lambda  = 2./[0.35, 0.15, 0.075, 0.025 ]; 
    update_lSigma  = [0.8 0.85 0.9 0.9];
    update_hSigma  = [0.025 0.025 0.025 0.025];
end

if(Color == 1) 
  SSIM_CODE = 0;
  lambda  = 2./[0.4, 0.25, 0.15, 0.075 ]; 
    update_lSigma  = [0.8 0.85 0.9 0.9];
    update_hSigma  = [0.025 0.025 0.025 0.025]/2;
end

pars = irntvInputPars('tv_poisson_adapt');


  pars.adapt_epsR   = 1;
  pars.epsR_cutoff  = 0.1;
  pars.adapt_epsF   = 1;
  pars.epsF_cutoff  = 0.15;

  pars.loops_NQP    = 45;
  pars.loops        = 1;
  pars.gamma_NQP    = 5e-3;
  pars.alpha_NQP    = 5e-1;


for m=1:length(lambda)

%      for d = 1:Ndims,
%        pars.U0(:,:,d) = conv2( Ig_p{m}(:,:,d), ones(2,2)/4, 'same');
%      end

  pars.U0 = [];

  tic;

  T = lambda(m);
  IRNNQP_Ig = irntv(Ig_p{m}+0.01, {}, T, pars)-0.01;
    snr(M(m)*Ig, IRNNQP_Ig)


  for ln=2:8,



    [sig1, cv1, mu1, siginv1, cvinv1, S1, lmu1] = estNoise(IRNNQP_Ig, 2*3+1, 100, 1);
    [sig2, cv2, mu2, siginv2, cvinv2, S2, lmu2] = estNoise(lmu1, 2*3+1, 100, 1);

%      if(ln==2)
%        pars.U0 = lmu2;
%        IRNNQP_Ig = irntv(Ig_p{m}+0.1, {}, T, pars)-0.1;
%        snr(M(m)*Ig, IRNNQP_Ig)
%        continue
%      else
%        pars.U0 = IRNNQP_Ig;
%      end

      pars.U0 = IRNNQP_Ig;
      vmin = min( pars.U0(:) );
      if(vmin < 0) pars.U0 = pars.U0 - vmin + 0.01; end

    alpha = 1.0;
    [Nrows Ncols Ndims] = size(IRNNQP_Ig);


    mask = zeros(size(S2)); 

    for d = 1:Ndims,
      maskTmp = zeros(Nrows*Ncols,1);

      vS1 = S2(:,:,d); vS1 = vS1(:);

if 1
%        p1 = find(vS1 < sig1(d));
%        maskTmp(p1) = update_lSigma(m);
%  
%        tau = unimodal(S1(:,:,d));
%        p2 = find( (vS1 >= sig1(d)).*(vS1 < alpha*tau) );
%        maskTmp(p2) = (1 - Normalize( vS1(p2) ) )*(1-update_lSigma(m)) + update_lSigma(m);
%  
%        p3 = find(vS1 >= alpha*tau);
%        maskTmp(p3) = Normalize( vS1(p3) )*update_hSigma(m) + 1.0;


%---------------------------------------------------------------


      tau = unimodal(S1(:,:,d));
      p2 = find( vS1 < alpha*tau );
      maskTmp(p2) = (1 - Normalize( vS1(p2) ) )*(1-update_lSigma(m)) + update_lSigma(m);

      p3 = find(vS1 >= alpha*tau);
      maskTmp(p3) = 1.0;


else

      p1 = find(vS1 < sig1(d)*alpha);
      maskTmp(p1) = (1 - Normalize( vS1(p1) ) )*(1-update_lSigma(m)) + update_lSigma(m);
%        maskTmp(p1) = Normalize( vS1(p1)  )*0.15 + 0.85;

      p2 = find(vS1 >= sig1(d)*alpha);
      maskTmp(p2) = Normalize( vS1(p2) )*update_hSigma(m) + 1.0;
%        maskTmp(p2) = Normalize( vS1(p2) )*0.025 + 1.0;
%        maskTmp(p2) = (1-Normalize( vS1(p2) ))*0.025 + 0.975;
end

      mask(:,:,d) = reshape(maskTmp, [Nrows Ncols]);
      
    end

    T = T.*mask;


    IRNNQP_Ig = irntv(Ig_p{m}+0.1, {}, T, pars)-0.1;

  end

%    figure; imagesc( Normalize(1./T) ); colormap gray; axis image;
  tirn_deconv(m, loops) = toc;
  Ts{m} = T;


  irn_snr(m, loops) = snr(M(m)*Ig, IRNNQP_Ig); 
  irn_isnr(m, loops) = isnr(M(m)*Ig, IRNNQP_Ig, Ig_p{m}); 

%    irn_mae(m, loops) = norm(Ig(:)*M(m) - IRNNQP_Ig(:), 1) / Ndata;

        if(SSIM_CODE)
          irn_mae(m, loops) = ssim_index(255*Normalize(Ig*M(m)), 255*Normalize(IRNNQP_Ig));
        else
          irn_mae(m, loops) = NaN;
        end

  if(SSIM_CODE)
    irn_ssim(m, loops) = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig));
  else
    irn_ssim(m,loops) = NaN;
  end

  if(SHOW_IMGS)

    figure; imagesc( Normalize(IRNNQP_Ig) ); 
    colormap gray; axis image; axis off;
%      title(sprintf('Denoised Image - IRN-NQP.\n SNR: %2.1fdB. ISNR : %2.1fdB. SSIM: %1.3fdB. MAE: %3.3f', ...
%            irn_snr(m, loops),  irn_isnr(m, loops), irn_ssim(m, loops),irn_mae(m, loops)));

  end % _END_ IF(SHOW_IMGS)

end % _END_ FOR(m)

if(loops == Nloops)

  mae_mean  = mean(irn_mae,2);
  time_mean = mean(tirn_deconv, 2);
  snr_mean = mean(irn_snr, 2);
  isnr_mean = mean(irn_isnr, 2);

  if(k==1) str_all = sprintf(' \n');
  else str = sprintf(' \n'); str_all = [ str_all str ];
  end


  if(Color == 0) str = sprintf('Img: grayscale %s (denoising)\n', lower(kImg));
  else str = sprintf('Img: color %s (denoising)\n', lower(kImg));
  end

  str_all = [ str_all str ];

  str = sprintf(' M \t MAE \t SNR \t ISNR    Time\n');
  str_all = [ str_all str ];

  for m=1:length(M)
    str = sprintf('%d \t%1.2f \t%1.2f \t%1.2f    %1.2f \n', M(m), mae_mean(m), snr_mean(m), isnr_mean(m), time_mean(m) );
    str_all = [ str_all str ];
  end

  str = sprintf(' \n');
  str_all = [ str_all str ];

  if( (all==1) && (k==NImgs) ) disp(str_all); end
  if( all==0 ) disp(str_all); end


end % _END_ IF(Nloops)



  end % _END_ if( strcmp(example,'l1deconv') )


end % _END_ FOR(loops=1:Nloops)


end % _END_ FOR(k=1:Nimgs)

