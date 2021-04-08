
% Simulation code for [1], presented in EUSIPCO'11 (http://www.eusipco2011.org/)
%
% [1] Paul Rodriguez "Non-Convex Total Variation Restoration of Speckled 
%                     Vector-Valued Images via a Non-negative Quadratic 
%                     Programming algorithm"
%                     submitted to the 2011 European Signal Processing 
%                     Conference (EUSIPCO'11), Barcelona, Spain, 2011
%  
% 
%
% Legal:
%   Eusipco.m is based on NUMIPAD (http://numipad.sf.net). NUMIPAD is free 
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

%  example = 'eusipco11'
%  example = 'speckle_deconv';
example = 'speckle_denoise';

%  NOISE = 'imnoise';
NOISE = 'gamrnd';

Img = 'lena';
%  Img = 'tank';
%  Img = 'cman512';
%  Img = 'Clena';

%  Img = 'Cboats';

%  Img = 'peppers';
%  Img = 'cman';
%  Img = 'goldhl.';
%  Img = 'satell.';
%  Img = 'all';

Color = 0;  % grayscale images
%  Color = 1;  % color images

Nloops = 1; % This value will be change to 10 if example == eusipco11

nmpdef;

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


if( strcmp(example, 'eusipco11') )

  disp('  Running the "eusipco11" example...');
  disp('  This script reproduces the Exprimental Results of ');
  disp('  "Non-Convex Total Variation Restoration of Speckled');
  disp('   Vector-Valued Images via a Non-negative Quadratic');
  disp('   Programming algorithmNon-negative Quadratic Programming');
  disp('   Total Variation"');
  disp(' ')
  disp('  This script takes sometime (ten-trial average)');
  disp(' ')
  Nloops = 10;

  if( ~strcmp(Img, 'all') ) 

    disp('NOTE:');
    disp('  If you run the "eusipco11" example, you should set Img = "all" ...');
    disp('  setting Img = "all" ');
    Img = 'all';
  end

end


str_all = [];

if( strcmp(lower(Img), 'all') )

  all = 1;
  if( strcmp(example, 'eusipco11') )
    Img = {'lena', 'cman512', 'Clena' };
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

    case{'tank'}
      Ig = double( imread('gray_imgs/tank_512.tiff') ) / 255;

    case{'peppers'}
      Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;
%        Ig = double( imread('gray_imgs/peppers_gray.png') );

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
      disp('Exiting Icip11 code...');
      return;

    otherwise
      error('Not a valid image\n');
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting Icip11 code...');
      
  end % _END_ SWITCH(Img)

  else
  
  SSIM_CODE = 0;
  switch lower(Img)

    case{'clena'}
      Ig = double( imread('color_imgs/lena_color_512.png') ) / 255;

    case{'cpeppers'}
      Ig = double( imread('color_imgs/peppers_color.png') ) / 255;

    case{'cboats'}
      Ig = double( imread('color_imgs/boats_color.png') ) / 255;

    case{'cmandrill'}
      Ig = double( imread('color_imgs/mandrill_color.png') ) / 255;

    case{'none'}
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting Icip11 code...');
      return;

    otherwise
      error('Not a valid image\n');
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting Icip11 code...');

  end % _END_ SWITCH(Img)

  end % _END_ IF(color)

end

%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));

%%%%%%%%%%%%%%%%%%%%%%%%
%       kernels        %
%%%%%%%%%%%%%%%%%%%%%%%%


for k = 1:NImgs,

  if( all == 1 )
    kImg = Img{k}

    switch lower(example)

    case{'eusipco11'}

      switch lower(kImg);
  
        case{'cman512'}
          Ig = double( imread('gray_imgs/cameraman_512x512.tiff') ) / 255;
          Color = 0;
          SSIM_CODE = 1;

        case{'lena'}
          Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;
          Color = 0;
          SSIM_CODE = 1;

        case{'clena'}
          Ig = double( imread('color_imgs/lena_color_512.png') ) / 255;
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
  
        case{'clena'}
          Ig = double( imread('color_imgs/lena_color_512.png') ) / 255;
          Color = 1;
          SSIM_CODE = 0;

        case{'cpeppers'}
          Ig = double( imread('color_imgs/peppers_color.png') ) / 255;
          Color = 1;
          SSIM_CODE = 0;

        case{'cmandrill'}
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
  Color
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

%      case{'speckle_deconv', 'eusipco11'}
    case{'speckle_deconv'}

      if(strcmp(NOISE,'imnoise'))
        M = [ 0.05, 0.1 ];

        for m=1:length(M)
          IgBlur = K(Ig*M(m));
          vmax(m) = max(abs(IgBlur(:)));
          Ig_p{m} = imnoise( IgBlur/vmax(m), 'speckle', M(m) );
        end
      end

      if(strcmp(NOISE,'gamrnd'))

        if(strcmp(kImg,'cman512')) M = [ 13 3 ];
        else M = [ 33 10 5 ];
        end

        for m=1:length(M)
          IgBlur = K(Ig);
          vmax(m) = max(abs(IgBlur(:)));
          gam = gamrnd(M(m), 1/M(m), [Nrows, Ncols, Ndims]);
          Ig_p{m} = IgBlur.*gam;
        end
      end


    case{'speckle_denoise','eusipco11'}

      if(strcmp(NOISE,'imnoise'))
        M = [ 0.05, 0.1 ];
        for m=1:length(M)
          vmax(m) = max(abs(Ig(:)));
          Ig_p{m} = imnoise( Ig, 'speckle', M(m) );
        end
      end

      if(strcmp(NOISE,'gamrnd'))

        if(strcmp(kImg,'cman512')) M = [ 13 3 ];
        else M = [ 33 10 5 ];
        end

        for m=1:length(M)
          vmax(m) = max(abs(Ig(:)));
          gam = gamrnd(M(m), 1/M(m), [Nrows, Ncols, Ndims]);
          Ig_p{m} = Ig.*gam;
        end
      end


  end % _END_ switch lower(example)


      if(SHOW_IMGS)

      for m=1:length(M)

%          p_snr(m) = snr(Normalize(Ig), Normalize(Ig_p{m})); 
%          p_psnr(m) = psnr(Normalize(Ig), Normalize(Ig_p{m})); 
%          p_reErr(m) = reErr(Normalize(Ig), Normalize(Ig_p{m}));


        p_snr(m) = snr(Ig, Ig_p{m}); 
        p_psnr(m) = psnr(Ig, Ig_p{m}); 
        p_reErr(m) = reErr(Ig, Ig_p{m});

        if(SSIM_CODE)
          p_ssim(m) = ssim_index(255*Normalize(Ig), 255*Normalize(Ig_p{m}));
        else
          p_ssim(m) = NaN;
        end

        figure; imagesc( Normalize(Ig_p{m}) ); 
        colormap gray; axis image; axis off;
        title(sprintf('Noisy/Blurred Image. L %d\n SNR: %2.1fdB.  PSNR: %2.1fdB. \nSSIM: %1.3fdB, reErr: %1.3f', ...
              M(m), p_snr(m), p_psnr(m), p_ssim(m), p_reErr(m)));

      end

      end % _END_ IF(SHOW_IMGS)



%-----------------------------------------------------------------------------


  if( strcmp(example,'speckle_deconv') || strcmp(example,'all') || ( strcmp(example, 'eusipco11') ) )
%    if( strcmp(example,'speckle_deconv') || strcmp(example,'all') )

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %     Speckle Deconvolved     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noise level: max pixel value=30

pars = irntvInputPars('tv_speckle');

if(strcmp(NOISE,'gamrnd'))
  if(Color == 0) 
      if(strcmp(kImg,'cman512')) lambda  = [0.2 0.4]; 
      else lambda  = [0.15, 0.2 0.3]; 
      end
  end

  if(Color == 1) 
    lambda  = [0.125, 0.175, 0.4]; 
  end
end


  pars.loops        = 8;

  pars.adapt_epsF   = 0;
  pars.epsF_cutoff  = 1e-4;
  pars.adapt_epsR   = 1;
  pars.epsR_cutoff  = 0.05;

  pars.loops_NQP    = 50;
  pars.gamma_NQP    = 5e-3;
  pars.alpha_NQP    = 5e-1;

  pars.epsF_cutoff  = 1e-4;


for m=1:length(lambda)


  pars.M = M(m);

  tic;
  IRNNQP_Ig = irntv(Ig_p{m}+0.1, KC, lambda(m), pars)-0.1;
  tirn_deconv(m, loops) = toc;

  irn_snr(m, loops) = snr(Ig, IRNNQP_Ig); 
  irn_isnr(m, loops) = isnr(Ig, IRNNQP_Ig, Ig_p{m}); 
  irn_psnr(m, loops) = psnr(Ig, IRNNQP_Ig); 
  irn_reErr(m, loops) = reErr(Ig, IRNNQP_Ig); 

%    irn_mae(m, loops) = norm(Ig(:) - IRNNQP_Ig(:), 1) / Ndata;

  if(SSIM_CODE)
    irn_ssim(m, loops) = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig));
  else
    irn_ssim(m,loops) = NaN;
  end

  if(SHOW_IMGS)

    figure; imagesc( Normalize(IRNNQP_Ig) ); 
    colormap gray; axis image; axis off;
%      title(sprintf('Denoised Image - IRN-NQP.\n SNR: %2.1fdB. PSNR: %2.1fdB.\nISNR : %2.1fdB. SSIM: %1.3fdB. reErr: %3.3f', ...
%            irn_snr(m, loops), irn_snr(m, loops), irn_isnr(m, loops), irn_ssim(m, loops),irn_reErr(m, loops)));

  end % _END_ IF(SHOW_IMGS)

end % _END_ FOR(m)

if(loops == Nloops)

  time_mean = mean(tirn_deconv, 2);
  snr_mean = mean(irn_snr, 2);
  psnr_mean = mean(irn_psnr, 2);
  isnr_mean = mean(irn_isnr, 2);
  reErr_mean = mean(irn_reErr, 2);
  ssim_mean = mean(irn_ssim, 2);

  if(k==1) str_all = sprintf(' \n');
  else str = sprintf(' \n'); str_all = [ str_all str ];
  end

  str_deconv = sprintf('(blurred by a 7x7 out-of-focus kernel)');

  if(Color == 0) str = sprintf('Img: grayscale %s %s\n', lower(kImg), str_deconv);
  else str = sprintf('Img: color %s %s\n', lower(kImg), str_deconv);
  end

  str_all = [ str_all str ];


  str = sprintf(' M \t SNR \t SSIM \t PSNR \t ReErr    Time\n');
  str_all = [ str_all str ];

  for m=1:length(M)
    str = sprintf('%1.2f \t%1.2f \t%1.3f \t%1.2f \t%1.3f    %1.2f \n', M(m), ...
                  snr_mean(m), ssim_mean(m), psnr_mean(m), reErr_mean(m), time_mean(m) );
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


%    if( strcmp(example,'speckle_denoise') || ( strcmp(example, 'eusipco11') ) ) 
  if( strcmp(example,'speckle_denoise')  ) 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %      Speckle Denoising      %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%noise level: max pixel value=30

if(strcmp(NOISE,'imnoise'))
  if(Color == 0) 
    lambda  = [0.15, 0.2 ]; 
  end

  if(Color == 1) 
    lambda  = [0.15, 0.2 ]; 
  end
end


if(strcmp(NOISE,'gamrnd'))
  if(Color == 0) 
      if(strcmp(kImg,'cman512')) lambda  = [0.2 0.4]; 
%        else lambda  = [0.15, 0.2 0.3]; 
      else lambda  = [0.15, 0.2 0.3]; 
      end
  end

  if(Color == 1) 
    lambda  = [0.175, 0.25, 0.35]; 
  end
end

pars = irntvInputPars('tv_speckle');


  pars.adapt_epsR   = 1;
  pars.epsR_cutoff  = 0.1;
  pars.adapt_epsF   = 1;
  pars.epsF_cutoff  = 1e-4;

  pars.loops_NQP    = 40;
  pars.loops        = 10;
  pars.gamma_NQP    = 5e-3;
  pars.alpha_NQP    = 5e-1;


for m=1:length(lambda)


    off = 0.1;
    pars.speckle_variant = NMP_TV_SPECKLE_AA;
    pars.vmax_NQP = vmax(m)+off;
    alpha = 1.0;

  tic;
  IRNNQP_Ig = irntv(Ig_p{m}+off, {}, lambda(m), pars)-off;
  tirn_deconv(m, loops) = toc;

  irn_snr(m, loops) = snr(Ig, IRNNQP_Ig); 
  irn_psnr(m, loops) = psnr(Ig, IRNNQP_Ig); 
  irn_isnr(m, loops) = isnr(Ig, IRNNQP_Ig, Ig_p{m}); 
  irn_reErr(m, loops) = reErr(Ig, IRNNQP_Ig); 

  if(SSIM_CODE)
    irn_ssim(m, loops) = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig));
  else
    irn_ssim(m,loops) = NaN;
  end

  if(SHOW_IMGS)

    figure; imagesc( Normalize(IRNNQP_Ig) ); 
    colormap gray; axis image; axis off;
%      title(sprintf('Denoised Image - IRN-NQP.\n SNR: %2.1fdB. PSNR: %2.1fdB. \nISNR : %2.1fdB. SSIM: %1.3fdB. reErr: %3.3f', ...
%            irn_snr(m, loops), irn_psnr(m, loops), irn_isnr(m, loops), irn_ssim(m, loops),irn_reErr(m, loops)));

  end % _END_ IF(SHOW_IMGS)

end % _END_ FOR(m)

if(loops == Nloops)

  time_mean = mean(tirn_deconv, 2);
  snr_mean = mean(irn_snr, 2);
  psnr_mean = mean(irn_psnr, 2);
  isnr_mean = mean(irn_isnr, 2);
  reErr_mean = mean(irn_reErr, 2);
  ssim_mean = mean(irn_ssim, 2);

  if(k==1) str_all = sprintf(' \n');
  else str = sprintf(' \n'); str_all = [ str_all str ];
  end


  if(Color == 0) str = sprintf('Img: grayscale %s (denoising)\n', lower(kImg));
  else str = sprintf('Img: color %s (denoising)\n', lower(kImg));
  end

  str_all = [ str_all str ];

  str = sprintf(' M \t SNR \t SSIM \t PSNR \t ReErr    Time\n');
  str_all = [ str_all str ];

  for m=1:length(M)
    str = sprintf('%1.2f \t%1.2f \t%1.3f \t%1.2f \t%1.3f    %1.2f \n', M(m), ...
                  snr_mean(m), ssim_mean(m), psnr_mean(m), reErr_mean(m), time_mean(m) );
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

