
% Simulation code for [1], presented in ICIP'10 (http://www.icip2010.org/)
%
% [1] Paul Rodriguez "Multiplicative Updates Algorithm to Minimize the 
%     Generalized Total Variation Functional with a Non-negativity Constraint"
%     in Proceedings of IEEE International Conference on Image Processing (ICIP),
%     (Hong Kong), pp. 2509--2512, Sep. 2010
%  
% 
%
% Legal:
%   Icip10.m is based on NUMIPAD (http://numipad.sf.net). NUMIPAD is free 
%   software, you can redistribute it and/or modify it under the terms of 
%   the GNU General Public License (version 2).
%
%   The NUMIPAD library is being developed under U.S. Government contract
%   W-7405-ENG-36 for Los Alamos National Laboratory.
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov



SHOW_IMGS = false;
%  SHOW_IMGS = true;


%  example = 'l1deconv';
%  example = 'l2deconv';
%  example = 'l1denoise';
%  example = 'l2denoise';
example = 'icip10'
%  example = 'all';


Img = 'lena';
%  Img = 'peppers';
%  Img = 'goldhl.';
%  Img = 'satell.';
%  Img = 'all';


TVBC_CODE = exist('irntvc');
if( (TVBC_CODE == 0) && (strcmp(example, 'icip10')) )
  disp('NOTE:');
  disp('  The function irntvc (code for [TVBC]) is not in your path...');
  disp(sprintf('  disabling TVBC simulations for %s.\n',example));
  disp('  [TVBC] R. Chartrand and B. Wohlberg');
  disp('         "Total-variation regularization with bound constraints" ');
  disp('         Proceedings of ICASSP 2010, pp. 766-769, Mar 2010');
end


NNCGM_CODE = exist('run_cgm');
if( (NNCGM_CODE == 0) && (strcmp(example, 'l2deconv')) && (strcmp(example, 'icip10')) )
  disp('NOTE:');
  disp('  The function run_cgm (code for [NNCGM]) is not in your path...');
  disp(sprintf('  disabling NNCGM simulations for %s.\n',example));
  disp('  [NNCGM]  D. Krishnan, P. Lin, and A. Yip');
  disp('         "A primal-dual active-set method for non-negativity');
  disp('         constrained total variation deblurring problems" ');
  disp('         IEEE Transactions on Image Processing, 2007, 16:11(2766-2777).');
end

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

if( strcmp(example, 'icip10') && ~strcmp(Img, 'all') ) 

  disp('NOTE:');
  disp('  If you run the "icip10" example, you should set Img = "all" ...');
  disp('  setting Img = "all" ');
  Img = 'all';
end


if( strcmp(lower(Img), 'all') )

  all = 1;
  if( strcmp(example, 'icip10') )
    Img = {'lena', 'satell.'};
  else
    Img = {'lena', 'peppers', 'goldhl.', 'satell.'};
  end

  NImgs = length(Img);
  str_all = [];

else

  all = 0;
  NImgs = 1;

  switch lower(Img)

    case{'lena'}
      Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

    case{'peppers'}
      Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;

    case{'goldhl.'}
      Ig = double( imread('gray_imgs/goldhill_gray.png') ) / 255;

    case{'satell.'}
      Ig = double( imread('gray_imgs/satellite_128x128.tiff') ) / 255;

    case{'none'}
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting Icip10 code...');
      return;

    otherwise
      error('Not a valid image\n');

  end

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
    kImg = Img{k};
    switch lower(kImg)

      case{'lena'}
        Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

      case{'peppers'}
        Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;

      case{'goldhl.'}
        Ig = double( imread('gray_imgs/goldhill_gray.png') ) / 255;

      case{'satell.'}
        Ig = double( imread('gray_imgs/satellite_128x128.tiff') ) / 255;

    end % _END_ switch lower(Img)

  else

    kImg = Img;

  end

  disp(' ');
  disp(kImg);
  disp(' ');

  if( strcmp(lower(kImg), 'satell.') )
    kernel = fspecial('gaussian', [9 9], 20);
  else
    kernel = fspecial('disk',3.2);
  end

  K = @(x) imfilter(x, kernel, 'symmetric','conv');

  KT = @(x) K(x);
  KC = {K, KT};



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %       Blurred & noisy images        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  switch lower(example)

    case{'l1deconv'}
      IgBlur = K(Ig);
      Ig_01L1 = imnoise(IgBlur, 'salt & pepper', 0.1);
      Ig_03L1 = imnoise(IgBlur, 'salt & pepper', 0.3);

    case{'l1denoise'}
      Ig_01L1 = imnoise(Ig, 'salt & pepper', 0.1);
      Ig_03L1 = imnoise(Ig, 'salt & pepper', 0.3);

    case{'l2deconv'}
      IgBlur = K(Ig);                                             
      if( strcmp(lower(kImg), 'satell.') )
        nncgm_sig = 4.291670778829401/236; % to match nncgm example
        Ig_nnL2 = imnoise(IgBlur, 'gaussian', 0, (nncgm_sig*nncgm_sig)*max(IgBlur(:)) ); 
      else
        Ig_01L2 = imnoise(IgBlur, 'gaussian', 0, 0.01*max(IgBlur(:)) ); 
        Ig_001L2 = imnoise(IgBlur,'gaussian', 0, 0.0001*max(IgBlur(:)) );
      end

    case{'l2denoise'}
      Ig_01L2 = imnoise(Ig,'gaussian', 0, 0.01*max(Ig(:)) );
      Ig_005L2 = imnoise(Ig,'gaussian', 0, 0.0025*max(Ig(:)) ); 

    case{'icip10'}
      IgBlur = K(Ig);                                             
      if( strcmp(lower(kImg), 'satell.') )
        nncgm_sig = 4.291670778829401/236; % to match nncgm example
        Ig_nnL2 = imnoise(IgBlur, 'gaussian', 0, (nncgm_sig*nncgm_sig)*max(IgBlur(:)) ); 
      else
        Ig_01L1 = imnoise(IgBlur, 'salt & pepper', 0.1);
        Ig_03L1 = imnoise(IgBlur, 'salt & pepper', 0.3);
      end


  end % _END_ switch lower(example)



%-----------------------------------------------------------------------------


  if( strcmp(example,'l1deconv') || strcmp(example,'all') || ( strcmp(example, 'icip10') && strcmp(kImg, 'lena') ) )

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %       L1 Deconvolved        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %noise level: 0.1

  lambda = 0.45;

  pars = irntvInputPars('l1tv_nqp');

  pars.adapt_epsR   = 1;
  pars.epsR_cutoff  = 0.1;
  pars.adapt_epsF   = 1;
  pars.epsF_cutoff  = 0.15;

  pars.loops_NQP    = 25;
  pars.loops        = 10;
  pars.gamma_NQP    = 5e-3;
  pars.alpha_NQP    = 5e-1;

  
  tic;
  % The NQP solver need values diferent from zero to consider 
  % them as part of the solution (note the "+1" and "-1").
  IRNNQP_Ig_01L1 = irntv(Ig_01L1+1, KC, lambda, pars)-1;
  tirn_01L1 = toc;

  irn_01L1_snr = snr(Ig, IRNNQP_Ig_01L1); 

  if(SSIM_CODE)
    irn_01L1_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig_01L1));
  else
    irn_01L1_ssim = NaN;
  end

  if(SHOW_IMGS)

    figure; imagesc( Normalize(IRNNQP_Ig_01L1) ); 
    colormap gray; axis image; axis off;
    title(sprintf('Deconvolved Image - IRN-NQP.\n SNR: %2.1fdB.  SSIM: %1.3fdB.', ...
          irn_01L1_snr, irn_01L1_ssim));

  end % _END_ IF(SHOW_IMGS)

  
  if(TVBC_CODE)
  %-- TVBC (Total Variation Boundary Condition)

    lambda = 0.45;
    beta   = 0.1;

    tic;
    TVBC_Ig_01L1 = irntvc(Ig_01L1, KC, lambda, beta, 1, 1, 0, Inf, 1e-4, 1e-4, 10, Ig_01L1);
    ttvbc_01L1 = toc;

    tvbc_01L1_snr  = snr(Ig, TVBC_Ig_01L1);

    if(SSIM_CODE)
      tvbc_01L1_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(TVBC_Ig_01L1));
    else
      tvbc_01L1_ssim = 0;
    end

    if(SHOW_IMGS)
      figure; imagesc( Normalize(TVBC_Ig_01L1) ); 
      colormap gray; axis image; axis off;
      title(sprintf('Deconvolved Image - TV-BC.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
               tvbc_01L1_snr, tvbc_01L1_ssim));
    end

  else % IF(TVBC_CODE)

    tvbc_01L1_snr   = NaN;
    tvbc_01L1_ssim  = NaN;
    ttvbc_01L1      = NaN

  end % _END_ IF(TVBC_CODE)


%------------------
%------------------


%noise level: 0.3

  lambda = 0.95;

  pars = irntvInputPars('l1tv_nqp');

  pars.adapt_epsR   = 1;
  pars.epsR_cutoff  = 0.01;
  pars.adapt_epsF   = 1;
  pars.epsF_cutoff  = 0.05;

  pars.loops_NQP    = 25;
  pars.loops        = 10;
  pars.gamma_NQP    = 5e-3;
  pars.alpha_NQP    = 5e-1;

  tic;
  % The NQP solver need values diferent from zero to consider 
  % them as part of the solution (note the "+1" and "-1").
  IRNNQP_Ig_03L1 = irntv(Ig_03L1+1, KC, lambda, pars)-1;
  tirn_03L1 = toc;


  irn_03L1_snr = snr(Ig, IRNNQP_Ig_03L1);

  if(SSIM_CODE)
    irn_03L1_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig_03L1));
  else
    irn_03L1_ssim = 0;
  end


  if(SHOW_IMGS)

    figure; imagesc( Normalize(IRNNQP_Ig_03L1) ); 
    colormap gray; axis image; axis off;
    title(sprintf('Deconvolved Image - IRN-NQP.\n SNR: %2.1fdB.  SSIM: %1.3fdB.', ...
               irn_03L1_snr, irn_03L1_ssim));
  end


  if(TVBC_CODE)
  %-- TVBC (Total Variation Boundary Condition)

    lambda = 0.9;
    beta   = 0.1;

    tic;
    TVBC_Ig_03L1 = irntvc(Ig_03L1, KC, lambda, beta, 1, 1, 0, Inf, 1e-4, 1e-4, 10, Ig_03L1);
    ttvbc_03L1 = toc;

    tvbc_03L1_snr  = snr(Ig, TVBC_Ig_03L1);

    if(SSIM_CODE)
      tvbc_03L1_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(TVBC_Ig_03L1));
    else
      tvbc_03L1_ssim = 0;
    end

    if(SHOW_IMGS)
      figure; imagesc( Normalize(TVBC_Ig_03L1) ); 
      colormap gray; axis image; axis off;
      title(sprintf('Deconvolved Image - TV-BC.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
               tvbc_03L1_snr, tvbc_03L1_ssim));
    end

  else % IF(TVBC_CODE)

    tvbc_03L1_snr   = NaN;
    tvbc_03L1_ssim  = NaN;
    ttvbc_03L1      = NaN;

  end % _END_ IF(TVBC_CODE)




if(all == 0)

disp('l1 Deconvolve TV ');
disp(' ');
disp('                    SNR (db)          Time (s)             SSIM ');
disp(' Img    Noise     IRN-NQP  TVBC    IRN-NQP  TVBC     IRN-NQP     TVBC');
disp(' ');
disp(sprintf('%s\t 0.1       %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f', ...
      kImg, ...
      irn_01L1_snr, tvbc_01L1_snr, tirn_01L1, ttvbc_01L1, ...
      irn_01L1_ssim, tvbc_01L1_ssim));
disp(sprintf('        0.3       %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f', ...
      irn_03L1_snr, tvbc_03L1_snr, tirn_03L1, ttvbc_03L1, ...
      irn_03L1_ssim, tvbc_03L1_ssim));

else


  if( (k==1) || ( strcmp(example, 'icip10') && strcmp(kImg, 'lena') ) )

    str_all = [ str_all sprintf('\n') ];
    str_all = [ str_all sprintf('\n') ];
    str_all = [ str_all sprintf('l1 Deconvolve TV \n') ];
    str_all = [ str_all sprintf('\n') ];
    str_all = [ str_all sprintf('                    SNR (db)          Time (s)             SSIM\n') ];
    str_all = [ str_all, sprintf(' Img    Noise     IRN-NQP  TVBC    IRN-NQP  TVBC     IRN-NQP     TVBC\n') ];
    str_all = [ str_all, sprintf('\n') ];

  end

    str_all = [ str_all, sprintf('%s\t 0.1       %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f\n', ...
      Img{k}, ...
      irn_01L1_snr, tvbc_01L1_snr, tirn_01L1, ttvbc_01L1, ...
      irn_01L1_ssim, tvbc_01L1_ssim) ];

    str_all = [ str_all sprintf('    \t 0.3       %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f\n', ...
      irn_03L1_snr, tvbc_03L1_snr, tirn_03L1, ttvbc_03L1, ...
      irn_03L1_ssim, tvbc_03L1_ssim) ];


  if( k==NImgs ) 
      sprintf(str_all)
  end

end % _END_ IF(All)

end % _END_ if( strcmp(example,'l1deconv') || strcmp(example,'all') )

%-----------------------------------------------------------------------------


if( strcmp(example,'l2deconv') || strcmp(example,'all') || ( strcmp(example, 'icip10') && strcmp(kImg, 'satell.') ) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L2 Deconvolved        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  if( strcmp(kImg, 'satell.') ) % to match NNCGM example


    lambda = 8.5e-4;

    pars = irntvInputPars('l2tv_nqp');

    pars.adapt_epsR   = 1;
    pars.epsR_cutoff  = 0.01;
    pars.adapt_epsF   = 1;
    pars.epsF_cutoff  = 0.05;

    pars.loops_NQP    = 25;
    pars.loops        = 15;
    pars.gamma_NQP    = 1e-3;
    pars.alpha_NQP    = 5e-1;

    vmin = min(Ig_nnL2(:));
    if(vmin>=0) vmin = 0;
    else vmin = abs(vmin) + 1e-2; end

    tic
    % The NQP solver need values diferent from zero to consider 
    % them as part of the solution (note the "+vmin" and "-vmin").
    IRNNQP_Ig_nnL2 = irntv(Ig_nnL2+vmin, KC, lambda, pars) - vmin;
    tirn_nnL2 = toc;

    irn_nnL2_snr = snr(Ig, IRNNQP_Ig_nnL2);
    
    if(SSIM_CODE)
      irn_nnL2_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig_nnL2))
    else
      irn_nnL2_ssim = NaN;
    end

    if(SHOW_IMGS)
      figure; imagesc( Normalize(IRNNQP_Ig_nnL2) ); 
      colormap gray; axis image; axis off;
      title(sprintf('Deconvolved Image - IRN-NQP.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
            irn_nnL2_snr, irn_nnL2_ssim));
    end


    if(TVBC_CODE) %-- TVBC (Total Variation Boundary Condition)
   
      lambda = 8.5e-4;
      beta   = 1e-3;

      tic;
      TVBC_Ig_nnL2 = irntvc(Ig_nnL2, KC, lambda, beta, 2, 1, 0, Inf, 1e-4, 1e-4, 10, Ig_nnL2);
      ttvbc_nnL2 = toc;

      tvbc_nnL2_snr  = snr(Ig, TVBC_Ig_nnL2);
      tvbc_nnL2_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(TVBC_Ig_nnL2));

      if(SHOW_IMGS)
        figure; imagesc( Normalize(TVBC_Ig_nnL2) ); 
        colormap gray; axis image; axis off;
        title(sprintf('Deconvolved Image - TV-BC.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
               tvbc_nnL2_snr, tvbc_nnL2_ssim));
      end

    else % IF(TVBC_CODE)

      tvbc_nnL2_snr   = NaN;
      tvbc_nnL2_ssim  = NaN;
      ttvbc_nnL2      = NaN;

    end % _END_ IF(TVBC_CODE)


    if(NNCGM_CODE) %-- NNCGM (Non-Negatively Constrained Chan-Golub-Mulet algorithm)

      run_cgm

    else % if(NNCGM_CODE)

      nncgm_nnL2_snr   = NaN;
      nncgm_nnL2_ssim  = NaN;
      tnncgm_nnL2      = NaN;

    end % _END_ if(NNCGM_CODE)


    if(all == 0)

      disp('l2 Deconvolve TV ');
      disp(' ');
      disp('                    SNR (db)          Time (s)             SSIM ');
      disp(' Img    Noise     IRN-NQP  TVBC    IRN-NQP  TVBC     IRN-NQP     TVBC');
      disp(' ');
      disp(sprintf('%s\t %1.2f      %5.1f   %5.1f   %5.1f   %5.1f     %5.6f   %5.6f', ...
        kImg, nncgm_sig, ...
        irn_nnL2_snr, tvbc_nnL2_snr, tirn_nnL2, ttvbc_nnL2, ...
        irn_nnL2_ssim, tvbc_nnL2_ssim));

    else

      if( (k==1) || strcmp(example, 'icip10') )

        str_all = [ str_all sprintf('\n') ];
        str_all = [ str_all sprintf('\n') ];
        str_all = [ str_all sprintf('l2 Deconvolve TV \n') ];
        str_all = [ str_all sprintf('\n') ];
        str_all = [ str_all sprintf('                    SNR (db)          Time (s)             SSIM\n') ];
        str_all = [ str_all, sprintf(' Img    Noise     IRN-NQP  TVBC    IRN-NQP  TVBC     IRN-NQP     TVBC\n') ];
        str_all = [ str_all, sprintf('\n') ];

      end

      str_all = [ str_all, sprintf('%s\t %1.2f      %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f\n', ...
        Img{k}, nncgm_sig, ...
        irn_nnL2_snr, tvbc_nnL2_snr, tirn_nnL2, ttvbc_nnL2, ...
        irn_nnL2_ssim, tvbc_nnL2_ssim) ];

      if( k==NImgs ) 
        sprintf(str_all)
      end

    end

  % -------------------------
  else % not 'satellite' case
  % -------------------------

    %noise level: 0.01

    lambda = 1e-3;

    pars = irntvInputPars('l2tv_nqp');

    pars.adapt_epsR   = 1;
    pars.epsR_cutoff  = 0.01;
    pars.adapt_epsF   = 1;
    pars.epsF_cutoff  = 0.05;

    pars.loops_NQP    = 25;
    pars.loops        = 15;
    pars.gamma_NQP    = 1e-3;
    pars.alpha_NQP    = 5e-1;

    vmin = min(Ig_001L2(:));
    if(vmin>=0) vmin = 0;
    else vmin = abs(vmin) + 1e-2; end

    tic
    % The NQP solver need values diferent from zero to consider 
    % them as part of the solution (note the "+vmin" and "-vmin").
    IRNNQP_Ig_001L2 = irntv(Ig_001L2+vmin, KC, lambda, pars) - vmin;
    tirn_001L2 = toc;

    irn_001L2_snr = snr(Ig, IRNNQP_Ig_001L2);

    if(SSIM_CODE)
      irn_001L2_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig_001L2));
    else
      irn_001L2_ssim = NaN;
    end

    if(SHOW_IMGS)
      figure; imagesc( Normalize(IRNNQP_Ig_001L2) ); 
      colormap gray; axis image; axis off;
      title(sprintf('Deconvolved Image - IRN-NQP.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
            irn_001L2_snr, irn_001L2_ssim ));
    end

    % -------------------

    %noise level: 0.1

    lambda = 5e-2;

    pars = irntvInputPars('l2tv_nqp');

    pars.adapt_epsR   = 1;
    pars.epsR_cutoff  = 0.01;
    pars.adapt_epsF   = 1;
    pars.epsF_cutoff  = 0.05;

    pars.loops_NQP    = 25;
    pars.loops        = 15;
    pars.gamma_NQP    = 1e-3;
    pars.alpha_NQP    = 5e-1;

    vmin = min(Ig_01L2(:));
    if(vmin>=0) vmin = 0;
    else vmin = abs(vmin) + 1e-2; end

    tic
    % The NQP solver need values diferent from zero to consider 
    % them as part of the solution (note the "+vmin" and "-vmin").
    IRNNQP_Ig_01L2 = irntv(Ig_01L2+vmin, KC, lambda, pars) - vmin;
    tirn_01L2 = toc;

    irn_01L2_snr = snr(Ig, IRNNQP_Ig_01L2);
    if(SSIM_CODE)
      irn_01L2_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(IRNNQP_Ig_01L2));
    else
      irn_01L2_ssim = NaN;
    end

    if(SHOW_IMGS)
      figure; imagesc( Normalize(IRNNQP_Ig_01L2) ); 
      colormap gray; axis image; axis off;
      title(sprintf('Deconvolved Image - IRN-NQP.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
            irn_01L2_snr, irn_01L2_ssim ));
    end

    if(TVBC_CODE) %-- TVBC (Total Variation Boundary Condition)
    
      %noise level: 0.01

      lambda = 1e-3;
      beta   = 0.01;

      tic;
      TVBC_Ig_001L2 = irntvc(Ig_001L2, KC, lambda, beta, 2, 1, 0, Inf, 1e-4, 1e-4, 10, Ig_001L2);
      ttvbc_001L2 = toc;

      tvbc_001L2_snr  = snr(Ig, TVBC_Ig_001L2);
      tvbc_001L2_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(TVBC_Ig_001L2));

      if(SHOW_IMGS)
        figure; imagesc( Normalize(TVBC_Ig_001L2) ); 
        colormap gray; axis image; axis off;
        title(sprintf('Deconvolved Image - TV-BC.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
               tvbc_001L2_snr, tvbc_001L2_ssim));
      end

      % -------------------

      %noise level: 0.1

      lambda = 5e-2;
      beta   = 0.01;

      tic;
      TVBC_Ig_01L2 = irntvc(Ig_01L2, KC, lambda, beta, 2, 1, 0, Inf, 1e-4, 1e-4, 10, Ig_01L2);
      ttvbc_01L2 = toc;

      tvbc_01L2_snr  = snr(Ig, TVBC_Ig_01L2);
      tvbc_01L2_ssim = ssim_index(255*Normalize(Ig), 255*Normalize(TVBC_Ig_01L2));

      if(SHOW_IMGS)
        figure; imagesc( Normalize(TVBC_Ig_01L2) ); 
        colormap gray; axis image; axis off;
        title(sprintf('Deconvolved Image - TV-BC.\n SNR: %4.1fdB.   SSIM: %1.3fdB.', ...
               tvbc_01L2_snr, tvbc_01L2_ssim));
      end

    else  % IF(TVBC_CODE)

      tvbc_001L2_snr   = NaN;
      tvbc_001L2_ssim  = NaN;
      ttvbc_001L2      = NaN

      tvbc_01L2_snr   = NaN;
      tvbc_01L2_ssim  = NaN;
      ttvbc_01L2      = NaN

    end % _END_ IF(TVBC_CODE)


    if(all == 0)

      disp('l2 Deconvolve TV ');
      disp(' ');
      disp('                    SNR (db)          Time (s)             SSIM ');
      disp(' Img    Noise     IRN-NQP  TVBC    IRN-NQP  TVBC     IRN-NQP     TVBC');
      disp(' ');
      disp(sprintf('%s\t 0.01       %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f', ...
        kImg, ...
        irn_001L2_snr, tvbc_001L2_snr, tirn_001L2, ttvbc_001L2, ...
        tvbc_001L2_ssim, tvbc_001L2_ssim));

      disp(sprintf('     \t 0.1        %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f', ...
        irn_01L2_snr, tvbc_01L2_snr, tirn_01L2, ttvbc_01L2, ...
        tvbc_01L2_ssim, tvbc_01L2_ssim));

    else

      if(k==1) 

        str_all = sprintf('l2 Deconvolve TV \n');
        str_all = [ str_all sprintf('\n') ];
        str_all = [ str_all sprintf('                    SNR (db)          Time (s)             SSIM\n') ];
        str_all = [ str_all, sprintf(' Img    Noise     IRN-NQP  TVBC    IRN-NQP  TVBC     IRN-NQP     TVBC\n') ];
        str_all = [ str_all, sprintf('\n') ];

      end

      str_all = [ str_all, sprintf('%s\t 0.01      %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f\n', ...
        Img{k}, ...
        irn_001L2_snr, tvbc_001L2_snr, tirn_001L2, ttvbc_001L2, ...
        tvbc_001L2_ssim, tvbc_001L2_ssim) ];

      str_all = [ str_all sprintf('    \t 0.1       %5.1f   %5.1f    %5.1f   %5.1f    %5.6f   %5.6f\n', ...
        irn_01L2_snr, tvbc_01L2_snr, tirn_01L2, ttvbc_01L2, ...
        tvbc_01L2_ssim, tvbc_01L2_ssim) ];

      if(k==NImgs)
        sprintf(str_all)
      end

end


  end % _END_ if( strcmp(kImg, 'satell.') )



end % _END_ if( strcmp(example,'l2deconv') || strcmp(example,'all') )

%-----------------------------------------------------------------------------

if( strcmp(example,'l1denoise') || strcmp(example,'all') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%




end % _END_ if( strcmp(example,'l1denoise') || strcmp(example,'all') )


%-----------------------------------------------------------------------------


if( strcmp(example,'l2denoise') || strcmp(example,'all') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L2 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%



end % _END_ if( strcmp(example,'l2denoise') || strcmp(example,'all') )

end % _END_ for(k)

%-----------------------------------------------------------------------------


