
% Simulation code for [1], presented in ICIP'09 (http://www.icip2009.org/)
%
% [1] Paul Rodriguez and Brendt Wohlberg, "A Generalized Vector-Valued Total 
%     Variation Algorithm", in Proceedings of IEEE International Conference 
%     on Image Processing (ICIP), (Cairo, Egypt), doi:10.1109/ICIP.2009.5413587 , 
%     pp. 1309--1312, Nov 2009
%
% Legal:
%   irnIcip09.m is based on NUMIPAD (http://numipad.sf.net). NUMIPAD is free
%   software, you can redistribute it and/or modify it under the terms of 
%   the GNU General Public License (version 2).
%
%   The NUMIPAD library is being developed under U.S. Government contract
%   W-7405-ENG-36 for Los Alamos National Laboratory.
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov



clear all; 
%close all;

SHOW_IMGS = false;
%  SHOW_IMGS = true;

%  example = 'icip09'
%  example = 'l1deconv';
%  example = 'l2deconv';
example = 'l1denoise';
%  example = 'l2denoise';


extreme = 0;


nmppath;

BKS_CODE = exist('MainRestoration');
if( (BKS_CODE == 0) && ( strcmp(example, 'l1deconv') || strcmp(example, 'l2deconv') || strcmp(example, 'icip09') ) )
  disp('NOTE:');
  disp('  The function MainRestoration (code for [BKS]) is not in your path...');
  disp(sprintf('  disabling BKS simulations for %s.\n',example));
  disp('  [BKS] L. Bar, A. Brook, N. Sochen and N. Kiryati');
  disp('        "Deblurring of Color Images Corrupted by Impulsive Noise" ');
  disp('        IEEE Transactions on Image Processing, 16 (1101-1111), 2007');
end

BnG_CODE = exist('tvdenoise');
if( (BnG_CODE == 0) && ( strcmp(example, 'l2denoise') || strcmp(example, 'icip09') ) )
  disp('NOTE:');
  disp('  The function tvdenoise (code for [BnG]) is not in your path...');
  disp(sprintf('  disabling BnG simulations for %s.\n',example));
  disp('  [BnG] an implementation of the fast dual minimization of VTV [1].');
  disp('        Code may be downloaded from:');
  disp('        http://www.mathworks.fr/matlabcentral/fileexchange/16236');
  disp('  [1] X. Bresson and T. Chan');
  disp('      "Fast dual minimization of the vectorial total variation norm and');
  disp('       applications to color image processing"');
  disp('      Journal of Inverse Problems and Imaging, 2:4(455--484), 2008');
end


if( exist('color_imgs/peppers_color.png') && exist('color_imgs/mandrill_color.png') ...
    && exist('color_imgs/lena_color_256.png') )
  disp(sprintf('\nRunning %s simulation... \n', example));
else
  disp(' ');
  disp('One or more test images are not in the current directory...');
  disp('You may download them from:');
  disp('http://sites.google.com/a/istec.net/prodrig/Home/en/pubs');
  disp('look for "test images" under "A Generalized Vector-Valued Total ');
  disp('Variation Algorithm"');
  disp(' ');
  disp('Exiting simulation code...');
  return;
end

if( strcmp(example, 'icip09') ) 

  str_all = [];
  ICIP_FLAG = 1;

else

  ICIP_FLAG = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pepImg = imread('peppers_color.png');
pepImg = double(pepImg)/255.0;

lenImg = imread('lena_color_256.png');
lenImg = double(lenImg)/255.0;

mdrilImg = imread('mandrill_color.png');
mdrilImg = double(mdrilImg)/255.0;


%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));


%%%%%%%%%%%%%%%%%%%%%%%%
%       kernels        %
%%%%%%%%%%%%%%%%%%%%%%%%

kernel_BSK = fspecial('disk',3.2);

K_BSK = @(x) imfilter(x, kernel_BSK, 'symmetric','conv');
KT_BSK = @(x) K_BSK(x);
KC_BSK = {K_BSK, KT_BSK};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Blurred & noisy images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lenImgBlur = K_BSK(lenImg);

lenImgBlur_01L1 = imnoise(lenImgBlur, 'salt & pepper', 0.1);
lenImgBlur_03L1 = imnoise(lenImgBlur, 'salt & pepper', 0.3);

% ---

lenImgBlur_01L2 = imnoise(lenImgBlur, 'gaussian', 0, ...
                          0.01*max(lenImgBlur(:)) ); %NOTE sigma^2 in imnoise

lenImgBlur_005L2 = imnoise(lenImgBlur,'gaussian', 0, ...
                           0.0025*max(lenImgBlur(:)) );

lenImgBlur_001L2 = imnoise(lenImgBlur,'gaussian', 0, ...
                           0.0001*max(lenImgBlur(:)) );
% -- BSK example ---:
lenImgBlur_bskL2 = imnoise(lenImgBlur,'gaussian', 0, ...
                           0.00001*max(lenImgBlur(:)) ); 

% ---

mdrilImg_01L1 = imnoise(mdrilImg, 'salt & pepper', 0.1);
mdrilImg_03L1 = imnoise(mdrilImg, 'salt & pepper', 0.3);

mdrilImg_01L2 = imnoise(mdrilImg,  'gaussian', 0, ...
                        0.01*max(mdrilImg(:)) ); %NOTE sigma^2 in imnoise
mdrilImg_005L2 = imnoise(mdrilImg, 'gaussian', 0, ...
                         0.0025*max(mdrilImg(:)) );

% ---

pepImg_01L1 = imnoise(pepImg, 'salt & pepper', 0.1);
pepImg_03L1 = imnoise(pepImg, 'salt & pepper', 0.3);

pepImg_01L2 = imnoise(pepImg, 'gaussian', 0, 0.01*max(pepImg(:)) );
pepImg_005L2 = imnoise(pepImg, 'gaussian', 0, 0.0025*max(pepImg(:)) );


%-----------------------------------------------------------------------------


if( strcmp(example,'l1deconv') || strcmp(example,'icip09') )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Deconvolved        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noise level: 0.1

%-- IRN
lambda = 0.035;
pars = irntvInputPars('l1tv');

pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 4;
%  pars.pcgtol_ini   = 1e-4;
pars.U0           = lenImgBlur_01L1;

tic;
IRN_lenImgBlur_01L1 = irntv(lenImgBlur_01L1, KC_BSK, lambda, pars);
tirn_01L1n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_lenImgBlur_01L1) );
  axis image; axis off;
  title(sprintf('Deconvolved Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, IRN_lenImgBlur_01L1), tirn_01L1n0));
end

if(BKS_CODE)
%-- BSK MSTV (Mumford-Shah TV)

Params = SetParams;
Params.beta=0.5;              % check - same value as bar-2007-deblurring
Params.alpha=0.1;             % check
Params.epsilon=0.1;           % check
Params.gamma = 2*10^(-3);     % as "\mu in bar-2007-deblurring table IV

%  [uh_mstv,V] = MainRestoration(z, kernel, 'L1', 'MSTV', Params);
tic;
[MSTV_lenImgBlur_01L1, V_MSTV_01L1] = MainRestoration(lenImgBlur_01L1, ...
                                            kernel_BSK, 'L1', 'MSTV', Params);
tbsk_01L1n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( MSTV_lenImgBlur_01L1 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, MSTV_lenImgBlur_01L1), tbsk_01L1n0));
end

%-- BSK 1^1-TV 
Params = SetParams;
Params.beta=0.1;              % check - same value as bar-2007-deblurring

tic;
[BKSL1_lenImgBlur_01L1, V_BKSL1_01L1] = MainRestoration(lenImgBlur_01L1, ...
                                              kernel_BSK, 'L1', 'L1', Params);
tbsk_01L1n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( BKSL1_lenImgBlur_01L1 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - BKSL1. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, BKSL1_lenImgBlur_01L1), tbsk_01L1n1));
end

end % _END_ if(BKS_CODE)
%------------------------------------------------------------------------------


%noise level: 0.3

lambda = 0.070;
pars = irntvInputPars('l1tv');

pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 6;
%  pars.pcgtol_ini   = 1e-4;
pars.U0           = lenImgBlur_03L1;


tic;
IRN_lenImgBlur_03L1 = irntv(lenImgBlur_03L1, KC_BSK, lambda, pars);
tirn_03L1n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_lenImgBlur_03L1) );
  axis image; axis off;
  title(sprintf('Deconvolved Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, IRN_lenImgBlur_03L1), tirn_03L1n0));
end


if(BKS_CODE)
%-- BSK MSTV (Mumford-Shah TV)

Params = SetParams;
Params.beta=1.1;            % as describe in bar-2007-deblurring, table IV
Params.alpha=0.5;           % as describe in bar-2007-deblurring, table IV
Params.epsilon=0.1;         % check
Params.gamma = 2*10^(-3);

%  [uh_mstv,V] = MainRestoration(z, kernel, 'L1', 'MSTV', Params);
tic;
[MSTV_lenImgBlur_03L1, V_MSTV_03L1] = MainRestoration(lenImgBlur_03L1, ...
                                            kernel_BSK, 'L1', 'MSTV', Params);
tbsk_03L1n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( MSTV_lenImgBlur_03L1 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, MSTV_lenImgBlur_03L1), tbsk_03L1n0));
end

%-- BSK L1 

Params = SetParams;
Params.beta=0.2;            % as describe in bar-2007-deblurring, table IV

tic;
[BKSL1_lenImgBlur_03L1, V_MSTV_03L1] = MainRestoration(lenImgBlur_03L1, ...
                                              kernel_BSK, 'L1', 'L1', Params);
tbsk_03L1n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( BKSL1_lenImgBlur_03L1 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, BKSL1_lenImgBlur_03L1), tbsk_03L1n1));
end

end % _END_  if(BKS_CODE)

if(ICIP_FLAG == 0)

  disp('L1 Deconvolve Vector TV');  
  disp(' ');
  disp('                          SNR (db)                         Time (s)');
  disp(' Img      Noise    VTV-IRN    MSTV    BKSL1          VTV IRN    MSTV    BKSL1');
  disp(' ');
  if(BKS_CODE)
  disp(sprintf('Lena       0.1      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_01L1), snr(lenImg, MSTV_lenImgBlur_01L1), snr(lenImg, BKSL1_lenImgBlur_01L1), ...
     tirn_01L1n0, tbsk_01L1n0, tbsk_01L1n1));
  disp(sprintf('           0.3      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_03L1), snr(lenImg, MSTV_lenImgBlur_03L1), snr(lenImg, BKSL1_lenImgBlur_03L1), ...
     tirn_03L1n0, tbsk_03L1n0, tbsk_03L1n1));
  else
  disp(sprintf('Lena       0.1      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_01L1), NaN, NaN, ...
     tirn_01L1n0, NaN, NaN));
  disp(sprintf('           0.3      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_03L1), NaN, NaN, ...
     tirn_03L1n0, NaN, NaN));

  end % _END_ if(BKS_CODE)

else % IF(icip)

  str = sprintf('L1 Deconvolve Vector TV\n\n');  
  str_all = [str_all str];
  str = sprintf('                          SNR (db)                         Time (s)\n');
  str_all = [str_all str];
  str = sprintf(' Img      Noise    VTV-IRN    MSTV    BKSL1          VTV IRN    MSTV    BKSL1\n\n');
  str_all = [str_all str];

  if(BKS_CODE)
  str = sprintf('Lena       0.1      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f\n', ...
     snr(lenImg, IRN_lenImgBlur_01L1), snr(lenImg, MSTV_lenImgBlur_01L1), snr(lenImg, BKSL1_lenImgBlur_01L1), ...
     tirn_01L1n0, tbsk_01L1n0, tbsk_01L1n1);
  str_all = [str_all str];
  str = sprintf('           0.3      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f\n\n\n', ...
     snr(lenImg, IRN_lenImgBlur_03L1), snr(lenImg, MSTV_lenImgBlur_03L1), snr(lenImg, BKSL1_lenImgBlur_03L1), ...
     tirn_03L1n0, tbsk_03L1n0, tbsk_03L1n1);
  str_all = [str_all str];
  else
  str = sprintf('Lena       0.1      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f\n', ...
     snr(lenImg, IRN_lenImgBlur_01L1), NaN, NaN, ...
     tirn_01L1n0, NaN, NaN);
  str_all = [str_all str];
  str = sprintf('           0.3      %5.1f     %5.1f    %5.1f          %5.1f     %5.1f    %5.1f\n\n\n', ...
     snr(lenImg, IRN_lenImgBlur_03L1), NaN, NaN, ...
     tirn_03L1n0, NaN, NaN);
  str_all = [str_all str];

  end % _END_ if(BKS_CODE)

end

end % _END_ if( strcmp(example,'l1deconv') || strcmp(example,'icip09') )

%-----------------------------------------------------------------------------


if( strcmp(example,'l2deconv') || strcmp(example,'icip09') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L2 Deconvolved        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(extreme == 1)

%noise level: 0.1

%-- IRN
lambda  = 0.04;
pars = irntvInputPars('l2tv');

pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 3;
pars.pcgtol_ini   = 1e-4;
%  pars.U0           = lenImgBlur_01L2;


tic;
IRN_lenImgBlur_01L2 = irntv(lenImgBlur_01L2, KC_BSK, lambda, pars);
tirn_01L2n0 = toc;


if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_lenImgBlur_01L2) );
  axis image; axis off;
  title(sprintf('Deconvolved Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, IRN_lenImgBlur_01L2), tirn_01L2n0));
end

if(BKS_CODE)
%-- BSK MSTV (Mumford-Shah TV)


Params      = SetParams;
Params.beta = 0.01*max(lenImgBlur(:));

tic;
[MSTV_lenImgBlur_01L2, V_MSTV_001L2] = MainRestoration(lenImgBlur_01L2, ...
                                              kernel_BSK, 'L2', 'MS', Params);
tbsk_01L2n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( MSTV_lenImgBlur_01L2 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, MSTV_lenImgBlur_01L2), tbsk_01L2n0));
end

end % _END_ if(BKS_CODE)
end % _END_ if(extreme == 1)
%-----------------------------------------------------------------------------

%noise level: 0.05

%-- IRN
lambda  = 0.01;
pars = irntvInputPars('l2tv');

pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 3;
%  pars.pcgtol_ini   = 1e-4;
pars.U0           = lenImgBlur_005L2;


tic;
IRN_lenImgBlur_005L2 = irntv(lenImgBlur_005L2, KC_BSK, lambda, pars);
tirn_005L2n0 = toc;


if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_lenImgBlur_005L2) );
  axis image; axis off;
  title(sprintf('Deconvolved Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, IRN_lenImgBlur_005L2), tirn_005L2n0));
end

if(BKS_CODE)
%-- BSK

Params      = SetParams;
Params.beta = 0.0025*max(lenImgBlur(:));

tic;
[MSTV_lenImgBlur_005L2, V_MSTV_05L2] = MainRestoration(lenImgBlur_005L2, ...
                                              kernel_BSK, 'L2', 'MS', Params);
tbsk_005L2n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( MSTV_lenImgBlur_005L2 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, MSTV_lenImgBlur_005L2), tbsk_005L2n0));
end

%-- BSK l2-VTV

Params      = SetParams;
Params.beta = 0.0025*max(lenImgBlur(:));

tic;
[BSKL2_lenImgBlur_005L2, V_MSTV_05L2] = MainRestoration(lenImgBlur_005L2, ...
                                              kernel_BSK, 'L2', 'L1', Params);
tbsk_005L2n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( BSKL2_lenImgBlur_005L2 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, BSKL2_lenImgBlur_005L2), tbsk_005L2n1));
end

end % _END_ if(BKS_CODE)

%-----------------------------------------------------------------------------


%noise level: 0.01

%-- IRN
lambda  = 0.0005;

pars = irntvInputPars('l2tv');

pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 3;
%  pars.pcgtol_ini   = 1e-4;
pars.U0           = lenImgBlur_001L2;

tic;
IRN_lenImgBlur_001L2 = irntv(lenImgBlur_001L2, KC_BSK, lambda, pars);
tirn_001L2n0 = toc;


if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_lenImgBlur_001L2) );
  axis image; axis off;
  title(sprintf('Deconvolved Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, IRN_lenImgBlur_001L2), tirn_001L2n0));
end

if(BKS_CODE)
%-- BSK MSTV (Mumford-Shah TV)


Params      = SetParams;
Params.beta = 0.0001;

tic;
[MSTV_lenImgBlur_001L2, V_MSTV_001L2] = ...
  MainRestoration(Normalize(lenImgBlur_001L2), kernel_BSK, 'L2', 'MS', Params);
tbsk_001L2n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( MSTV_lenImgBlur_001L2 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, MSTV_lenImgBlur_001L2), tbsk_001L2n0));
end

%-- BSK L2-VTV


Params      = SetParams;
Params.beta = 0.0001;

tic;
[BSKL2_lenImgBlur_001L2, V_MSTV_001L2] = ...
  MainRestoration(Normalize(lenImgBlur_001L2), kernel_BSK, 'L2', 'L1', Params);
tbsk_001L2n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( BSKL2_lenImgBlur_001L2 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, BSKL2_lenImgBlur_001L2), tbsk_001L2n1));
end

end % _END_ if(BKS_CODE)

%-----------------------------------------------------------------------------

%noise level: sqrt(1e-5)   BSK example

%-- IRN
lambda  = 0.0001;
pars = irntvInputPars('l2tv');

pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 3;
%  pars.pcgtol_ini   = 1e-4;
pars.U0           = lenImgBlur_bskL2;

tic;
IRN_lenImgBlur_bskL2 = irntv(lenImgBlur_bskL2, KC_BSK, lambda, pars);
tirn_bskL2n0 = toc;


if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_lenImgBlur_bskL2) );
  axis image; axis off;
  title(sprintf('Deconvolved Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, IRN_lenImgBlur_bskL2), tirn_bskL2n0));
end

if(BKS_CODE)

%-- BSK MSTV (Mumford-Shah TV)


Params      = SetParams;
Params.beta = 0.00001;

tic;
[MSTV_lenImgBlur_bskL2, V_MSTV_bskL2] = ...
  MainRestoration(Normalize(lenImgBlur_bskL2), kernel_BSK, 'L2', 'MS', Params);
tbsk_bskL2n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( MSTV_lenImgBlur_bskL2 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, MSTV_lenImgBlur_bskL2), tbsk_bskL2n0));
end

%-- BSK l2-VTV


Params      = SetParams;
Params.beta = 0.00001;

tic;
[BSKL2_lenImgBlur_bskL2, V_MSTV_bskL2] = ...
  MainRestoration(Normalize(lenImgBlur_bskL2), kernel_BSK, 'L2', 'L1', Params);
tbsk_bskL2n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( BSKL2_lenImgBlur_bskL2 );
  axis image; axis off;
  title(sprintf('Deconvolved Image - MSTV. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(lenImg, BSKL2_lenImgBlur_bskL2), tbsk_bskL2n1));
end

end % _END_ if(BKS_CODE)

%-----------------------------------------------------------------------------

if(ICIP_FLAG == 0)

  disp('L2 Deconvolve Vector TV');
  disp(' ');
  disp('                        SNR (db)                    Time (s)');
  disp(' Img      Noise    VTV-IRN  L2-MS(BSK)   BSK-L2          VTV IRN   L2-MS(BKS)  BSK-L2');
  disp(' ');
  if(BKS_CODE)
  disp(sprintf('Lenna      0.05      %5.1f     %5.1f     %5.1f           %5.1f      %5.1f    %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_005L2), snr(lenImg, MSTV_lenImgBlur_005L2), snr(lenImg, BSKL2_lenImgBlur_005L2), ...
     tirn_005L2n0, tbsk_005L2n0, tbsk_005L2n1));
  disp(sprintf('Lenna      0.01      %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_001L2), snr(lenImg, MSTV_lenImgBlur_001L2), snr(lenImg, BSKL2_lenImgBlur_001L2), ...
     tirn_001L2n0, tbsk_001L2n0, tbsk_001L2n1));
  disp(sprintf('Lenna   sqrt(1e-5)   %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_bskL2), snr(lenImg, MSTV_lenImgBlur_bskL2), snr(lenImg, BSKL2_lenImgBlur_bskL2), ...
     tirn_bskL2n0, tbsk_bskL2n0, tbsk_bskL2n1));
  else
  disp(sprintf('Lenna      0.05      %5.1f     %5.1f     %5.1f           %5.1f      %5.1f    %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_005L2), NaN, NaN, ...
     tirn_005L2n0, NaN, NaN));
  disp(sprintf('Lenna      0.01      %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_001L2), NaN, NaN, ...
     tirn_001L2n0, NaN, NaN));
  disp(sprintf('Lenna   sqrt(1e-5)   %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f', ...
     snr(lenImg, IRN_lenImgBlur_bskL2), NaN, NaN, ...
     tirn_bskL2n0, NaN, NaN));
  end % _END_ if(BKS_CODE)


else

  str = sprintf('L2 Deconvolve Vector TV\n\n');
  str_all = [str_all str];
  str = sprintf('                        SNR (db)                    Time (s)\n');
  str_all = [str_all str];
  str = sprintf(' Img      Noise    VTV-IRN  L2-MS(BSK)   BSK-L2          VTV IRN   L2-MS(BKS)  BSK-L2\n\n');
  str_all = [str_all str];

  if(BKS_CODE)
  str = sprintf('Lena       0.05      %5.1f     %5.1f     %5.1f           %5.1f      %5.1f    %5.1f\n', ...
     snr(lenImg, IRN_lenImgBlur_005L2), snr(lenImg, MSTV_lenImgBlur_005L2), snr(lenImg, BSKL2_lenImgBlur_005L2), ...
     tirn_005L2n0, tbsk_005L2n0, tbsk_005L2n1);
  str_all = [str_all str];
  str = sprintf('Lena       0.01      %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f\n', ...
     snr(lenImg, IRN_lenImgBlur_001L2), snr(lenImg, MSTV_lenImgBlur_001L2), snr(lenImg, BSKL2_lenImgBlur_001L2), ...
     tirn_001L2n0, tbsk_001L2n0, tbsk_001L2n1);
  str_all = [str_all str];
  str = sprintf('Lena    sqrt(1e-5)   %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f\n\n\n', ...
     snr(lenImg, IRN_lenImgBlur_bskL2), snr(lenImg, MSTV_lenImgBlur_bskL2), snr(lenImg, BSKL2_lenImgBlur_bskL2), ...
     tirn_bskL2n0, tbsk_bskL2n0, tbsk_bskL2n1);
  str_all = [str_all str];
  else
  str = sprintf('Lena       0.05      %5.1f     %5.1f     %5.1f           %5.1f      %5.1f    %5.1f\n', ...
     snr(lenImg, IRN_lenImgBlur_005L2), NaN, NaN, ...
     tirn_005L2n0, NaN, NaN);
  str_all = [str_all str];
  str = sprintf('Lena       0.01      %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f\n', ...
     snr(lenImg, IRN_lenImgBlur_001L2), NaN, NaN, ...
     tirn_001L2n0, NaN, NaN);
  str_all = [str_all str];
  str = sprintf('Lena    sqrt(1e-5)   %5.1f     %5.1f      %5.1f          %5.1f      %5.1f     %5.1f\n\n\n', ...
     snr(lenImg, IRN_lenImgBlur_bskL2), NaN, NaN, ...
     tirn_bskL2n0, NaN, NaN);
  str_all = [str_all str];
  end % _END_ if(BKS_CODE)

end

end % _END_ if( strcmp(example,'l2deconv') || strcmp(example,'icp09') )

%-----------------------------------------------------------------------------




if( strcmp(example,'l1denoise') || strcmp(example,'icip09') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noise level: 0.1

%-- IRN
lambda  = 1.1;

pars = irntvInputPars('l1tv');

pars.pcgtol_ini = 1e-4;
pars.epsF       = 1e-2;    
pars.epsR       = 1e-4;
pars.loops      = 2;


%-- peppers
tic;
IRN_pepImg_01L1 = irntv(pepImg_01L1, [], lambda, pars);
tirn_01L1n0 = toc;


if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_pepImg_01L1) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(pepImg, IRN_pepImg_01L1), tirn_01L1n0));
end

%-- mandrill
tic;
IRN_mdrilImg_01L1 = irntv(mdrilImg_01L1, [], lambda, pars);
tirn_01L1n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_mdrilImg_01L1) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(mdrilImg, IRN_mdrilImg_01L1), tirn_01L1n1));
end

%-----------------------------------------------------------------------------

%noise level: 0.3

%-- IRN
lambda  = 1.2;

pars = irntvInputPars('l1tv');

pars.pcgtol_ini = 1e-4;
pars.epsF       = 1e-2;    
pars.epsR       = 1e-4;
pars.loops      = 2;


%-- peppers
tic;
IRN_pepImg_03L1 = irntv(pepImg_03L1, [], lambda, pars);
tirn_03L1n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_pepImg_03L1) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(pepImg, IRN_pepImg_03L1), tirn_03L1n0));
end

%-- mandrill

tic;
IRN_mdrilImg_03L1 = irntv(mdrilImg_03L1, [], lambda, pars);
tirn_03L1n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_mdrilImg_03L1) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(mdrilImg, IRN_mdrilImg_03L1), tirn_03L1n1));
end

%-----------------------------------------------------------------------------

if(ICIP_FLAG == 0)

  disp('L1 Denoise Vector TV (Vector IRN algorithm)');
  disp(' ');
  disp('                   SNR (db)          Time (s)');
  disp(' Img      Noise    VTV IRN           VTV IRN');
  disp(' ');
  disp(sprintf('Peppers    0.1       %5.1f            %5.1f  ', ...
      snr(pepImg, IRN_pepImg_01L1), tirn_01L1n0));
  disp(sprintf('           0.3       %5.1f            %5.1f  ', ...
      snr(pepImg, IRN_pepImg_03L1), tirn_03L1n0));
  disp(sprintf('Mandrill   0.1       %5.1f            %5.1f  ', ...
      snr(mdrilImg, IRN_mdrilImg_01L1), tirn_01L1n1));
  disp(sprintf('           0.3       %5.1f            %5.1f  ', ...
      snr(mdrilImg, IRN_mdrilImg_03L1), tirn_03L1n1));

else

  str = sprintf('L1 Denoise Vector TV (Vector IRN algorithm)\n\n');
  str_all = [str_all str];
  str = sprintf('                   SNR (db)          Time (s)\n');
  str_all = [str_all str];
  str = sprintf(' Img      Noise    VTV IRN           VTV IRN\n\n');
  str_all = [str_all str];
  
  str = sprintf('Peppers    0.1       %5.1f            %5.1f  \n', ...
      snr(pepImg, IRN_pepImg_01L1), tirn_01L1n0);
  str_all = [str_all str];
  str = sprintf('           0.3       %5.1f            %5.1f  \n', ...
      snr(pepImg, IRN_pepImg_03L1), tirn_03L1n0);
  str_all = [str_all str];
  str = sprintf('Mandrill   0.1       %5.1f            %5.1f  \n', ...
      snr(mdrilImg, IRN_mdrilImg_01L1), tirn_01L1n1);
  str_all = [str_all str];
  str = sprintf('           0.3       %5.1f            %5.1f  \n\n\n', ...
      snr(mdrilImg, IRN_mdrilImg_03L1), tirn_03L1n1);
  str_all = [str_all str];

end

end % _END_ if( strcmp(example,'l1denoise') || strcmp(example,'icip09') )


%-----------------------------------------------------------------------------


if( strcmp(example,'l2denoise') || strcmp(example,'icip09') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L2 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%noise level: 0.1

%-- IRN
lambda  = 0.3;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini = 1e-4;
%  pars.epsF       = 1e-1;    
%  pars.epsR       = 1e-2;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops      = 1;


tic;
IRN_pepImg_01L2 = irntv(pepImg_01L2, [], lambda, pars);
tirn_01L2n0 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_pepImg_01L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(pepImg, IRN_pepImg_01L2), tirn_01L2n0));
end

if(BnG_CODE)
%% http://www.mathworks.fr/matlabcentral/fileexchange/16236
%%  Pascal Getreuer based on X. Bresson and T.F. Chan, 
%%  "Fast Minimization of the Vectorial Total Variation Norm and Applications 
%%  to Color Image Processing", CAM Report 07-25.

lambda  = 0.15;

tic;
BnG_pepImg_01L2 = tvdenoise(pepImg_01L2, 1/lambda);
tirn_01L2n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(BnG_pepImg_01L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - BnG. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(pepImg, BnG_pepImg_01L2), tirn_01L2n1));
end

end % _END_ if(BnG_CODE)

%------------

%-- IRN
lambda  = 0.3;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini = 1e-4;
%  pars.epsF       = 1e-1;    
%  pars.epsR       = 1e-2;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops      = 1;

tic;
IRN_mdrilImg_01L2 = irntv(mdrilImg_01L2, [], lambda, pars);
tirn_01L2n2 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_mdrilImg_01L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(mdrilImg, IRN_mdrilImg_01L2), tirn_01L2n2));
end


if(BnG_CODE)
%% http://www.mathworks.fr/matlabcentral/fileexchange/16236
%%  Pascal Getreuer based on X. Bresson and T.F. Chan, 
%%  "Fast Minimization of the Vectorial Total Variation Norm and Applications 
%%  to Color Image Processing", CAM Report 07-25.

lambda  = 0.15;

tic;
BnG_mdrilImg_01L2 = tvdenoise(mdrilImg_01L2, 1/lambda);
tirn_01L2n3 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(BnG_mdrilImg_01L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - BnG. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(mdrilImg, BnG_mdrilImg_01L2), tirn_01L2n3));
end

end % _END_ if(BnG_CODE)

%-----------------------------------------------------------------------------


%noise level: 0.05

%-- IRN
lambda  = 0.1;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini = 1e-4;
%  pars.epsF       = 1e-1;    
%  pars.epsR       = 1e-2;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops      = 1;


tic;
IRN_pepImg_005L2 = irntv(pepImg_005L2, [], lambda, pars);
tirn_005L2n0 = toc;


if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_pepImg_005L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(pepImg, IRN_pepImg_005L2), tirn_005L2n0));
end

if(BnG_CODE)
%% http://www.mathworks.fr/matlabcentral/fileexchange/16236
%%  Pascal Getreuer based on X. Bresson and T.F. Chan, 
%%  "Fast Minimization of the Vectorial Total Variation Norm and Applications 
%%  to Color Image Processing", CAM Report 07-25.

lambda  = 0.05;

tic;
BnG_pepImg_005L2 = tvdenoise(pepImg_005L2, 1/lambda);
tirn_005L2n1 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(BnG_pepImg_005L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - BnG. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(pepImg, BnG_pepImg_005L2), tirn_005L2n1));
end


end % _END_ if(BnG_CODE)

%------------

%-- IRN
lambda  = 0.1;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini = 1e-4;
%  pars.epsF       = 1e-1;    
%  pars.epsR       = 1e-2;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops      = 1;


tic;
IRN_mdrilImg_005L2 = irntv(mdrilImg_005L2, [], lambda, pars);
tirn_005L2n2 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(IRN_mdrilImg_005L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - Vector IRN. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(mdrilImg, IRN_mdrilImg_005L2), tirn_005L2n2));
end

if(BnG_CODE)
%% http://www.mathworks.fr/matlabcentral/fileexchange/16236
%%  Pascal Getreuer based on X. Bresson and T.F. Chan, 
%%  "Fast Minimization of the Vectorial Total Variation Norm and Applications 
%%  to Color Image Processing", CAM Report 07-25.

lambda  = 0.05;

tic;
BnG_mdrilImg_005L2 = tvdenoise(mdrilImg_005L2, 1/lambda);
tirn_005L2n3 = toc;

if(SHOW_IMGS)
  figure; imagesc( Normalize(BnG_mdrilImg_005L2) );
  axis image; axis off;
  title(sprintf('Denoised Image - BnG. SNR: %4.1fdB.\n Time %4.1f sec', ...
               snr(mdrilImg, BnG_mdrilImg_005L2), tirn_005L2n3));
end


end % _END_ if(BnG_CODE)

%-----------------------------------------------------------------------------

if(ICIP_FLAG == 0)

  disp('L2 Denoise Vector TV (Vector IRN algorithm)');
  disp(' ');
  disp('                         SNR (db)                    Time (s)');
  disp(' Img      Noise   VTV-IRN   bresson-2008-fast     VTV IRN   bresson-2008-fast');
  disp(' ');
  if(BnG_CODE)
  disp(sprintf('Peppers,   0.1      %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(pepImg, IRN_pepImg_01L2), snr(pepImg, BnG_pepImg_01L2), ...
      tirn_01L2n0, tirn_01L2n1));
  disp(sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(pepImg, IRN_pepImg_005L2), snr(pepImg, BnG_pepImg_005L2), ...
      tirn_005L2n0, tirn_005L2n1));
  disp(sprintf('Mandrill,  0.1      %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(mdrilImg, IRN_mdrilImg_01L2), snr(mdrilImg, BnG_mdrilImg_01L2), ...
      tirn_01L2n2, tirn_01L2n3));
  disp(sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(mdrilImg, IRN_mdrilImg_005L2), snr(mdrilImg, BnG_mdrilImg_005L2), ...
      tirn_005L2n2, tirn_005L2n3));
  else
  disp(sprintf('Peppers,   0.1      %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(pepImg, IRN_pepImg_01L2), NaN, ...
      tirn_01L2n0, NaN));
  disp(sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(pepImg, IRN_pepImg_005L2), NaN, ...
      tirn_005L2n0, NaN));
  disp(sprintf('Mandrill,  0.1      %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(mdrilImg, IRN_mdrilImg_01L2), NaN, ...
      tirn_01L2n2, NaN));
  disp(sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f', ...
      snr(mdrilImg, IRN_mdrilImg_005L2), NaN, ...
      tirn_005L2n2, NaN));
  end % _END_ if(BnG_CODE)

else

  str = sprintf('L2 Denoise Vector TV (Vector IRN algorithm)\n\n');
  str_all = [str_all str];
  str = sprintf('                         SNR (db)                    Time (s)\n');
  str_all = [str_all str];
  str = sprintf(' Img      Noise   VTV-IRN   bresson-2008-fast     VTV IRN   bresson-2008-fast\n\n');
  str_all = [str_all str];

  if(BnG_CODE)
  str = sprintf('Peppers,   0.1      %5.1f      %5.1f                %5.1f         %5.1f\n', ...
      snr(pepImg, IRN_pepImg_01L2), snr(pepImg, BnG_pepImg_01L2), ...
      tirn_01L2n0, tirn_01L2n1);
  str_all = [str_all str];
  str = sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f\n', ...
      snr(pepImg, IRN_pepImg_005L2), snr(pepImg, BnG_pepImg_005L2), ...
      tirn_005L2n0, tirn_005L2n1);
  str_all = [str_all str];
  str = sprintf('Mandrill,  0.1      %5.1f      %5.1f                %5.1f         %5.1f\n', ...
      snr(mdrilImg, IRN_mdrilImg_01L2), snr(mdrilImg, BnG_mdrilImg_01L2), ...
      tirn_01L2n2, tirn_01L2n3);
  str_all = [str_all str];
  str = sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f\n\n\n', ...
      snr(mdrilImg, IRN_mdrilImg_005L2), snr(mdrilImg, BnG_mdrilImg_005L2), ...
      tirn_005L2n2, tirn_005L2n3);
  str_all = [str_all str];
  else
  str = sprintf('Peppers,   0.1      %5.1f      %5.1f                %5.1f         %5.1f\n', ...
      snr(pepImg, IRN_pepImg_01L2), NaN, ...
      tirn_01L2n0, NaN);
  str_all = [str_all str];
  str = sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f\n', ...
      snr(pepImg, IRN_pepImg_005L2), NaN, ...
      tirn_005L2n0, NaN);
  str_all = [str_all str];
  str = sprintf('Mandrill,  0.1      %5.1f      %5.1f                %5.1f         %5.1f\n', ...
      snr(mdrilImg, IRN_mdrilImg_01L2), NaN, ...
      tirn_01L2n2, NaN);
  str_all = [str_all str];
  str = sprintf('           0.05     %5.1f      %5.1f                %5.1f         %5.1f\n\n\n', ...
      snr(mdrilImg, IRN_mdrilImg_005L2), NaN, ...
      tirn_005L2n2, NaN);
  str_all = [str_all str];
  end % _END_ if(BnG_CODE)


end

end % _END_ if( strcmp(example,'l2denoise') || strcmp(example,'icp09') )

%-----------------------------------------------------------------------------

if(ICIP_FLAG == 1)

  disp(str_all);
end
