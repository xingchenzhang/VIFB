
% Test code for irntv.m
%
% Legal:
%   irnTest.m is part of NUMIPAD (http://numipad.sf.net). 
%
%   The NUMIPAD library is being developed under U.S. Government contract
%   W-7405-ENG-36 for Los Alamos National Laboratory.
%  
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov



%  example = 'l1deconv';
%  example = 'l2deconv';
%  example = 'l1denoise';
example = 'l2denoise';

%  Img = 'lena';
%  Img = 'peppers';
%  Img = 'goldhill';
Img = 'lena';
%  nmppath;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(Img)

  case{'lena'}
    Ig = double( imread('lena_gray_512.png') ) / 255;
    Ic = double( imread('lena_color_512.png') ) / 255;

  case{'peppers'}
    Ig = double( imread('peppers_gray.png') ) / 255;
    Ic = double( imread('peppers_color.png') ) / 255;

  case{'goldhill'}
    Ig = double( imread('goldhill_gray.png') ) / 255;
    Ic = double( imread('goldhill_color.png') ) / 255;

  case{'none'}
    disp(' ');
    disp('Select a test image (see code for an example)...');
    disp(' ');
    disp('Exiting irnTest code...');
    return;

  otherwise
    error('Not a valid image\n');

end



%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));

%%%%%%%%%%%%%%%%%%%%%%%%
%       kernels        %
%%%%%%%%%%%%%%%%%%%%%%%%

kernel = fspecial('disk',3.2);

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

    IcBlur = K(Ic);
    Ic_01L1 = imnoise(IcBlur, 'salt & pepper', 0.1);
    Ic_03L1 = imnoise(IcBlur, 'salt & pepper', 0.3);

  case{'l1denoise'}

    Ig_01L1 = imnoise(Ig, 'salt & pepper', 0.1);
    Ig_03L1 = imnoise(Ig, 'salt & pepper', 0.3);

    Ic_01L1 = imnoise(Ic, 'salt & pepper', 0.1);
    Ic_03L1 = imnoise(Ic, 'salt & pepper', 0.3);

  case{'l2deconv'}

    IgBlur = K(Ig);
                                             %NOTE sigma^2 in imnoise
    Ig_01L2 = imnoise(IgBlur, 'gaussian', 0, 0.01*max(IgBlur(:)) ); 
    Ig_001L2 = imnoise(IgBlur,'gaussian', 0, 0.0001*max(IgBlur(:)) );

    IcBlur = K(Ic);
    Ic_01L2 = imnoise(IcBlur,'gaussian', 0, 0.01*max(IgBlur(:)) );
    Ic_001L2 = imnoise(IcBlur,'gaussian', 0, 0.0001*max(IgBlur(:)) ); 

  case{'l2denoise'}
    Ig_01L2 = imnoise(Ig,'gaussian', 0, 0.01*max(Ig(:)) );
    Ic_01L2 = imnoise(Ic,'gaussian', 0, 0.01*max(Ig(:)) ); 

    Ig_005L2 = imnoise(Ig,'gaussian', 0, 0.0025*max(Ig(:)) ); 
    Ic_005L2 = imnoise(Ic,'gaussian', 0, 0.0025*max(Ig(:)) ); 

end



%-----------------------------------------------------------------------------


if( strcmp(example,'l1deconv') || strcmp(example,'all') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Deconvolved        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noise level: 0.1

lambda = 0.45;

pars = irntvInputPars('l1tv');

%  pars.pcgtol_ini   = 1e-4;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 4;

%-----------------%
% -- Grayscale -- %

pars.U0           = Ig_01L1;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ig_01L1 = irntv(Ig_01L1, KC, lambda, pars);

figure; imagesc( Normalize(IRN_Ig_01L1) ); 
colormap gray; axis image; axis off;
title(sprintf('Deconvolved Image - Scalar IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ig, IRN_Ig_01L1), snr(Ig, Ig_01L1)));

%-------------%
% -- Color -- %

pars.U0           = Ic_01L1;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ic_01L1 = irntv(Ic_01L1, KC, lambda, pars);

figure; imagesc( Normalize(IRN_Ic_01L1) ); axis image; axis off;
title(sprintf('Deconvolved Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ic, IRN_Ic_01L1), snr(Ic, Ic_01L1)));


%------------------
%------------------

%noise level: 0.3

lambda = 0.95;

pars = irntvInputPars('l1tv');

%  pars.pcgtol_ini   = 1e-4;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 4;

%-----------------%
% -- Grayscale -- %

pars.U0           = Ig_03L1;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ig_03L1 = irntv(Ig_03L1, KC, lambda, pars);

figure; imagesc( Normalize(IRN_Ig_03L1) ); 
colormap gray; axis image; axis off;
title(sprintf('Deconvolved Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n', ...
               snr(Ig, IRN_Ig_03L1), snr(Ig, Ig_03L1)));


%-------------%
% -- Color -- %

pars.U0           = Ic_03L1;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ic_03L1 = irntv(Ic_03L1, KC, lambda, pars);

figure; imagesc( Normalize(IRN_Ic_03L1) ); axis image; axis off;
title(sprintf('Deconvolved Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB). \n', ...
               snr(Ic, IRN_Ic_03L1), snr(Ic, Ic_03L1)));


end % _END_ if( strcmp(example,'l1deconv') || strcmp(example,'all') )

%-----------------------------------------------------------------------------


if( strcmp(example,'l2deconv') || strcmp(example,'all') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L2 Deconvolved        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sigma (noise level) : 0.01

lambda  = 0.001;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini   = 1e-4;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 3;

%-----------------%
% -- Grayscale -- %

pars.U0           = Ig_001L2;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ig_001L2 = irntv(Ig_001L2, KC, lambda, pars);


figure; imagesc( Normalize(IRN_Ig_001L2) ); 
colormap gray; axis image; axis off;
title(sprintf('Deconvolved Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ig, IRN_Ig_001L2), snr(Ig, Ig_001L2)));


%-------------%
% -- Color -- %

pars.U0           = Ic_001L2;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ic_001L2 = irntv(Ic_001L2, KC, lambda, pars);


figure; imagesc( Normalize(IRN_Ic_001L2) ); axis image; axis off;
title(sprintf('Deconvolved Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ic, IRN_Ic_001L2), snr(Ic, Ic_001L2)));


% sigma (noise level) : 0.1

lambda  = 0.05;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini   = 1e-4;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.05;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.01;
pars.loops        = 3;


%-----------------%
% -- Grayscale -- %

pars.U0           = Ig_01L2;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ig_01L2 = irntv(Ig_01L2, KC, lambda, pars);


figure; imagesc( Normalize(IRN_Ig_01L2) ); 
colormap gray; axis image; axis off;
title(sprintf('Deconvolved Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ig, IRN_Ig_01L2), snr(Ig, Ig_01L2)));


%-------------%
% -- Color -- %

pars.U0           = Ic_01L2;  % initial solution

% >> Deconvolution via IRN <<
IRN_Ic_01L2 = irntv(Ic_01L2, KC, lambda, pars);


figure; imagesc( Normalize(IRN_Ic_01L2) ); 
colormap gray; axis image; axis off;
title(sprintf('Deconvolved Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ic, IRN_Ic_01L2), snr(Ic, Ic_01L2)));


end % _END_ if( strcmp(example,'l2deconv') || strcmp(example,'all') )

%-----------------------------------------------------------------------------




if( strcmp(example,'l1denoise') || strcmp(example,'all') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%noise level: 10%

%-- IRN
lambda  = 1.1;

pars = irntvInputPars('l1tv');

pars.pcgtol_ini = 1e-4;
pars.epsf       = 1e-2;    
pars.epsr       = 1e-4;
pars.loops      = 2;

%-----------------%
% -- Grayscale -- %

% >> Denoising via IRN <<
IRN_Ig_01L1 = irntv(Ig_01L1, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ig_01L1) ); 
colormap gray; axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ig, IRN_Ig_01L1), snr(Ig, Ig_01L1)));


%-------------%
% -- Color -- %

% >> Denoising via IRN <<
IRN_Ic_01L1 = irntv(Ic_01L1, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ic_01L1) ); 
colormap gray; axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ic, IRN_Ic_01L1), snr(Ic, Ic_01L1)));


%-----------------
%-----------------

%noise level: 30%

lambda  = 1.2;

pars = irntvInputPars('l1tv');

pars.pcgtol_ini   = 1e-4;
pars.adapt_epsR   = 1;
pars.epsR_cutoff  = 0.01;
pars.adapt_epsF   = 1;
pars.epsF_cutoff  = 0.05;
pars.loops        = 3;


%-----------------%
% -- Grayscale -- %

% >> Denoising via IRN <<
IRN_Ig_03L1 = irntv(Ig_03L1, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ig_03L1) ); 
colormap gray; axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ig, IRN_Ig_03L1), snr(Ig, Ig_03L1)));


%-------------%
% -- Color -- %

% >> Denoising via IRN <<
IRN_Ic_03L1 = irntv(Ic_03L1, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ic_03L1) ); axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ic, IRN_Ic_03L1), snr(Ic, Ic_03L1)));


end % _END_ if( strcmp(example,'l1denoise') || strcmp(example,'all') )


%-----------------------------------------------------------------------------


if( strcmp(example,'l2denoise') || strcmp(example,'all') )

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L2 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sigma (noise level) : 0.05

lambda  = 0.05;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini = 1e-4;
pars.epsf       = 1e-2;    
pars.epsr       = 1e-5;    
pars.loops      = 2;

%-----------------%
% -- Grayscale -- %

% >> Denoising via IRN <<
IRN_Ig_005L2 = irntv(Ig_005L2, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ig_005L2) ); 
colormap gray; axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ig, IRN_Ig_005L2), snr(Ig, Ig_005L2)));

%-------------%
% -- Color -- %

% >> Denoising via IRN <<
IRN_Ic_005L2 = irntv(Ic_005L2, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ic_005L2) ); axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ic, IRN_Ic_005L2), snr(Ic, Ic_005L2)));

%-----------------
%-----------------


% sigma (noise level) : 0.1

lambda  = 0.1;

pars = irntvInputPars('l2tv');

pars.pcgtol_ini = 1e-4;
pars.epsf       = 1e-2;    
pars.epsr       = 1e-5;    
pars.loops      = 2;

%-----------------%
% -- Grayscale -- %

% >> Denoising via IRN <<
IRN_Ig_01L2 = irntv(Ig_01L2, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ig_01L2) ); 
colormap gray; axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ig, IRN_Ig_01L2), snr(Ig, Ig_01L2)));

%-------------%
% -- Color -- %

% >> Denoising via IRN <<
IRN_Ic_01L2 = irntv(Ic_01L2, [], lambda, pars);


figure; imagesc( Normalize(IRN_Ic_01L2) ); axis image; axis off;
title(sprintf('Denoised Image - Vector IRN. \nSNR: %4.1fdB (noisy SNR: %4.1fdB).\n ', ...
               snr(Ic, IRN_Ic_01L2), snr(Ic, Ic_01L2)));


end % _END_ if( strcmp(example,'l2denoise') || strcmp(example,'all') )

%-----------------------------------------------------------------------------


