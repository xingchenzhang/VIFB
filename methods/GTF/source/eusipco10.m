

% Simulation code for [1], presented in EUSIPCO'10 (http://www.eusipco2010.org/)
%
% [1] Paul Rodriguez "A Non-negative Quadratic Programming approach to Minimize 
%                     the Generalized Vector-Valued Total Variation Functional",  
%                     Proceedings of the European Signal Processing Conference 
%                     (EUSIPCO), (Aalborg, Dinamarca), pp. 314--318, Agosto, 2010 
%  
%  
% NOTE (Jan. 2011) :
%  Results from this simulation are slightly different (computational performance
%  and/or reconstruction quality) from those listed in [1] due to changes 
%  (improvements) in the irntv base code and better selection of parameters
%
% Legal:
%   eusipco.m is based on NUMIPAD (http://numipad.sf.net). NUMIPAD is free 
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

%  example = 'eusipco10';
example = 'l1deconv';
%  example = 'l2deconv';
%  example = 'l1denoise';
%  example = 'l2denoise';
%  example = 'all';



%  Img = 'lena';
%  Img = 'peppers';
Img = 'goldhl.';
%  Img = 'mandrill';
%  Img = 'none';

extreme = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( strcmp(example, 'eusipco10') )

  disp('  Running the "eusipco10" example...');
  disp('  This script reproduces the Exprimental Results reported in');
  disp('  "A Non-negative Quadratic Programming approach to Minimize');
  disp('  the Generalized Vector-Valued Total Variation Functional"');
  disp(' ')

  if( ~strcmp(Img, 'all') ) 

    disp('NOTE:');
    disp('  If you run the "eusipco10" example, you should set Img = "all" ...');
    disp('  setting Img = "all" ');
    Img = 'all';
  end

end



str_all = [];
kernel_sp = [];

if( strcmp(lower(Img), 'all') )

  all = 1;
  if( strcmp(example, 'eusipco10') )
    Img = {'lena', 'goldhl.'};
  else
    Img = {'lena', 'peppers', 'goldhl.', 'mandrill', 'SIPI101', 'SIPI102', 'SIPI05', 'SIPI106', 'SIPI107', ...
           'SIPI202', 'SIPI205', 'SIPI206', 'SIPIhouse', 'SIPI221'};
  end

  NImgs = length(Img);


else

  all = 0;
  NImgs = 1;


  switch lower(Img)

    case{'lena'}
      Ic = double( imread('color_imgs/lena_color_512.png') ) / 255.0;


    case{'goldhl.'}
%        Ic = double( imread('color_imgs/goldhill_color.png') ) / 255.0;
      Ic = double( imread('gray_imgs/goldhill_gray.png') ) / 255.0;

    case{'peppers'}
      Ic = double( imread('color_imgs/peppers_color.png') ) / 255.0;


    case{'p13ang'}
      kernel_sp = fspecial('disk',3.0);
      Ic = double( imread('color_imgs/p13ang_1296x864.png') ) / 255.0;


    case{'fortaleza'}
      kernel_sp = fspecial('disk',3.0);
      Ic = double( imread('color_imgs/saccsayhuman_1296x864.png') ) / 255.0;


    case{'buda'}
      Ic = double( imread('color_imgs/buda.png') ) / 255.0;
      Ig = imadjust(rgb2gray(Ic));

    case{'none'}
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting irnTest code...');
      return;

    otherwise
      error('Not a valid image\n');

  end % _END_ SWITCH(Img)


end % _END_ IF( strcmp(lower(Img), 'all') )


%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));


%%%%%%%%%%%%%%%%%%%%%%%%
%       kernels        %
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
%       kernels        %
%%%%%%%%%%%%%%%%%%%%%%%%

str_all = []; % string to report results


for k = 1:NImgs,


  if( all == 1 )

    kImg = Img{k};

    switch lower(example)

    case{'eusipco10'}

      switch lower(kImg)

        case{'lena'}
          Ic = double( imread('color_imgs/lena_color_512.png') ) / 255.0;


        case{'goldhl.'}
          Ic = double( imread('color_imgs/goldhill_color.png') ) / 255.0;

      end % _END_ switch lower(kImg)

    otherwise
      % FIXME: complete

    end % _END_ switch lower(example)

  else
  
    kImg = Img;

  end % if( all == 1 )


  if( isempty(kernel_sp) ) kernel_len = fspecial('disk',3.2);
  else kernel_len = kernel_sp;
  end

  K = @(x) imfilter(x, kernel_len, 'symmetric','conv');

  KT = @(x) K(x);
  KC = {K, KT};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Blurred & noisy images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  noise_L1deconv = [ 0.1, 0.3, 0.5, 0.7 ];
%    noise_L2deconv = [ 0.05, 0.1 ];
  noise_L2deconv = [ 0.1, 0.2 ];
  noise_L1den = [ 0.1, 0.3, 0.5, 0.7 ];
  noise_L2den = [ 0.05, 0.1 ];


  switch lower(example)

    case{'eusipco10'}

     IcBlur = K(Ic);

      for m=1:length(noise_L1deconv)
        IcBlur_l1{m} = imnoise(IcBlur, 'salt & pepper', noise_L1deconv(m));
      end

      for m=1:length(noise_L2deconv)
        IcBlur_l2{m} = imnoise(IcBlur,'gaussian', 0, (noise_L2deconv(m)*noise_L2deconv(m))*max(IcBlur(:)) );
      end

      for m=1:length(noise_L1den)
        Ic_l1{m} = imnoise(Ic, 'salt & pepper', noise_L1den(m));
      end

      for m=1:length(noise_L2den)
        Ic_l2{m} = imnoise(Ic,'gaussian', 0, (noise_L2den(m)*noise_L2den(m))*max(Ic(:)) );
      end


    case{'l1deconv'}

      IcBlur = K(Ic);

      for m=1:length(noise_L1deconv)
        IcBlur_l1{m} = imnoise(IcBlur, 'salt & pepper', noise_L1deconv(m));
      end


    case{'l2deconv'}

      IcBlur = K(Ic);

      for m=1:length(noise_L2deconv)
        IcBlur_l2{m} = imnoise(IcBlur,'gaussian', 0, (noise_L2deconv(m)*noise_L2deconv(m))*max(IcBlur(:)) );
      end


    case{'l1denoise'}

      for m=1:length(noise_L1den)
        Ic_l1{m} = imnoise(Ic, 'salt & pepper', noise_L1den(m));
      end


    case{'l2denoise'}

      for m=1:length(noise_L2den)
        Ic_l2{m} = imnoise(Ic,'gaussian', 0, (noise_L2den(m)*noise_L2den(m))*max(Ic(:)) );
      end

  end % _END_ switch lower(example)
  




if( strcmp(example,'l1deconv') || strcmp(example,'eusipco10') )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Deconvolved        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% >> Deconvolution via IRN NQP (non-negative quadractic programming) <<

%                  0.1   0.3   0.5   0.7
lambda_irnnqp = [  0.1   0.2  0.3   0.5 ];
lambda_irn    = [  0.1   0.2  0.3   0.5 ];

% ---------- Setup for IRN-NQP  ------------------
pars_nqp = irntvInputPars('l1tv_nqp');

pars_nqp.adapt_epsR   = 1;
pars_nqp.epsR_cutoff  = 0.01;
pars_nqp.adapt_epsF   = 1;
pars_nqp.epsF_cutoff  = 0.1; 

%  pars_nqp.loops_NQP    = 15;
%  pars_nqp.loops        = 15;
inner_loops           = [ 25, 15, 15, 15];
outer_loops           = [ 20, 15, 15, 20];

pars_nqp.gamma_NQP    = 0.5e-3;
pars_nqp.alpha_NQP    = 5e-1;
pars_nqp.vmax_NQP     = 1+1;            % maximum value of output image


% ---------- Setup for IRN (original)  ------------------

pars_irn = irntvInputPars('l1tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 10;


for m=1:length(lambda_irnnqp)


  % --- IRN-NQP ---
  pars_nqp.loops_NQP = inner_loops(m);
  pars_nqp.loops = outer_loops(m);

  
  tic;
  NQP_Ic = irntv(IcBlur_l1{m}+1, KC, lambda_irnnqp(m), pars_nqp)-1;
  tnqp_deconv(m,k) = toc;

  nqp_snr_deconv(m,k) = snr(Ic, NQP_Ic);

  % --- IRN ---

  tic;
  IRN_Ic = irntv(IcBlur_l1{m}, KC, lambda_irn(m), pars_irn);
  tirn_deconv(m,k) = toc;

  irn_snr_deconv(m,k) = snr(Ic, IRN_Ic);

  % ---    ---

  if(SHOW_IMGS)
    figure; imagesc( Normalize(NQP_Ic) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Deconvolve Image - IRN NQP. SNR: %4.1fdB.\n ', ...
                  nqp_snr_deconv(m,k)));

    figure; imagesc( Normalize(IRN_Ic) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Deconvolve Image - IRN. SNR: %4.1fdB.\n ', ...
                  irn_snr_deconv(m,k)));

  end

end


    tnqp_den_mean  = mean(tnqp_deconv,2);
    tirn_den_mean  = mean(tirn_deconv,2);  % FIXME

    nqp_snr_mean = mean(nqp_snr_deconv, 2);
    irn_snr_mean = mean(irn_snr_deconv, 2);  % FIXME


    if( ((k==1) && (all==0)) || (all==1) )
        str = sprintf('\n\nl1 Deconvolve \n\n'); str_all = [ str_all str ];
        str = sprintf('                    SNR (db)          Time (s)     \n');
        str_all = [ str_all str ];
        str = sprintf(' Img    Noise     IRN   IRN-NQP      IRN    IRN-NQP     \n');
        str_all = [ str_all str ];
    end

    for m=1:length(lambda_irnnqp)
      if(m==1) str = sprintf('\n'); str_all = [ str_all str ]; imgname = kImg; else imgname = '     '; end
      str = sprintf('%s\t %1.2f    %5.1f  %5.1f       %5.1f   %5.1f \n', ...
                    imgname, noise_L1deconv(m), ...
                    irn_snr_mean(m), nqp_snr_mean(m), tirn_den_mean(m), tnqp_den_mean(m));
      str_all = [ str_all str ];
    end

    if(all == 0)
      disp(str_all);
    end



end % _END_ if( strcmp(example,'l1deconv') || strcmp(example,'eusipco10') )




if( strcmp(example,'l2deconv') || strcmp(example,'eusipco10') )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L2 Deconvolved        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% >> Deconvolution via IRN NQP (non-negative quadractic programming) <<

%noise level:      0.05   0.1

%  lambda_irnnqp = [  0.01   0.04  ];
%  lambda_irn    = [  0.01   0.04  ];

lambda_irnnqp = [  0.04   0.06  ];
lambda_irn    = [  0.04   0.06  ];


% ---------- Setup for IRN-NQP  ------------------
pars_nqp = irntvInputPars('l2tv_nqp');

pars_nqp.adapt_epsR   = 1;
pars_nqp.epsR_cutoff  = 0.01;
pars_nqp.adapt_epsF   = 1;
pars_nqp.epsF_cutoff  = 0.1; 

inner_loops           = [ 20, 20 ];
outer_loops           = [ 6, 6 ];

pars_nqp.gamma_NQP    = 0.5e-3;
pars_nqp.alpha_NQP    = 5e-1;
pars_nqp.vmax_NQP     = 1+1;            % maximum value of output image


% ---------- Setup for IRN (original)  ------------------

pars_irn = irntvInputPars('l2tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 6;


for m=1:length(lambda_irnnqp)

%  if(m < 2) continue; end

  % --- IRN-NQP ---
  pars_nqp.loops_NQP = inner_loops(m);
  pars_nqp.loops = outer_loops(m);

  
  tic;
  NQP_Ic = irntv(IcBlur_l2{m}+1, KC, lambda_irnnqp(m), pars_nqp)-1;
  tnqp_deconv(m,k) = toc;

  nqp_snr_deconv(m,k) = snr(Ic, NQP_Ic);

  % --- IRN ---

  tic;
  IRN_Ic = irntv(IcBlur_l2{m}, KC, lambda_irn(m), pars_irn);
  tirn_deconv(m,k) = toc;

  sum(IRN_Ic(:) < 0)

  irn_snr_deconv(m,k) = snr(Ic, IRN_Ic);

  % ---    ---

  if(SHOW_IMGS)
    figure; imagesc( Normalize(NQP_Ic) ); 
    axis image; axis off; colormap gray;
%      title(sprintf('Deconvolve Image - IRN NQP. SNR: %4.1fdB.\n ', ...
%                    nqp_snr_deconv(m,k)));

    figure; imagesc( Normalize(IRN_Ic) ); 
    axis image; axis off; colormap gray;
%      title(sprintf('Deconvolve Image - IRN. SNR: %4.1fdB.\n ', ...
%                    irn_snr_deconv(m,k)));

  end

end


    tnqp_den_mean  = mean(tnqp_deconv,2);
    tirn_den_mean  = mean(tirn_deconv,2);  % FIXME

    nqp_snr_mean = mean(nqp_snr_deconv, 2);
    irn_snr_mean = mean(irn_snr_deconv, 2);  % FIXME


    if( ((k==1) && (all==0)) || (all==1) )
        str = sprintf('\n\nl2 Deconvolve \n\n'); str_all = [ str_all str ];
        str = sprintf('                    SNR (db)          Time (s)     \n');
        str_all = [ str_all str ];
        str = sprintf(' Img    Noise     IRN   IRN-NQP      IRN    IRN-NQP     \n');
        str_all = [ str_all str ];
    end

    for m=1:length(lambda_irnnqp)
      if(m==1) str = sprintf('\n'); str_all = [ str_all str ]; imgname = kImg; else imgname = '     '; end
      str = sprintf('%s\t %1.2f    %5.1f  %5.1f       %5.1f   %5.1f \n', ...
                    imgname, noise_L2deconv(m), ...
                    irn_snr_mean(m), nqp_snr_mean(m), tirn_den_mean(m), tnqp_den_mean(m));
      str_all = [ str_all str ];
    end

    if(all == 0)
      disp(str_all);
    end

end % _END_ if( strcmp(example,'l2deconv') || strcmp(example,'eusipco10') )



if( strcmp(example,'l1denoise') || strcmp(example,'eusipco10') )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         L1 Denoised         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% >> Denoising via IRN NQP (non-negative quadractic programming) <<

%                  0.1   0.3   0.5   0.7
lambda_irnnqp = [  0.85   1.1  1.2   1.4 ];
lambda_irn    = [  0.85   1.1  1.2   1.4 ];

% ---------- Setup for IRN-NQP  ------------------
pars_nqp = irntvInputPars('l1tv_nqp');

pars_nqp.adapt_epsR   = 1;
pars_nqp.epsR_cutoff  = 0.01;
pars_nqp.adapt_epsF   = 1;
pars_nqp.epsF_cutoff  = 0.1; 

%  pars_nqp.loops_NQP    = 15;
%  pars_nqp.loops        = 15;
inner_loops           = [  5,  5, 12, 12];
outer_loops           = [ 15, 15, 15, 20];

pars_nqp.gamma_NQP    = 0.5e-3;
pars_nqp.alpha_NQP    = 5e-1;
pars_nqp.vmax_NQP     = 1+1;            % maximum value of output image


% ---------- Setup for IRN (original)  ------------------

pars_irn = irntvInputPars('l1tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.1;

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 8;


for m=1:length(lambda_irnnqp)


  % --- IRN-NQP ---
  pars_nqp.loops_NQP = inner_loops(m);
  pars_nqp.loops = outer_loops(m);

  tic;
  NQP_Ic = irntv(Ic_l1{m}+1, {}, lambda_irnnqp(m), pars_nqp)-1;
  tnqp_denoise(m,k) = toc;

  nqp_snr_denoise(m,k) = snr(Ic, NQP_Ic);

  % --- IRN ---

  tic;
  IRN_Ic = irntv(Ic_l1{m}, {}, lambda_irn(m), pars_irn);
  tirn_denoise(m,k) = toc;

  irn_snr_denoise(m,k) = snr(Ic, IRN_Ic);

  % ---    ---

  if(SHOW_IMGS)
    figure; imagesc( Normalize(NQP_Ic) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN NQP. SNR: %4.1fdB.\n ', ...
                  nqp_snr_denoise(m,k)));

    figure; imagesc( Normalize(IRN_Ic) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN. SNR: %4.1fdB.\n ', ...
                  irn_snr_denoise(m,k)));

  end

end


    tnqp_den_mean  = mean(tnqp_denoise,2);
    tirn_den_mean  = mean(tirn_denoise,2);  % FIXME

    nqp_snr_mean = mean(nqp_snr_denoise, 2);
    irn_snr_mean = mean(irn_snr_denoise, 2);  % FIXME


    if( ((k==1) && (all==0)) || (all==1) )
        str = sprintf('\n\nl1 Denoise \n\n'); str_all = [ str_all str ];
        str = sprintf('                    SNR (db)          Time (s)     \n');
        str_all = [ str_all str ];
        str = sprintf(' Img    Noise     IRN   IRN-NQP      IRN    IRN-NQP     \n');
        str_all = [ str_all str ];
    end

    for m=1:length(lambda_irnnqp)
      if(m==1) str = sprintf('\n'); str_all = [ str_all str ]; imgname = kImg; else imgname = '     '; end
      str = sprintf('%s\t %1.2f    %5.1f  %5.1f       %5.1f   %5.1f \n', ...
                    imgname, noise_L1den(m), ...
                    irn_snr_mean(m), nqp_snr_mean(m), tirn_den_mean(m), tnqp_den_mean(m));
      str_all = [ str_all str ];
    end

    if(all == 0)
      disp(str_all);
    end


end % _END_ if( strcmp(example,'l1deconv') || strcmp(example,'eusipco10') )



if( strcmp(example,'l2denoise') || strcmp(example,'eusipco10') )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         L2 Denoised         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% >> Denoising via IRN NQP (non-negative quadractic programming) <<

%noise level:      0.05   0.1

lambda_irnnqp = [  0.035   0.09  ];
lambda_irn    = [  0.035   0.09  ];



% ---------- Setup for IRN-NQP  ------------------
pars_nqp = irntvInputPars('l2tv_nqp');

pars_nqp.adapt_epsR   = 1;
pars_nqp.epsR_cutoff  = 0.01;
pars_nqp.adapt_epsF   = 1;
pars_nqp.epsF_cutoff  = 0.1; 

%  pars_nqp.loops_NQP    = 15;
%  pars_nqp.loops        = 15;
inner_loops           = [ 5, 5 ];
outer_loops           = [ 5, 5 ];

pars_nqp.gamma_NQP    = 5e-3;
pars_nqp.alpha_NQP    = 5e-1;
pars_nqp.vmax_NQP     = 1+1;            % maximum value of output image


% ---------- Setup for IRN (original)  ------------------

pars_irn = irntvInputPars('l2tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;

pars_irn.pcgtol_ini = 1e-4;

irn_loops           = [3, 3];
%  pars_irn.loops      = 1;


for m=1:length(lambda_irnnqp)


  % --- IRN-NQP ---
  pars_nqp.loops_NQP = inner_loops(m);
  pars_nqp.loops = outer_loops(m);

  tic;
  NQP_Ic = irntv(Ic_l2{m}+1, {}, lambda_irnnqp(m), pars_nqp)-1;
  tnqp_denoise(m,k) = toc;

  nqp_snr_denoise(m,k) = snr(Ic, NQP_Ic);

  % --- IRN ---
  pars_irn.loops = irn_loops(m);
  tic;
  IRN_Ic = irntv(Ic_l2{m}, {}, lambda_irn(m), pars_irn);
  tirn_denoise(m,k) = toc;

  irn_snr_denoise(m,k) = snr(Ic, IRN_Ic);

  % ---    ---

  if(SHOW_IMGS)
    figure; imagesc( Normalize(NQP_Ic) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN NQP. SNR: %4.1fdB.\n ', ...
                  nqp_snr_denoise(m,k)));

    figure; imagesc( Normalize(IRN_Ic) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN. SNR: %4.1fdB.\n ', ...
                  irn_snr_denoise(m,k)));

  end

end


    tnqp_den_mean  = mean(tnqp_denoise,2);
    tirn_den_mean  = mean(tirn_denoise,2);  % FIXME

    nqp_snr_mean = mean(nqp_snr_denoise, 2);
    irn_snr_mean = mean(irn_snr_denoise, 2);  % FIXME


    if( ((k==1) && (all==0)) || (all==1) )
        str = sprintf('\n\nl2 Denoise \n\n'); str_all = [ str_all str ];
        str = sprintf('                    SNR (db)          Time (s)     \n');
        str_all = [ str_all str ];
        str = sprintf(' Img    Noise     IRN   IRN-NQP      IRN    IRN-NQP     \n');
        str_all = [ str_all str ];
    end

    for m=1:length(lambda_irnnqp)
      if(m==1) str = sprintf('\n'); str_all = [ str_all str ]; imgname = kImg; else imgname = '     '; end
      str = sprintf('%s\t %1.2f    %5.1f  %5.1f       %5.1f   %5.1f \n', ...
                    imgname, noise_L2den(m), ...
                    irn_snr_mean(m), nqp_snr_mean(m), tirn_den_mean(m), tnqp_den_mean(m));
      str_all = [ str_all str ];
    end

    if(all == 0)
      disp(str_all);
    elseif(k==NImgs)
      disp(str_all)
    end



end % _END_ if( strcmp(example,'l2deconv') || strcmp(example,'eusipco10') )




end % _END_ FOR(k=1:NImgs)

