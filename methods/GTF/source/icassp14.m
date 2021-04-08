function[] = icassp14()



SHOW_IMGS = false;
%  SHOW_IMGS = true;


%  example = 'icassp14'
%  example = 'l2TV';
example = 'l1TV';
%  example = 'l1l2TV';


%  Img = 'Clena';
%  Img = 'peppers';
Img = 'lena';
%  Img = 'barb';
%  Img = 'cman';
%  Img = 'Cboats';

%  Img = 'all';


Nloops = 1; % This value will be change to 10 if example == icassp12



%%%%%%%%%%%%%%%%%%%%%%%%
%   Initial  setup     %
%%%%%%%%%%%%%%%%%%%%%%%%

[SSIM_CODE, Nloops, Img, NImgs all] = initial_setup(example, Nloops, Img);

str_all = [];







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          loop over Img
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NImgs

for t = 1:NImgs,

  t

  if( all == 1 )
    tImg = Img{t}
  else
    tImg = lower(Img)
  end

  [Ig, Color, SSIM_CODE] = readImage( tImg, SSIM_CODE );

  disp(' '); disp(tImg);
  Color 
  disp(' ');


  [Nrows Ncols Ndims] = size(Ig);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   Set noise level and Reg. Parameters     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch lower(example)

      % ----------
      case{'l2tv'}

        sigmaGauss = [0.05, 0.1, 0.2 ];

        if(Ndims == 1)
          lambda_l2 = [ 0.035, 0.05, 0.15];
        else
          lambda_l2 = 1.5*[ 0.075, 0.1, 0.25];
        end


      % ----------
      case{'l1tv'}

        sigmaSP = [0.1, 0.3, 0.5 ];

        if(Ndims == 1)
          lambda_l1 = [ 0.85, 1.2, 1.4];
        else
          lambda_l1 = 1.5*[ 0.85, 1.2, 1.4];
        end

      % ------------
      case{'l1l2tv'}


    end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  loops for average results    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for loops=1:Nloops


    %%%%%%%%%%%%%%%%%%%%%%%%
    %  Generate noisy images  
    %%%%%%%%%%%%%%%%%%%%%%%%
        
    switch lower(example)

      % ----------
      case{'l2tv'}

        clear InCell;
        NGauss = length(sigmaGauss);
    
        for k=1:NGauss,
    
          InCell{k} = imnoise(Ig,'gaussian', 0, (sigmaGauss(k)^2) );          

        end

      % ----------
      case{'l1tv'}

        clear InCell;
        Nsp = length(sigmaSP);
    
        for k=1:Nsp,
    
          InCell{k} = imnoise(Ig,'salt & pepper', sigmaSP(k) );          

        end

      % ----------
      case{'l1l2tv'}

    end % _END_ SWITCH

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    icassp14_defs

    % loop over noisy images
    for k=1:length(InCell)


      switch lower(example)

        % ----------
        case{'l2tv'}
%            [tv_time(:,loops), tv_snr(:,loops), tv_psnr(:,loops), tv_ssim(:,loops)] = ...
%                               denoise_l2tv(InCell{k}, k, lambda_l2(k), sigmaGauss(k), Ig, SSIM_CODE);
%  
%            str_all = report(str_all, sigmaGauss, tv_snr, tv_ssim, tv_psnr, tv_time, tImg, k, loops, Nloops, Color, 0);

          [tv_time(:,:,loops), tv_snr(:,:,loops), tv_psnr(:,:,loops), tv_ssim(:,:,loops)] = ...
                             denoise_l2tv_loop(InCell{k}, k, lambda_l2(k), sigmaGauss(k), Ig, SSIM_CODE);
	  figure(k);
	  plot( tv_time(TVL2_WTHRESH_1,:,loops), tv_snr(TVL2_WTHRESH_1,:,loops), 'r-*' ); hold on;
	  plot( tv_time(TVL2_WTHRESH_2,:,loops), tv_snr(TVL2_WTHRESH_2,:,loops), 'b-*' ); 
	  plot( tv_time(TVL2_WMIL_1,:,loops), tv_snr(TVL2_WMIL_1,:,loops), 'g-*' ); 
	  plot( tv_time(TVL2_WMIL_2,:,loops), tv_snr(TVL2_WMIL_2,:,loops), 'c-*' ); 
	  plot( tv_time(TVL2_SB_1,:,loops), tv_snr(TVL2_SB_1,:,loops), 'k-*' ); 

	  title(sprintf('L2-TV, sigma = %1.2f', sigmaGauss(k) ));
	  legend(TVL2_CELLNAMES{1}, TVL2_CELLNAMES{2}, TVL2_CELLNAMES{3}, TVL2_CELLNAMES{4}, TVL2_CELLNAMES{5}, ...
		 'Location', 'SouthEast');

        % ----------
        case{'l1tv'}

%            [tv_time(:,loops), tv_snr(:,loops), tv_psnr(:,loops), tv_ssim(:,loops)] = ...
%                               denoise_l1tv(InCell{k}, k, lambda_l1(k), sigmaSP(k), Ig, SSIM_CODE);
%  
%            str_all = report(str_all, sigmaSP, tv_snr, tv_ssim, tv_psnr, tv_time, tImg, k, loops, Nloops, Color, 1);


          [tv_time(:,:,loops), tv_snr(:,:,loops), tv_psnr(:,:,loops), tv_ssim(:,:,loops)] = ...
                             denoise_l1tv_loop(InCell{k}, k, lambda_l1(k), sigmaSP(k), Ig, SSIM_CODE);
	  figure(k);
	  plot( tv_time(TVL1_WTHRESH_SUB_1,:,loops), tv_snr(TVL1_WTHRESH_SUB_1,:,loops), 'm-*' ); hold on;
	  plot( tv_time(TVL1_WTHRESH_SUB_2,:,loops), tv_snr(TVL1_WTHRESH_SUB_2,:,loops), 'y-*' ); 
	  plot( tv_time(TVL1_WTHRESH_1,:,loops), tv_snr(TVL1_WTHRESH_1,:,loops), 'r-*' ); 
	  plot( tv_time(TVL1_WTHRESH_2,:,loops), tv_snr(TVL1_WTHRESH_2,:,loops), 'b-*' ); 
	  plot( tv_time(TVL1_WMIL_1,:,loops), tv_snr(TVL1_WMIL_1,:,loops), 'g-*' ); 
	  plot( tv_time(TVL1_WMIL_2,:,loops), tv_snr(TVL1_WMIL_2,:,loops), 'c-*' ); 
	  plot( tv_time(TVL1_SB_1,:,loops), tv_snr(TVL1_SB_1,:,loops), 'k-*' ); 

	  title(sprintf('L1-TV, sigma = %1.2f', sigmaSP(k) ));
	  legend(TVL1_CELLNAMES{1}, TVL1_CELLNAMES{2}, TVL1_CELLNAMES{3}, TVL1_CELLNAMES{4}, TVL1_CELLNAMES{5}, ...
		 TVL1_CELLNAMES{6}, TVL1_CELLNAMES{7}, 'Location', 'SouthEast');



        % ----------
        case{'l1l2tv'}


      end % _END_ SWITCH(example)




    end % _END_ FOR(length(InCell))


  end % FOR(loops)

  disp(str_all)

end % FOR(t)  -- test Image

  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[str_all] = report(str_all, sigmaGauss, tv_snr, tv_ssim, tv_psnr, tv_time, tImg, k, loops, Nloops, Color, tvcase)

icassp14_defs

switch(tvcase)
  case{0}
    TV_CELLNAMES = TVL2_CELLNAMES;
%  
  case{1}
    TV_CELLNAMES = TVL1_CELLNAMES;
%  
end

  if(loops == Nloops)

      if(k==1)

        if(Color == 0) str = sprintf('Img: grayscale %s (denoising)\n', lower(tImg));
        else str = sprintf('Img: color %s (denoising)\n', lower(tImg));
        end

        str_all = [str_all str];

        str = sprintf(' sigma \t variant \t\t\t SNR \t SSIM \t PSNR    Time\n');
        str_all = [str_all str];

      end % _END_ IF(k==1)

      str = sprintf(' \n');
      str_all = [str_all str];

      str = sprintf(' %0.2f     -- %s \t\t %1.2f \t %1.2f \t %1.2f    %1.2f\n', ...
		    sigmaGauss(k), TV_CELLNAMES{1}, ...
                    mean( tv_snr(1,:) ), mean( tv_ssim(1,:) ), ...
                    mean( tv_psnr(1,:) ), mean( tv_time(1,:) ) );
      str_all = [ str_all str ];

      for m=2:length( tv_snr(:,Nloops) )

        str = sprintf('          -- %s \t\t %1.2f \t %1.2f \t %1.2f    %1.2f\n', ...
		      TV_CELLNAMES{m}, ...
                      mean( tv_snr(m,:) ), mean( tv_ssim(m,:) ), ...
                    mean( tv_psnr(m,:) ), mean( tv_time(m,:) ) );
        str_all = [ str_all str ];
      end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[SSIM_CODE, Nloops, Img, NImgs, all] = initial_setup(example, Loops, myImg);


% Check for SSIM
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

% Setup simulation
if( strcmp(example, 'icassp14') )

  disp('  Running the "icassp14" example...');
  disp('  This script reproduces the Exprimental Results of ');
  disp('  "paper Title"');
  disp(' ')
  disp('  This script takes sometime (ten-trial average)');
  disp(' ')
  Nloops = 10;

  if( ~strcmp(myImg, 'all') ) 

    disp('NOTE:');
    disp('  If you run the "icassp14" example, you should set Img = "all" ...');
    disp('  setting Img = "all" ');
    Img = 'all';
  end

else
  Nloops = Loops;
  Img    = myImg;
end


if( strcmp(lower(Img), 'all') )

  all = 1;
  if( strcmp(example, 'icassp14') )
    Img = {'lena', 'goldhl', 'cman512', 'barb', 'Clena', 'Cbarb', 'Cboats', 'Cgoldhl' };
  else
    Img = {'lena', 'peppers', 'goldhl.', 'cman512', 'Clena', 'Cpeppers', 'Cmandrill'};
  end

  NImgs = length(Img);

else

  all = 0;
  NImgs = 1;

end  % _END_ IF(strcmp(lower(Img))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Ig, Color, SSIM_CODE] = readImage( tImg, mySSIM_CODE )

  switch lower(tImg)

    case{'lena'}
      Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;
      Color = 0;
      SSIM_CODE = mySSIM_CODE;

    case{'peppers'}
      Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;
      Color = 0;
      SSIM_CODE = mySSIM_CODE;

    case{'cman'}
      Ig = double( imread('gray_imgs/cameraman_256x256.tiff') ) / 255;
      Color = 0;
      SSIM_CODE = mySSIM_CODE;

    case{'cman512'}
      Ig = double( imread('gray_imgs/cameraman_512x512.tiff') ) / 255;
      Color = 0;
      SSIM_CODE = mySSIM_CODE;

    case{'goldhl'}
      Ig = double( imread('gray_imgs/goldhill_gray.png') ) / 255;
      Color = 0;
      SSIM_CODE = mySSIM_CODE;

    case{'barb'}
      Ig = double( imread('gray_imgs/barbara_gray512.png') ) / 255;
      Color = 0;
      SSIM_CODE = mySSIM_CODE;

    case{'clena'}
      Ig = double( imread('color_imgs/lena_color_512.png') ) / 255;
      SSIM_CODE = 0;
      Color = 1;

    case{'cbarb'}
      Ig = double( imread('color_imgs/barbara_color.png') ) / 255;
      SSIM_CODE = 0;
      Color = 1;

    case{'cboats'}
      Ig = double( imread('color_imgs/boats_color.png') ) / 255;
      SSIM_CODE = 0;
      Color = 1;

    case{'cgoldhl'}
      Ig = double( imread('color_imgs/goldhill_color.png') ) / 255;
      SSIM_CODE = 0;
      Color = 1;

    case{'none'}
      disp(' ');
      disp('Select a test image (see code for an example)...');
      disp(' ');
      disp('Exiting ICASSP14 test code...');
      return;

    otherwise
      error('Not a valid image\n');

  end  % _END_ SWITCH(Img)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[time_l2, snr_l2, psnr_l2, ssim_l2] = denoise_l2tv(InCell, k, lambda, sigma, Inoiseless, SSIM_CODE)

% General setup for l2-TV

nmpdef
icassp14_defs


pars_irn = irntvInputPars('l2tv');

pars_irn.loops  = 5;
pars_irn.U0     = InCell;



time_l2 = zeros(4,1);
snr_l2  = zeros(4,1);
psnr_l2 = zeros(4,1);
ssim_l2 = zeros(4,1);


% =========================
% >>> WEIGHTS_THRESHOLD <<<
% =========================


% -----------------------------------------------
% --- Adapt cutoff value, adapt PCG tolerance ---

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.adaptPCGtol  = 1;

pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;



t = tic;
  I_L2 = irntv_nv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WTHRESH_1, 1) = toc(t);


  snr_l2(TVL2_WTHRESH_1, 1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WTHRESH_1, 1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WTHRESH_1, 1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WTHRESH_1, 1) = NaN;
  end




% -----------------------------------------------
% --- Adapt cutoff value, fixed PCG tolerance ---

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;



t = tic;
  I_L2 = irntv_nv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WTHRESH_2, 1) = toc(t);

  snr_l2(TVL2_WTHRESH_2, 1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WTHRESH_2, 1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WTHRESH_2, 1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WTHRESH_2, 1) = NaN;
  end




% ===================
% >>> WEIGHTS_MIL <<<   -- Matrix Inversion Lemma
% ===================

% ---------------------------
% --- adapt PCG tolerance ---

pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;
pars_irn.adaptPCGtol   = 1;



t = tic;
I_L2 = irntv_nv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WMIL_1,1) = toc(t);

  snr_l2(TVL2_WMIL_1,1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WMIL_1,1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WMIL_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WMIL_1,1) = NaN;
  end




% ---------------------------
% --- fixed PCG tolerance ---

pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;

pars_irn.adaptPCGtol   = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance


t = tic;
  I_L2 = irntv_nv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WMIL_2,1) = toc(t);

  snr_l2(TVL2_WMIL_2,1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WMIL_2,1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WMIL_2,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WMIL_2,1) = NaN;
  end




% ===================
% >>>    SB-TV    <<<   -- Split-Breagman / ADMM
% ===================

% ----------------------------------
% --- fixed gamma (twice lambda) ---


t = tic;
  I_L2 = sbtv2_BW(InCell, {}, lambda, 2*lambda, pars_irn.loops, pars_irn.U0);
time_l2(TVL2_SB_1,1) = toc(t);

  snr_l2(TVL2_SB_1,1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_SB_1,1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_SB_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_SB_1,1) = NaN;
  end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[time_l2, snr_l2, psnr_l2, ssim_l2] = denoise_l2tv_loop(InCell, k, lambda, sigma, Inoiseless, SSIM_CODE)

% General setup for l2-TV

nmpdef
icassp14_defs


pars_irn = irntvInputPars('l2tv');

% NOTE: this should be a parameter
Loops               = 6;
%  pars_irn.loops      = 5;

pars_irn.lambda_ini = 2*lambda; 	% force IRN and SB to have the same initial sol.


time_l2 = zeros(4,Loops);
snr_l2  = zeros(4,Loops);
psnr_l2 = zeros(4,Loops);
ssim_l2 = zeros(4,Loops);


% =========================
% >>> WEIGHTS_THRESHOLD <<<
% =========================


% -----------------------------------------------
% --- Adapt cutoff value, adapt PCG tolerance ---

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.adaptPCGtol  = 1;

pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;


% 1st iter (one outer loop)
pars_irn.loops  = 0;
pars_irn.U0     = InCell;

t = tic;
  I_L2 = irntv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WTHRESH_1, 1) = toc(t);


  snr_l2(TVL2_WTHRESH_1, 1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WTHRESH_1, 1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WTHRESH_1, 1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WTHRESH_1, 1) = NaN;
  end


for m=2:Loops,

  pars_irn.loops  = m-1;		


  t = tic;
    I_L2 = irntv(InCell, {}, lambda, pars_irn);
  time_l2(TVL2_WTHRESH_1, m) = toc(t);


  snr_l2(TVL2_WTHRESH_1, m)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WTHRESH_1, m) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WTHRESH_1, m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WTHRESH_1, m) = NaN;
  end


end


% -----------------------------------------------
% --- Adapt cutoff value, fixed PCG tolerance ---

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;


% 1st iter (one outer loop)
pars_irn.loops = 0;
pars_irn.U0     = InCell;

t = tic;
  I_L2 = irntv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WTHRESH_2, 1) = toc(t);

  snr_l2(TVL2_WTHRESH_2, 1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WTHRESH_2, 1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WTHRESH_2, 1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WTHRESH_2, 1) = NaN;
  end


for m = 2: Loops,

  pars_irn.loops = m-1;

  t = tic;
    I_L2 = irntv(InCell, {}, lambda, pars_irn);
  time_l2(TVL2_WTHRESH_2, m) = toc(t);

  snr_l2(TVL2_WTHRESH_2, m)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WTHRESH_2, m) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WTHRESH_2, m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WTHRESH_2, m) = NaN;
  end


end


% ===================
% >>> WEIGHTS_MIL <<<   -- Matrix Inversion Lemma
% ===================

% ---------------------------
% --- adapt PCG tolerance ---

pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;
pars_irn.adaptPCGtol   = 1;


% 1st iter (one outer loop)
pars_irn.loops  = 0;
pars_irn.U0     = InCell;

t = tic;
I_L2 = irntv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WMIL_1,1) = toc(t);

  snr_l2(TVL2_WMIL_1,1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WMIL_1,1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WMIL_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WMIL_1,1) = NaN;
  end


for m = 2:Loops,

  pars_irn.loops  = m-1;


  t = tic;
  I_L2 = irntv(InCell, {}, lambda, pars_irn);
  time_l2(TVL2_WMIL_1,m) = toc(t);

  snr_l2(TVL2_WMIL_1,m)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WMIL_1,m) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WMIL_1,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WMIL_1,m) = NaN;
  end

end


% ---------------------------
% --- fixed PCG tolerance ---

pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;

pars_irn.adaptPCGtol   = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

% 1st iter (one outer loop)
pars_irn.loops  = 0;
pars_irn.U0     = InCell;

t = tic;
  I_L2 = irntv(InCell, {}, lambda, pars_irn);
time_l2(TVL2_WMIL_2,1) = toc(t);

  snr_l2(TVL2_WMIL_2,1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WMIL_2,1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WMIL_2,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WMIL_2,1) = NaN;
  end


for m=2:Loops,

  pars_irn.loops  = m-1;		


  t = tic;
    I_L2 = irntv(InCell, {}, lambda, pars_irn);
  time_l2(TVL2_WMIL_2,m) = toc(t);

  snr_l2(TVL2_WMIL_2,m)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_WMIL_2,m) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_WMIL_2,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_WMIL_2,m) = NaN;
  end


end


% ===================
% >>>    SB-TV    <<<   -- Split-Breagman / ADMM
% ===================

% ----------------------------------
% --- fixed gamma (twice lambda) ---

% 1st iter (one outer loop)
pars_irn.loops = 0;
pars_irn.U0    = InCell;

t = tic;
  I_L2 = sbtv2_BW(InCell, {}, lambda, 2*lambda, pars_irn.loops, pars_irn.U0);
time_l2(TVL2_SB_1,1) = toc(t);

  snr_l2(TVL2_SB_1,1)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_SB_1,1) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_SB_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_SB_1,1) = NaN;
  end


for m=2:Loops

  pars_irn.loops = m-1;


  t = tic;
    I_L2 = sbtv2_BW(InCell, {}, lambda, 2*lambda, pars_irn.loops, pars_irn.U0);
  time_l2(TVL2_SB_1,m) = toc(t);

  snr_l2(TVL2_SB_1,m)  = snr(Inoiseless, I_L2);
  psnr_l2(TVL2_SB_1,m) = psnr(Inoiseless, I_L2);

  if(SSIM_CODE)
    ssim_l2(TVL2_SB_1,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L2));
  else
    ssim_l2(TVL2_SB_1,m) = NaN;
  end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[time_l1, snr_l1, psnr_l1, ssim_l1] = denoise_l1tv(InCell, k, lambda, sigma, Inoiseless, SSIM_CODE)

% General setup for l1-TV

nmpdef
icassp14_defs


pars_irn = irntvInputPars('l1tv');

pars_irn.loops      = 5;
pars_irn.U0         = InCell;

time_l2 = zeros(7,1);
snr_l2  = zeros(7,1);
psnr_l2 = zeros(7,1);
ssim_l2 = zeros(7,1);



% =============================
% >>> FIDELITY_SUBSTITUTION <<<
% =============================


% -----------------------------------------------
% --- Adapt cutoff value, fixed PCG tolerance ---

pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance


t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_SUB_1) = toc(t);


  snr_l1(TVL1_WTHRESH_SUB_1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_SUB_1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_SUB_1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_SUB_1) = NaN;
  end

% -----------------------------------------------
% --- Adapt cutoff value, adapt PCG tolerance ---

pars_irn = irntvInputPars('l1tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 5;
pars_irn.U0         = InCell;

pars_irn.variant       = NMP_TV_SUBSTITUTION;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_SUB_2) = toc(t);


  snr_l1(TVL1_WTHRESH_SUB_2)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_SUB_2) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_SUB_2) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_SUB_2) = NaN;
  end



% =========================
% >>> WEIGHTS_THRESHOLD <<<
% =========================

% -----------------------------------------------
% --- Adapt cutoff value, fixed PCG tolerance ---

pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_1) = toc(t);


  snr_l1(TVL1_WTHRESH_1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_1) = NaN;
  end

% -----------------------------------------------
% --- Adapt cutoff value, adapt PCG tolerance ---

pars_irn.variant       = NMP_TV_STANDARD;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;
pars_irn.U0            = [];    % Initial sol. is empty: In this case, this give better results

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_2) = toc(t);


  snr_l1(TVL1_WTHRESH_2)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_2) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_2) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_2) = NaN;
  end


% ===================
% >>> WEIGHTS_MIL <<<   -- Matrix Inversion Lemma
% ===================

% ---------------------------
% --- adapt PCG tolerance ---

pars_irn = irntvInputPars('l1tv');

pars_irn.pcgtol_ini   = 1e-4;
pars_irn.loops        = 5;
pars_irn.U0           = InCell;
pars_irn.sbstflg      = 0;
pars_irn.adapt_epsR   = 0;
pars_irn.adapt_epsF   = 0;


pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WMIL_1) = toc(t);


  snr_l1(TVL1_WMIL_1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WMIL_1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WMIL_1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WMIL_1) = NaN;
  end

% ---------------------------
% --- fixed PCG tolerance ---


pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WMIL_2) = toc(t);


  snr_l1(TVL1_WMIL_2)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WMIL_2) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WMIL_2) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WMIL_2) = NaN;
  end


% ===================
% >>>    SB-TV    <<<   -- Split-Breagman / ADMM
% ===================

% ----------------------------------
% --- fixed gamma (twice lambda) ---

%  sbOpts.MaxMainIter = pars_irn.loops;
%  sbOpts.BoundRange = [-Inf Inf];
%  sbOpts.gamma = 2*lambda;

t = tic;
%    I_L1 = sbl1tv2_BW(InCell, {}, lambda, sbOpts);
  I_L1 = sbl1tv2_BW(InCell, {}, lambda);
time_l1(TVL1_SB_1) = toc(t);

  snr_l1(TVL1_SB_1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_SB_1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_SB_1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_SB_1) = NaN;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[time_l1, snr_l1, psnr_l1, ssim_l1] = denoise_l1tv_loop(InCell, k, lambda, sigma, Inoiseless, SSIM_CODE)

% General setup for l1-TV

nmpdef
icassp14_defs


pars_irn = irntvInputPars('l1tv');

Loops    	    = 6;
%  pars_irn.loops      = 5;

pars_irn.U0         = InCell;
pars_irn.lambda_ini = 2*lambda;


time_l2 = zeros(7,Loops);
snr_l2  = zeros(7,Loops);
psnr_l2 = zeros(7,Loops);
ssim_l2 = zeros(7,Loops);



% =============================
% >>> FIDELITY_SUBSTITUTION <<<
% =============================


% -----------------------------------------------
% --- Adapt cutoff value, fixed PCG tolerance ---

pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

pars_irn.loops = 0;

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_SUB_1,1) = toc(t);


  snr_l1(TVL1_WTHRESH_SUB_1,1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_SUB_1,1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_SUB_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_SUB_1,1) = NaN;
  end


for m=2:Loops

pars_irn.loops = m-1;

  t = tic;
    I_L1 = irntv(InCell, {}, lambda, pars_irn);
  time_l1(TVL1_WTHRESH_SUB_1,m) = toc(t);


  snr_l1(TVL1_WTHRESH_SUB_1,m)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_SUB_1,m) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_SUB_1,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_SUB_1,m) = NaN;
  end

end


% -----------------------------------------------
% --- Adapt cutoff value, adapt PCG tolerance ---

pars_irn = irntvInputPars('l1tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 5;
pars_irn.U0         = InCell;

pars_irn.variant       = NMP_TV_SUBSTITUTION;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;


pars_irn.loops = 0;
t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_SUB_2,1) = toc(t);


  snr_l1(TVL1_WTHRESH_SUB_2,1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_SUB_2,1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_SUB_2,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_SUB_2,1) = NaN;
  end


for m=2:Loops,

  pars_irn.loops = m-1;

  t = tic;
    I_L1 = irntv(InCell, {}, lambda, pars_irn);
  time_l1(TVL1_WTHRESH_SUB_2,m) = toc(t);


  snr_l1(TVL1_WTHRESH_SUB_2,m)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_SUB_2,m) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_SUB_2,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_SUB_2,m) = NaN;
  end

end


% =========================
% >>> WEIGHTS_THRESHOLD <<<
% =========================

% -----------------------------------------------
% --- Adapt cutoff value, fixed PCG tolerance ---

pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

pars_irn.loops = 0;

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_1,1) = toc(t);


  snr_l1(TVL1_WTHRESH_1,1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_1,1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_1,1) = NaN;
  end


for m=2:Loops,

  pars_irn.loops = m-1;

  t = tic;
    I_L1 = irntv(InCell, {}, lambda, pars_irn);
  time_l1(TVL1_WTHRESH_1,m) = toc(t);


  snr_l1(TVL1_WTHRESH_1,m)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_1,m) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_1,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_1,m) = NaN;
  end


end


% -----------------------------------------------
% --- Adapt cutoff value, adapt PCG tolerance ---

pars_irn.variant       = NMP_TV_STANDARD;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;
%  pars_irn.U0            = [];    % Initial sol. is empty: In this case, this give better results


pars_irn.loops = 0;

t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WTHRESH_2,1) = toc(t);


  snr_l1(TVL1_WTHRESH_2,1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_2,1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_2,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_2,1) = NaN;
  end


for m = 2:Loops

  pars_irn.loops = m-1;

  t = tic;
    I_L1 = irntv(InCell, {}, lambda, pars_irn);
  time_l1(TVL1_WTHRESH_2,m) = toc(t);


  snr_l1(TVL1_WTHRESH_2,m)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WTHRESH_2,m) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WTHRESH_2,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WTHRESH_2,m) = NaN;
  end

end

% ===================
% >>> WEIGHTS_MIL <<<   -- Matrix Inversion Lemma
% ===================

% ---------------------------
% --- adapt PCG tolerance ---

pars_irn = irntvInputPars('l1tv');

pars_irn.pcgtol_ini   = 1e-4;
pars_irn.loops        = 5;
pars_irn.U0           = InCell;
pars_irn.sbstflg      = 0;
pars_irn.adapt_epsR   = 0;
pars_irn.adapt_epsF   = 0;


pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;


pars_irn.loops = 0;
t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WMIL_1,1) = toc(t);


  snr_l1(TVL1_WMIL_1,1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WMIL_1,1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WMIL_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WMIL_1,1) = NaN;
  end


for m=2:Loops

  pars_irn.loops = m-1;
  t = tic;
    I_L1 = irntv(InCell, {}, lambda, pars_irn);
  time_l1(TVL1_WMIL_1,m) = toc(t);


  snr_l1(TVL1_WMIL_1,m)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WMIL_1,m) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WMIL_1,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WMIL_1,m) = NaN;
  end



end

% ---------------------------
% --- fixed PCG tolerance ---


pars_irn.adaptPCGtol  = 0;
pars_irn.pcgtol_ini   = 1e-1;   % fixed PCG tolerance

pars_irn.loops = 0;
t = tic;
  I_L1 = irntv(InCell, {}, lambda, pars_irn);
time_l1(TVL1_WMIL_2,1) = toc(t);


  snr_l1(TVL1_WMIL_2,1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WMIL_2,1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WMIL_2,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WMIL_2,1) = NaN;
  end


for m =1:Loops

  pars_irn.loops = m-1;
  t = tic;
    I_L1 = irntv(InCell, {}, lambda, pars_irn);
  time_l1(TVL1_WMIL_2,m) = toc(t);


  snr_l1(TVL1_WMIL_2,m)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_WMIL_2,m) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_WMIL_2,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_WMIL_2,m) = NaN;
  end


end


% ===================
% >>>    SB-TV    <<<   -- Split-Breagman / ADMM
% ===================

% ----------------------------------
% --- fixed gamma (twice lambda) ---

%  sbOpts.MaxMainIter = pars_irn.loops;
%  sbOpts.BoundRange = [-Inf Inf];
%  sbOpts.gamma = 2*lambda;

sbOpts.MaxMainIter = 0;
t = tic;
%    I_L1 = sbl1tv2_BW(InCell, {}, lambda, sbOpts);
  I_L1 = sbl1tv2_BW(InCell, {}, lambda);
time_l1(TVL1_SB_1,1) = toc(t);

  snr_l1(TVL1_SB_1,1)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_SB_1,1) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_SB_1,1) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_SB_1,1) = NaN;
  end


for m=2:Loops

  sbOpts.MaxMainIter = m-1;
  t = tic;
%    I_L1 = sbl1tv2_BW(InCell, {}, lambda, sbOpts);
  I_L1 = sbl1tv2_BW(InCell, {}, lambda);
  time_l1(TVL1_SB_1,m) = toc(t);

  snr_l1(TVL1_SB_1,m)  = snr(Inoiseless, I_L1);
  psnr_l1(TVL1_SB_1,m) = psnr(Inoiseless, I_L1);

  if(SSIM_CODE)
    ssim_l1(TVL1_SB_1,m) = ssim_index(255*Normalize(Inoiseless), 255*Normalize(I_L1));
  else
    ssim_l1(TVL1_SB_1,m) = NaN;
  end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

function[y] = Normalize(x)

 y = (x - min(x(:)))/(max(x(:)) - min(x(:)));
