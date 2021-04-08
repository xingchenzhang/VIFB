function[str_all, Irec] = icassp12(example, Img, Nloops, SHOW_IMGS);


if nargin < 4
  SHOW_IMGS = false;
  if nargin < 3
    Nloops = 1;
    if nargin < 2
      Img = lena;
      if nargin < 1;
	example = 'gauss_sp_denoise';
      end
    end
  end
end


%  SHOW_IMGS = true;

%  example = 'icassp12'
%  example = 'gauss_sp_denoise';
%  example = 'gauss_random_denoise';


%  Img = 'Clena';
%  Img = 'peppers';
%  Img = 'lena';
%  Img = 'barb';
%  Img = 'cman';
%  Img = 'Cboats';

%  Img = 'all';


%  Nloops = 1; % This value will be change to 10 if example == icassp12



%%%%%%%%%%%%%%%%%%%%%%%%
%        Setup         %
%%%%%%%%%%%%%%%%%%%%%%%%

  % Check for SSIM
  SSIM_CODE = check4ssim();

  % Setup simulation
  [Nloops, Img] = setupSim(example, Nloops, Img);


  % open string (for reports)
  str_all = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( strcmp(lower(Img), 'all') )

  all = 1;
  if( strcmp(example, 'icassp12') )
    Img = {'lena', 'goldhl', 'cman512', 'barb', 'Clena', 'Cbarb', 'Cboats', 'Cgoldhl' };
  else
    Img = {'lena', 'peppers', 'goldhl.', 'cman512', 'Clena', 'Cpeppers', 'Cmandrill'};
  end

  NImgs = length(Img);


else

  all = 0;
  NImgs = 1;


end  % _END_ IF(strcmp(lower(Img))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          loop over Img
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t = 1:NImgs,

  if( all == 1 )
    tImg = Img{t};
  else
    tImg = lower(Img);
  end

  [Ig Color] = selectImage(tImg);

  if Color == 1
    SSIM_CODE = 0;
  end

  disp(' ');
  disp(sprintf('%s Color %d', tImg, Color));
  disp(' ');


  [Nrows Ncols Ndims] = size(Ig);


  % Regularization parameter of the L2 / L1
  [lambdaS_l2 Lambda_O Lambda_update] = setRegParameter(example, Ndims);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  loops for average results    %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for loops=1:Nloops


    %%%%%%%%%%%%%%%%%
    %  add noise    %
    %%%%%%%%%%%%%%%%%

    % Noise level
    sigmaGauss = [5/255, 10/255, 15/255 ];
    NGauss = length(sigmaGauss);
        
    spnoise = [0.3 0.5 0.7];
    Nspnoise = length(spnoise);

    randomnoise = [0.1 0.2 0.3];
    NRnoise = length(randomnoise);

    % noisy images
    clear InCell;
    [InCell noiseName vNoise nNoise] = addNoise(Ig, example, sigmaGauss, NGauss, spnoise, Nspnoise, randomnoise, NRnoise);
    nNoisyIm = length(InCell);



    % loop over noisy images
    for k=1:length(InCell)

      In = InCell{k};

%^^^^^% ^^^^^^^^^^^^^^^^^^^^^
%     %   Clock is ticking
%     % ---------------------
      tstart = tic;


      % First estimate the pixels that are corrupted
      switch lower(example)
    
        case{'gauss_sp_denoise'}
          S0 = adaptMedian_colfilt(In);		% estimate S&P corrupted pixels

        case{'gauss_random_denoise'}	% estimate random value (impulse noise) corrupted pixels
          S0 = estimateRValImpulse(In);

      end	% _END_ SWITCH


      % Set lambda (regularization) matrix
      lambda_l1l2 = zeros(Nrows, Ncols, 2*Ndims);

      for d=1:Ndims,
          lambda_l1l2(:,:,d)  =  Lambda_O*(S0(:,:,d) ~= 0);
          lambda_l1l2(:,:,Ndims+d) = lambdaS_l2(k)*(S0(:,:,d)==0);
      end


      % =================================================
      % ---------- Setup for IRN_L1L2  ------------------

      pars_irn_adapt = irntvInputPars('l1tv_l2tv');

      pars_irn_adapt.adapt_epsR   = 1;
      pars_irn_adapt.epsR_cutoff  = 0.01;
      pars_irn_adapt.adapt_epsF   = 1;
      pars_irn_adapt.epsF_cutoff  = 0.05;

      pars_irn_adapt.pcgtol_ini = 1e-4;

      switch lower(example)    
        case{'gauss_sp_denoise'}
          pars_irn_adapt.loops      = 1;
          IRNloops = 7;
        case{'gauss_random_denoise'}
          pars_irn_adapt.loops      = 2;
          IRNloops = 4;
      end

      pars_irn_adapt.U0      = In;   % necessary the for adapt case

      % =================================================


      % 1st call
      Ig_L1L2 = irntv(In, {}, lambda_l1l2, pars_irn_adapt);


      % loop
  
      for m = 2:IRNloops
   
        for d=1:Ndims
          lambda_l1l2(:,:,d)  =  Lambda_update*lambda_l1l2(:,:,d);
        end

        pars_irn_adapt.U0      = Ig_L1L2;	% reset initial sol. to previous sol.

        Ig_L1L2 = irntv(In, {}, lambda_l1l2, pars_irn_adapt);

      end % _END_ FOR(m)

      % denoised image
      Irec{k} = Ig_L1L2;



%^^^^^% ^^^^^^^^^^^^^^^^^^^^^
%     %    Time performance
      l1l2_time(k, loops) = toc(tstart);
%     % ---------------------



      if(SHOW_IMGS)
        figure; imagesc( Normalize(Irec{k}) ); axis image; axis off; colormap gray; 
      end


      % Reconstruction quality assessment
      [psnr_l1l2(k, loops), snr_l1l2(k, loops), ssim_l1l2(k, loops)] = qualityAssessment(Ig, Irec{k}, SSIM_CODE);
      



      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %           Generate reports

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      if(k==1)

        if(Color == 0) str = sprintf('Img: grayscale %s (denoising)\n', lower(tImg));
        else str = sprintf('Img: color %s (denoising)\n', lower(tImg));
        end

        str_all = [ str_all str ];

        str = sprintf(' sigma/%s \t SNR \t SSIM \t PSNR    Time\n', noiseName);
        str_all = [ str_all str ];

      end

        str = sprintf(' %1.0f/255 -- %1.1f \t %1.2f \t %1.2f \t %1.2f    %1.2f\n', ...
                   255*sigmaGauss( rem(k,NGauss) + NGauss*(rem(k,NGauss)==0) ), vNoise( floor((k-1)/nNoise)+1 ), ...
                   snr_l1l2(k, loops), ssim_l1l2(k, loops), psnr_l1l2(k, loops), l1l2_time(k, loops) );

%        if(k<=nNoisyIm)
%          str = sprintf(' %1.0f/255 -- %1.1f \t %1.2f \t %1.2f \t %1.2f    %1.2f\n', ...
%                     255*sigmaGauss( rem(k,NGauss) + NGauss*(rem(k,NGauss)==0) ), vNoise( floor((k-1)/nNoise)+1 ), ...
%                     snr_l1l2(k, loops), ssim_l1l2(k, loops), psnr_l1l2(k, loops), l1l2_time(k, loops) );
%        else
%          str = sprintf(' %1.0f/255 -- %1.1f \t %1.2f \t %1.2f \t %1.2f    %1.2f\n', ...
%                     255*sigmaGauss( rem(k,NGauss) + NGauss*(rem(k,NGauss)==0) ), vNoise( floor((k-1-nNoisyIm)/nNoise)+1 ), ...
%                     snr_l1l2(k, loops), ssim_l1l2(k, loops), psnr_l1l2(k, loops), l1l2_time(k, loops) );
%        end
      str_all = [ str_all str ];


      if(k==length(InCell)) 
        str = sprintf(' \n\n');
      end


    end % _END_ FOR(length(InCell))


  end % FOR(loops)

  disp(str_all)

end % FOR(t)  -- test Image


end


  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[SSIM_CODE] = check4ssim()

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


end


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function[Nloops, Img] = setupSim(example, Nloops, Img)

  tmpN = Nloops;
  tmpI = Img;

  if( strcmp(example, 'icassp12') )

    disp('  Running the "icassp11" example...');
    disp('  This script reproduces the Exprimental Results of ');
    disp('  "Mixed Gaussian-impulse Noise Image Restoration via');
    disp('    Total Variation"');
    disp(' ')
    disp('  This script takes sometime (ten-trial average)');
    disp(' ')
    tmpN = 10;

    if( ~strcmp(tmpI, 'all') ) 

      disp('NOTE:');
      disp('  If you run the "icassp12" example, you should set Img = "all" ...');
      disp('  setting Img = "all" ');
      tmpI = 'all';
    end

  end

  Nloops = tmpN;
  Img    = tmpI;

end


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%


function[Ig Color] = selectImage(tImg)

  switch lower(tImg)

    case{'lena'}
      Ig = double( imread('gray_imgs/lena_gray_512.png') ) / 255;
      Color = 0;

    case{'peppers'}
      Ig = double( imread('gray_imgs/peppers_gray.png') ) / 255;
      Color = 0;

    case{'cman'}
      Ig = double( imread('gray_imgs/cameraman_256x256.tiff') ) / 255;
      Color = 0;

    case{'cman512'}
      Ig = double( imread('gray_imgs/cameraman_512x512.tiff') ) / 255;
      Color = 0;

    case{'goldhl'}
      Ig = double( imread('gray_imgs/goldhill_gray.png') ) / 255;
      Color = 0;

    case{'barb'}
      Ig = double( imread('gray_imgs/barbara_gray512.png') ) / 255;
      Color = 0;

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
      disp('Exiting ICASSP12 test code...');
      return;

    otherwise
      error('Not a valid image\n');

  end  % _END_ SWITCH(Img)


end




%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function[lambdaS_l2, Lambda_O, Lambda_update] = setRegParameter(example, Ndims)

    switch lower(example)
    
    case{'gauss_sp_denoise'}
      if(Ndims == 1)
        lambdaS_l2 = [ 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025 ];
      else
        lambdaS_l2 = 1.5*[ 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025 ];
      end

      Lambda_O = 5;
      Lambda_update = 1.15;

    case{'gauss_random_denoise'}
      if(Ndims == 1)
        lambdaS_l2 = [ 1.25*0.0075, 0.0185, 0.025, ...
                 1.25*0.0075, 0.0185, 0.025, ...
                 1.25*0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025 ];
%                   0.0075, 0.0185, 0.025 ]/8;
      else
        lambdaS_l2 = 1.5*[ 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025, ...
                 0.0075, 0.0185, 0.025 ]/8;
      end


      Lambda_O = 5;
      Lambda_update = 1.15;

    end

end


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function[InCell, noiseName, vNoise, nNoise] = addNoise(Ig, example, sigmaGauss, NGauss, spnoise, Nspnoise, randomnoise, NRnoise)

    % Add Gaussian noise
    Ig_05L2 = imnoise(Ig,'gaussian', 0, (sigmaGauss(1)^2) );
    Ig_10L2 = imnoise(Ig,'gaussian', 0, (sigmaGauss(2)^2) );
    Ig_15L2 = imnoise(Ig,'gaussian', 0, (sigmaGauss(3)^2) );
      
    % Add impulse noise (S&P or random value)
    switch lower(example)
    
      case{'gauss_sp_denoise'}
	
	noiseName = 'S&P';
	vNoise    = spnoise;
	nNoise    = Nspnoise;

	for k=1:Nspnoise,
      
	  InCell{NGauss*(k-1)+1} = imnoise(Ig_05L2, 'salt & pepper', spnoise(k));
	  InCell{NGauss*(k-1)+2} = imnoise(Ig_10L2, 'salt & pepper', spnoise(k));
	  InCell{NGauss*(k-1)+3} = imnoise(Ig_15L2, 'salt & pepper', spnoise(k));    
	end


      case{'gauss_random_denoise'}

	noiseName = 'RandVal';
	vNoise    = randomnoise;
	nNoise    = NRnoise;

	[Nrows Ncols Ndims] = size(Ig);
	dummy = 0.5*ones(Nrows, Ncols, Ndims);

	for k=1:NRnoise,

	  mask = (imnoise(dummy, 'salt & pepper', randomnoise(k)) ~= 0.5);
	  InCell{NGauss*(k-1)+1} = Ig_05L2.*(1-mask) + mask.*rand(Nrows, Ncols, Ndims);

	  mask = (imnoise(dummy, 'salt & pepper', randomnoise(k)) ~= 0.5);
	  InCell{NGauss*(k-1)+2} = Ig_10L2.*(1-mask) + mask.*rand(Nrows, Ncols, Ndims);

	  mask = (imnoise(dummy, 'salt & pepper', randomnoise(k)) ~= 0.5);
	  InCell{NGauss*(k-1)+3} = Ig_15L2.*(1-mask) + mask.*rand(Nrows, Ncols, Ndims);
	end


    end % _END_ SWITCH

end


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function[S0] = estimateRValImpulse(In)

          thres_diffmed= 0.2;
          window_side= 2;

          road_size= 7;
          thres_ranked= 0.05;

          [Nrows Ncols Ndims] = size(In);

          if(Ndims == 1)

            mask_diff_median = diff_median_blk( In, window_side, thres_diffmed);
            mask_ranked_mean = ranked_mean( In, window_side, road_size, thres_ranked);
            S0 = mask_ranked_mean & mask_diff_median;

          else
            
            for d=1:Ndims,
              mask_diff_median = diff_median_blk( In(:,:,d), window_side, thres_diffmed);
              mask_ranked_mean = ranked_mean( In(:,:,d), window_side, road_size, thres_ranked);
              S0(:,:,d) = mask_ranked_mean & mask_diff_median;
            end

          end

end


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function[psnr_l1l2, snr_l1l2, ssim_l1l2] = qualityAssessment(Ig, Ig_L1L2, SSIM_CODE)
      

      psnr_l1l2 = psnr(Ig, Ig_L1L2);
      snr_l1l2  = snr(Ig, Ig_L1L2); 
      if(SSIM_CODE)
        ssim_l1l2 = ssim_index(255*Normalize(Ig), 255*Normalize(Ig_L1L2));
      else
        ssim_l1l2 = NaN;
      end

end

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function[y] = Normalize(x)

  y = (x - min(x(:)))/(max(x(:)) - min(x(:)));

end




