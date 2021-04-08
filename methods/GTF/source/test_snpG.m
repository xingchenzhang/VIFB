


HEURISTICS = 0;


%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%

Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input images        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Img = 'lena';
Img = 'clena';

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
  %       noisy images        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      Ig_01L1 = imnoise(Ig, 'salt & pepper', 0.1);
      Ig_02L1 = imnoise(Ig, 'salt & pepper', 0.2);
      Ig_03L1 = imnoise(Ig, 'salt & pepper', 0.3);
      Ig_05L1 = imnoise(Ig, 'salt & pepper', 0.5);
      Ig_07L1 = imnoise(Ig, 'salt & pepper', 0.7);
      Ig_09L1 = imnoise(Ig, 'salt & pepper', 0.9);



      Ig_05L2 = imnoise(Ig,'gaussian', 0, ((5/255)^2) );
      Ig_10L2 = imnoise(Ig,'gaussian', 0, ((10/255)^2) );
      Ig_20L2 = imnoise(Ig,'gaussian', 0, ((20/255)^2) );
      Ig_30L2 = imnoise(Ig,'gaussian', 0, ((30/255)^2) );

      InCell{1} = imnoise(Ig_10L2, 'salt & pepper', 0.2);
      InCell{2} = imnoise(Ig_20L2, 'salt & pepper', 0.2);
      InCell{3} = imnoise(Ig_30L2, 'salt & pepper', 0.2);
      InCell{4} = imnoise(Ig_10L2, 'salt & pepper', 0.4);

      % -----------

      InCell{5} = imnoise(Ig_05L2, 'salt & pepper', 0.3);
      InCell{6} = imnoise(Ig_05L2, 'salt & pepper', 0.5);

      InCell{7} = imnoise(Ig_10L2, 'salt & pepper', 0.3);
      InCell{8} = imnoise(Ig_10L2, 'salt & pepper', 0.5);

%        snr_noisy12 = snr(Ig, In);

%      figure; imagesc( Normalize(In) ); 
%      axis image; axis off; colormap gray;
%      title(sprintf('Noisy - L1/L2. SNR: %4.1fdB.\n ', ...
%                    snr_noisy12));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing

InputDims = size(Ig);

if(length(InputDims) == 2)
  lambdaS_l2 = [ 0.0185, 0.03, 0.045, 0.0185, 0.0075, 0.0075, 0.0185, 0.0185 ];
else
  lambdaS_l2 = 1.5*[ 0.0185, 0.03, 0.045, 0.0185, 0.0075, 0.0075, 0.0185, 0.0185 ];
end


for k=1:length(InCell)

In = InCell{k};

  tstart = tic;

S0 = adaptMedian(In);

if(length(InputDims) == 2)

  if(HEURISTICS)
    lambda_l1l2(:,:,1)  =  0.85*(S0 == 3) + ...
                    1.1*(S0 == 5) + 1.2*(S0 == 7) + 1.4*(S0 == 9);
  else
    lambda_l1l2(:,:,1)  = 1.5*(S0 ~= 0);
  end

  lambda_l1l2(:,:,2) = lambdaS_l2(k)*(S0==0);

  planes=1;
else

  if(HEURISTICS)

    for d=1:InputDims(3),
    lambda_l1l2(:,:,d)  =  1.05*(S0(:,:,d) == 3) + ...
                    1.3*(S0(:,:,d) == 5) + 1.4*(S0(:,:,d) == 7) + 1.6*(S0(:,:,d) == 9);
    end
  else

    for d=1:InputDims(3),
      lambda_l1l2(:,:,d)  =  1.75*(S0(:,:,d) ~= 0);
    end
  end

  for d=1:InputDims(3),
    lambda_l1l2(:,:,3+d) = lambdaS_l2(k)*(S0(:,:,d)==0);
  end

  planes=3;
end

% ---------- Setup for IRN_L1L2  ------------------

pars_irn_adapt = irntvInputPars('l1tv_l2tv');

pars_irn_adapt.adapt_epsR   = 1;
pars_irn_adapt.epsR_cutoff  = 0.01;
pars_irn_adapt.adapt_epsF   = 1;
pars_irn_adapt.epsF_cutoff  = 0.05;

pars_irn_adapt.pcgtol_ini = 1e-4;

pars_irn_adapt.loops      = 1;

pars_irn_adapt.U0      = In;   % necessary the for adapt case
%  pars_irn_adapt.U0      = lmu;   % necessary the for adapt case


  Ig_L1L2 = irntv(In, {}, lambda_l1l2, pars_irn_adapt);


loops = 6;
%  
for m = 2:loops
   
   for d=1:InputDims(3)
      lambda_l1l2(:,:,3+d) = (0.975^(m-1))*lambda_l1l2(:,:,3+d);
   end
   Ipass = Ig_L1L2; 
   dI = Ipass - In; 
   lambda_l1l2(:,:,1:planes) = adaptLambda(dI,  lambda_l1l2(:,:,1:planes), 0.5, S0, 0.8); 
   pars_irn_adapt.U0      = Ipass;
   Ig_L1L2 = irntv(In, {}, lambda_l1l2, pars_irn_adapt);

end
l1l2_time(k) = toc(tstart)
%  
figure; imagesc( Normalize(Ig_L1L2) ); axis image; axis off; colormap gray; 

psnr_l1l2(k) = psnr(Ig, Ig_L1L2)


end

psnr_bing = [32.84, 30.23, 27.33, 29.96]; 
psnr_xiao = [36.20, 33.93, 33.19, 31.51];

[ psnr_l1l2' [psnr_bing'; psnr_xiao']]