function[S IRN_Ig Adapt_Ig] = test_adaptL1(Img, noisy)


if nargin < 2,
  noisy = 3;    % 50 percent SnP noise
  if nargin < 1
    Img='goldhl.';
  end
end


SHOW_IMGS = true;
loops = 2;

nmpdef;

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

      Ig_01L1 = imnoise(Ig, 'salt & pepper', 0.1);
      Ig_03L1 = imnoise(Ig, 'salt & pepper', 0.3);
      Ig_05L1 = imnoise(Ig, 'salt & pepper', 0.5);
      Ig_07L1 = imnoise(Ig, 'salt & pepper', 0.7);


      Ig_noisy = {Ig_01L1; Ig_03L1; Ig_05L1; Ig_07L1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       L1 Denoise        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_irn    = 1.2;

% pre-processing

S0 = adaptMedian_colfilt(Ig_05L1);


lambda_adapt  = 1e-4*(S0 == 0) +  5*(S0 > 0);


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

pars_irn_adapt.loops      = 8;

pars_irn.variant       = NMP_TV_LamdaAdapt;


pars_irn_adapt.U0      = Ig_05L1;   % necessary the for adapt case

% ---------- 
% ---------- 

  tic;
  IRN_Ig = irntv(Ig_05L1, {}, lambda_irn, pars_irn);
  tirn_denoise = toc;

  irn_snr_denoise = snr(Ig, IRN_Ig);



  tic;
  Adapt_Ig = irntv(Ig_05L1, {}, lambda_adapt, pars_irn_adapt);
  tirn_denoise = toc;

  adapt_snr_denoise = snr(Ig, Adapt_Ig);


  if(SHOW_IMGS)
    snr_noisy = snr(Ig_noisy{noisy}, IRN_Ig);

    figure; imagesc( Normalize(Ig_noisy{noisy}) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Original Image - SNR: %4.1fdB.\n ', ...
                  irn_snr_denoise));

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
%  (not really needed; here for "historical" reasons)

for m = 2:loops

  Ipass = Adapt_Ig;

  S = adaptMedian_colfilt(Ipass);

  D1 = ( S0 - (S0 > 0).*(S0 - S) ) > 0;
  D2 = (S0>0).*(D1==0);

  lambda_adapt = lambda_adapt.*(1.1*D1) + lambda_adapt.*(0.9*D2);

  pars_irn_adapt.U0      = Ipass;   % necessary the for adapt case

  tic;
  Adapt_Ig = irntv(Ipass, {}, lambda_adapt, pars_irn_adapt);
  tirn_denoise = toc;

  adapt_snr_denoise = snr(Ig, Adapt_Ig);


  if(SHOW_IMGS)

    figure; imagesc( Normalize(Adapt_Ig) ); 
    axis image; axis off; colormap gray;
    title(sprintf('Denoised Image - IRN Adapt(2nd). SNR: %4.1fdB.\n ', ...
                  adapt_snr_denoise));

  end

  S0 = S;

end

%==========================================================

function[S] = adaptMedian(I)

[Nrows Ncols] = size(I);

S = zeros(Nrows, Ncols);

vmin = 0;
vmax = 255;

w = 1;
w_max = 9;

for n = 1:Ncols,
  
  for k = 1:Nrows,

    cond = 1;

    while(cond)

    [p_min, p_max, p_med] = compute(I, w, k, n, Nrows, Ncols);

    if( (p_min < p_med) && (p_med < p_max) )

        if( ~( (p_min < I(k,n)) && (I(k,n) < p_max) ) ) 
%            I(k,n) = p_med;
          S(k,n) = w+2;
        end

        cond = 0;

    else

      w = w+2;
      if(w>w_max) 
%          I(k,n) = p_med;
        S(n,k) = w_max+2;
        cond = 0;
      end

    end

    end % _END_ WHILE(cond)

    w = 1;

  end % _END_ FOR(k)

end % _END_ FOR(n)



% ----------------------------------
function[p_min, p_max, p_med] = compute(I, w, k, n, Nrows, Ncols)

    k_min = k - w; if(k_min <= 0) k_min = 1; end
    k_max = k + w; if(k_max > Nrows) k_max = Nrows; end

    n_min = n - w; if(n_min <= 0) n_min = 1; end
    n_max = n + w; if(n_max > Ncols) n_max = Ncols; end

    p     = I(k_min:k_max, n_min:n_max);
    p_min = min(p(:));
    p_max = max(p(:));
    p_med = median(p(:));

