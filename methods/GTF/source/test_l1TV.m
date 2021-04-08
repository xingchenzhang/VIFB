

nmpdef;


I = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

In = imnoise(I, 'salt & pepper', 0.2);

figure; imagesc(In); colormap gray; axis off; axis image

% --------------------------------------------
% --- Adapt cutoff value plus substitution ---

pars_irn = irntvInputPars('l1tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 5;
pars_irn.U0         = In;

pars_irn.variant       = NMP_TV_SUBSTITUTION;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;

t = tic;
I_Threshold = irntv(In, {}, 1.25, pars_irn);
toc(t)

snr(I, I_Threshold)


% -----------------------------------------------
% --- Adapt cutoff value without substitution ---

pars_irn.variant       = NMP_TV_STANDARD;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;
pars_irn.U0            = [];

t = tic;
I_Threshold = irntv(In, {}, 1.25, pars_irn);
toc(t)

snr(I, I_Threshold)


% --------------------------
% --- Fixed cutoff value ---


pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;



% ------------------------
% ---   Add epsilon    ---


pars_irn.weight_scheme = NMP_WEIGHTS_EPSILON;

% ------------------------------
% --- Matrix inversion lemma ---


pars_irn = irntvInputPars('l1tv');

pars_irn.pcgtol_ini   = 1e-4;
pars_irn.loops        = 5;
pars_irn.U0           = In;
pars_irn.sbstflg      = 0;
pars_irn.adapt_epsR   = 0;
pars_irn.adapt_epsF   = 0;


pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;
pars_irn.pcgtol_ini    = 1e-2;
pars_irn.adaptPCGtol   = 1;


t = tic;
I_mil = irntv(In, {}, 1.25, pars_irn);
toc(t)

figure; imagesc(I_mil.*(I_mil>=0).*(I_mil<=1)); colormap gray; axis off; axis image
%  figure; imagesc(I_mil); colormap gray; axis off; axis image

snr(I, I_mil)


