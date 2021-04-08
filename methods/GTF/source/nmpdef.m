

% TV problem
% -----------------------------

NMP_L1TV                = 1;    % L1-TV (e.g. salt & pepper noise).
NMP_L2TV                = 2;    % L2-TV (e.g. Gaussian noise).
NMP_LPTV                = 3;    % 
NMP_L1L2TV              = 4;    % Mixed-noise TV (impulsive plus Gaussian noise).
NMP_TV_POISSON          = 5;    % Poisson noise (sol. is non-negative).
NMP_TV_SPECKLE          = 6;    % Speckel (Gamma) noise (sol. is non-negative).



% TV problem variant
% -----------------------------

NMP_TV_STANDARD         = 10;
NMP_TV_SUBSTITUTION     = 11;   % substitution for lp-TV (menaningful for p ~= 2)
NMP_TV_NQP              = 12;   % non-negative solution
NMP_TV_MIL              = 13;   % matrix inversion lemma

NMP_TV_LamdaAdapt       = 14;

NMP_TV_SPECKLE_AA       = 15;   % Aubert-Aujol        F(u) = \sum{  u_k/b_k + \log u_k  }
NMP_TV_SPECKLE_SO       = 16;   % Shi-Osher           z = log(u), F(z) = \sum{  z_k + b_k \exp -z_k  }
NMP_TV_SPECKLE_SO_FULL  = 17;   % Shi-Osher "full"    z = log(u), a1, a2 ctes.
                                %                     F(z) = \sum{  a1*b_k \exp -z_k + (a2/2)*b_k^2 \exp(-2 z_k)+ (a1+a2)*z_k  }

NMP_TV_SB               = 20;   % Split-bregman

NMP_TV_STANDARD_LamdaAdapt = 21;
NMP_TV_L1L2 = 22;
NMP_TV_CONSTRAINED = 23;
NMP_TV_POISSON_LamdaAdapt = 24;

% TV weighting scheme
% -----------------------------

NMP_WEIGHTS_THRESHOLD   = 30;   %      || grad u ||_1 = 0.5 || WDu ||_2^2,   W = 1/sqrt(| grad u |)   if  | grad u | > eps
                                %                                              = 1/sqrt( eps )        if  | grad u | <= eps

NMP_WEIGHTS_EPSILON     = 31;   %      || grad u ||_1 = 0.5 || sqrt( |grad u| + eps ) ||_2^2

NMP_WEIGHTS_HUBER       = 32;   %      || grad u ||_1 = Huber( | grad u |^2, eps )

NMP_WEIGHTS_MIL         = 33;   %       Related to NMP_TV_MIL variant (no thresholding)
