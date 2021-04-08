function [U, itstat] = irntv(S, KC, lambda, pars)
%  
% function [U, itstat] = irntv(S, KC, lambda, parameters)
%
% irntv --  Compute the minimum of a generalised TV functional 
%
%           T = F( S | K*U ) + lambda*|| sqrt( (Dx(U))^2 + (Dy(U))^2 ) ||^q
%               
%           where F( S | K*U ) depends on the noise model (eg. Guassian --> || K*U - S ||_2^2,
%           Salt & Pepper --> || K*U - S ||_1^1, Poisson --> sum( (K*U) - S.*log( (K*U) ) ),
%           etc. for grayscale / color (vector) images using the IRN [1,2,3,4,5,6,7] algorithm,
%           and q in (0 2].
%  
%
% Usage:
%       [U, itstat] = irntv(S, KC, lambda, parameters)
%
% Input:
%       S           Input image
%       KC          Cell array {K, KT} of function handles for
%                   forward linear operator and its
%                   transpose. These operators act on images rather
%                   than vectors, but the transpose should be
%                   considered to be in the sense of the operator
%                   on a vectorised image. Set to [] for a
%                   denoising problem.
%       lambda      Regularisation term weighting factor.
%       parameters  see the irntvInputPars.m file (>>help irntvInputPars.m)
%                   for setting p, q, etc. and/or to choose the enforcing
%                   of the non-negativity constrain.
%
% Output:
%       U           Output image
%       itstat      Iteration status
%  
%  
%
%           [1] "Efficient Minimization Method for a Generalized Total 
%                Variation Functional"
%               IEEE Transactions on Image Processing, 2009, 18:2(322-332).
%
%           [2] "A Generalized Vector-Valued Total Variation Algorithm"
%               Proceedings of the 2009 International Conference on
%               Image Processing (ICIP 2009), pp. 1309-1312, 2009
%  
%           [3] "A Non-negative Quadratic Programming approach to Minimize  
%               the Generalized Vector-Valued Total Variation Functional"
%               Proceedings of the 2010 European Signal Processing Conference
%                (EUSIPCO 2010), pp. 314-318, 2010
%  
%           [4] "Multiplicative Updates Algorithm To Minimize The Generalized 
%               Total Variation Functional With A Non-Negativity Constraint"
%               Proceedings of the 2010 International Conference on
%               Image Processing (ICIP 2010), pp. 2509-2512, 2010
%  
%           [5] "Spatially Adaptive Total Variation Image Denoising Under Salt 
%               And Pepper Noise", Proceedings of the European Signal Processing 
%               Conference (EUSIPCO 2011), pp 278-282, August, 2011
%  
%           [6] "Total Variation Regularization for Poisson Vector-Valued Image 
%               Restoration with a Spatially Adaptive Regularization Parameter
%               Selection",  Proceedings of the IEEE Symposium on Image and Signal 
%               Processing and Analysis (ISPA), pp. 402-407, September, 2011.
%  
%           [7] "Mixed Gaussian-Impulse Noise Image Restoration via Total Variation"
%               IEEE International Conference on Acoustics, Speech, and Signal 
%               Processing (ICASSP), pp. 1077--1080, March, 2012 
%  
%           irntv.m is part of NUMIPAD (http://numipad.sf.net). 
%
%           The NUMIPAD library is being developed under U.S. Government
%           contract W-7405-ENG-36 for Los Alamos National Laboratory.
%
%
%
% Authors
%   Paul Rodriguez    prodrig@pucp.edu.pe
%   Brendt Wohlberg   brendt@tmail.lanl.gov
%  



  % ---------------------------------
  % Simple check por input parameters
  % ---------------------------------

  if nargin < 4,
    pars = irntvInputPars('l1tv');
    disp('Model parameters has not been given... choosing the l1-TV model');
    if nargin < 3,
      lambda = 1;
      if nargin < 2,
        KC = [];
      end
    end
  end

  pars.lambda = lambda;


  % Define constants (see nmpdef.m)
  nmpdef;


  % Check Input/Output general parameters options
  [dflag, InputDims, pars.sbstflg, pars.lambda_ini] = check_irntv_options(nargout, S, KC, pars);
  

  % Set up status display
  iterdispopen(dflag);      


  % choose Regularization weighting scheme (valid for all cases)
  [fR fR_var4] = choose_Regularization_weighting_scheme(pars);


  % choose Fidelity weighting scheme (no need for POISSON or SPECKLE cases)
  if( ~((pars.problem == NMP_TV_POISSON) || (pars.problem == NMP_TV_SPECKLE)) )
      [fF fF_var3] = choose_Fidelity_weighting_scheme(pars);
  end


  %-----------------------------------------------------------------------
  %-----------------------------------------------------------------------
  % Choose between different TV problems

  switch( pars.problem )

    case {NMP_L1TV, NMP_L2TV, NMP_LPTV}

      % ------
      switch( pars.variant )

        case {NMP_TV_NQP}
          % call irntv_nqp
          [U, itstat] = irntv_nqp(S, KC, lambda, pars, InputDims, fF, fF_var3, fR, fR_var4, dflag);
      
        case {NMP_TV_STANDARD, NMP_TV_SUBSTITUTION, NMP_TV_MIL, NMP_TV_LamdaAdapt}
          % Call irntv_standard
          [U, itstat] = irntv_standard(S, KC, lambda, pars, InputDims, fF, fF_var3, fR, fR_var4, dflag);

        case {NMP_TV_SB}
          % Call irntv_sb
          error('to be implemented: Split-Bregma / ADMM TV.');

        otherwise
          error('Unknown L1-TV / L2-TV variant.');

      end % _END_SWITCH( variant )
      % ------

      %  -------------------------------------------------------------------------------------------

    case {NMP_L1L2TV}

      % Call irntv_l1l2
      [U, itstat] = irntv_l1l2(S, KC, lambda, pars, InputDims, fF, fF_var3, fR, fR_var4, dflag);

      %  -------------------------------------------------------------------------------------------

    case {NMP_TV_POISSON}

      % call irntv_poisson
      [U, itstat] = irntv_poisson(S, KC, lambda, pars, InputDims, fR, fR_var4, dflag);


    case {NMP_TV_SPECKLE}

      % call irntv_speckle
      [U, itstat] = irntv_speckle(S, KC, lambda, pars, InputDims, fR, fR_var4, dflag);

    otherwise
      error('Unknown TV problem.');

  end % _END_SWITCH( problem )



  %-----------------------------------------------------------------------

return



% ===========================================================================
% ===========================================================================

function[dflag, InputDims, sbstflg, lambda_ini] = check_irntv_options(irntv_Nargs, S, KC, pars)

  % Define constants (see nmpdef.m)
  nmpdef;

  dflag = 1;
  sbstflg = 1;
  lambda_ini = pars.lambda_ini;

  % Check number of outputs
  if irntv_Nargs > 1,
    dflag = 0;
  end

  % Check wether the input image is grayscale or color
  InputDims = size(S);

  if( (length(InputDims) < 2) && (length(InputDims) > 3) )
    disp( sprintf('\nInput data has %d dimensions...\n', length(InputDims) ) );
    error('Not a valid grayscale/color image');
  end


  % Check argument KC
  if ~isempty(KC),
    if ~(iscell(KC) && isa(KC{1}, 'function_handle') && isa(KC{2}, 'function_handle')),
      error('argument KC must be empty or a cell array {K,KT} of function handles');
    end

    sbstflg = 0; % No sbstflg for deconvolution

  end


  % check for sbstflg 
  if( (pars.problem == NMP_L2TV) || (pars.p==2) ) % Incongruent
    sbstflg = 0;
  end

  if(pars.variant == NMP_TV_STANDARD)
    sbstflg = 0;
  end

  
  % check for lambda_ini
  if( isempty(pars.lambda_ini) ) 
    lambda_ini = pars.lambda; 
  end

  


return

% ===========================================================================
% ===========================================================================

function[fF fF_var3] = choose_Fidelity_weighting_scheme(pars)
% Choose the weighting scheme for the Fidelity term (valid for l1/l2-TV)

  nmpdef;


  switch(pars.weight_scheme) 
 
    % Apply threshold to weights
    case{NMP_WEIGHTS_THRESHOLD}

      if( isempty(pars.adapt_epsF) || pars.adapt_epsF == 0 )    % No adaptation of the Fidelity wieghts

        if pars.sbstflg
          fF = @(var1, var2, var3) fF_fixed_substitute(var1, var2, var3);
        else
          fF = @(var1, var2, var3) fF_fixed(var1, var2, var3);
        end
        fF_var3 = pars.epsF;

      else % Adaptation of the Fidelity wieghts

        if pars.sbstflg
          fF = @(var1, var2, var3) fF_adapt_substitute(var1, var2, var3);
        else
          fF = @(var1, var2, var3) fF_adapt(var1, var2, var3);
        end
        fF_var3 = pars.epsF_cutoff;

      end


%      % Use  0.5 || sqrt( |u| + eps ) ||_2^2
%      case{NMP_WEIGHTS_EPSILON}
%        
%        disp('Fidelity NMP_WEIGHTS_EPSILON to be implemented');
%  
%  
%      % Use the Huber function
%      case{NMP_WEIGHTS_HUBER}
%  
%        disp('NMP_WEIGHTS_HUBER to be implemented');


    % Use the matrix inversion lemma to avoid thresholding
    case{NMP_WEIGHTS_MIL}

      fF = @(var1, var2, var3) fF_MatInvLemma(var1, var2, var3);
      fF_var3 = [];


  end


return

% ===========================================================================
% ===========================================================================


function[fR fR_var4] = choose_Regularization_weighting_scheme(pars)
% Choose the weighting scheme for the Regularization term (valid for all schemes)

  nmpdef;

  switch(pars.weight_scheme) 
 
    % Apply threshold to weights
    case{NMP_WEIGHTS_THRESHOLD}
      if( isempty(pars.adapt_epsR) || pars.adapt_epsR == 0 )    % Regularization
        fR = @(var1, var2, var3, var4) fR_fixed( var1, var2, var3, var4);
        fR_var4 = pars.epsR;
      else
        fR = @(var1, var2, var3, var4) fR_adapt( var1, var2, var3, var4);
        fR_var4 = pars.epsR_cutoff;
      end


    % Add epsilon to weights
    case{NMP_WEIGHTS_EPSILON}

      fR = @(var1, var2, var3, var4) fR_AddEPS( var1, var2, var3, var4);
      fR_var4 = pars.epsR;


    case{NMP_WEIGHTS_HUBER}

      disp('Regularization NMP_WEIGHTS_HUBER to be implemented');


    % Use the matrix inversion lemma to avoid thresholding
    case{NMP_WEIGHTS_MIL}

      fR = @(var1, var2, var3, var4) fR_MatInvLemma( var1, var2, var3, var4);
      fR_var4 = pars.lambda;


  end % _END_ SWITCH

return

% ===========================================================================
% ===========================================================================

function [U, itstat] = irntv_standard(S, KC, lambda, pars, InputDims, fF, fF_var3, fR, fR_var4, dflag)

  nmpdef;

  switch(pars.variant)

    case{NMP_TV_STANDARD, NMP_TV_SUBSTITUTION, NMP_TV_MIL} % FIXME: MIL should be here
      if isempty(KC)
        f_IDTD     = @(var1, var2, var3) IDTD(var1, var2, var3);
        f_IDTWD    = @(var1, var2, var3, var4, var5) IDTWD(var1, var2, var3, var4, var5);
        f_IWDTWDW  = @(var1, var2, var3, var4, var5) IWDTWDW(var1, var2, var3, var4, var5);
        f_KTKDTD   = [];
        f_KTWKDTWD = [];
      else
        f_IDTD     = [];
        f_IDTWD    = [];
        f_IWDTWDW  = [];
        f_KTKDTD   = @(var1, var2, var3, var4, var5) KTKDTD(var1, var2, var3, var4, var5); 
        f_KTWKDTWD = @(var1, var2, var3, var4, var5, var6, var7) KTWKDTWD(var1, var2, var3, var4, var5, var6, var7);
      end
  %  

    case{NMP_TV_LamdaAdapt}
      if isempty(KC)
        f_IDTD = @(var1, var2, var3) IDTD_lambdaMtx(var1, var2, var3);
        f_IDTWD = @(var1, var2, var3, var4, var5) IDTWD_lambdaMtx(var1, var2, var3, var4, var5);
        f_IWDTWDW = @(var1, var2, var3, var4, var5) IWDTWDW_lambdaMtx(var1, var2, var3, var4, var5);
        f_KTKDTD   = [];
        f_KTWKDTWD = [];
      else
        f_IDTD     = [];
        f_IDTWD    = [];
        f_IWDTWDW  = [];
        %        f_KTKDTD = @(var1, var2, var3) KTKDTD_lambdaMtx(var1, var2, var3, var4, var5); 
        %        f_KTWKDTWD = @(var1, var2, var3, var4, var5, var6, var7) KTWKDTWD_lambdaMtx(var1, var2, var3, var4, var5, var6, var7);
        f_KTKDTD = @(var1, var2, var3, var4, var5) KTKDTD(var1, var2, var3, var4, var5); 
        f_KTWKDTWD = @(var1, var2, var3, var4, var5, var6, var7) KTWKDTWD(var1, var2, var3, var4, var5, var6, var7);
      end

  end

  fHandles = { fF, fF_var3, fR, fR_var4, f_IDTD, f_IDTWD, f_IWDTWDW, f_KTKDTD, f_KTWKDTWD };


  if isempty(KC),  

    % Denoising problem
    [U, itstat] = irntv_lptv_denoising(S, lambda, pars, InputDims, dflag, fHandles);

  else 

    % General inverse problem with linear operator K
    [U, itstat] = irntv_standard_deconv(S, KC, lambda, pars, InputDims, dflag, fHandles);

  end


  % Final part of status display
  iterdispclose(dflag);

return


% ===========================================================================
% ===========================================================================

function[U, itstat] = irntv_lptv_denoising(S, lambda, pars, InputDims, dflag, fHandles)

  nmpdef;


  % Start clock
  tic;

  % Initialise iteration status vector
  itstat = [];

  % 1D representation of input image
  s = S(:); 


  fF         = fHandles{1};
  fF_var3    = fHandles{2};
  fR         = fHandles{3};
  fR_var4    = fHandles{4};

  f_IDTD     = fHandles{5};
  f_IDTWD    = fHandles{6};
  f_IWDTWDW  = fHandles{7};
  

  if(pars.weight_scheme == NMP_WEIGHTS_MIL)

     wfs = DxDy(S);
     wfs = wfs(:);

  end


  % =====================
  %   Initial Soluction
  % =====================

  if isempty(pars.U0),  % Construct an initial solution
    
    if(pars.weight_scheme ~= NMP_WEIGHTS_MIL)
      [u, ~, ~, pit] = pcg(@(x) f_IDTD(x, InputDims, pars.lambda_ini), s, ...
                             pars.pcgtol_ini, pars.pcgitn, [], [], s);
      U = reshape(u, InputDims);

    else
      [z, flg , ~, pit] = pcg(@(x) IDDT(x, InputDims, pars.lambda_ini), ...
                                wfs, pars.pcgtol_ini, pars.pcgitn, [], [], wfs);
      u = s - DxTDyT(z, InputDims);
      U = reshape(u, InputDims);

    end

  else  % Initial solution is given

    U = pars.U0;
    u = U(:);

  end


  % ===========
  %   Iterate
  % ===========

  for k = 1:pars.loops,

    % Set weight matrices

    if( (k==1) && (sum(abs(U(:)-S(:))) < 0.05/(InputDims(1)*InputDims(2)) ) )
      WF = (2/pars.p)*(0.1)^( (2-pars.p)/2 )*ones(InputDims);
    else
      WF = fF(U-S, pars.p, fF_var3);
    end

    WR = fR(U, InputDims, pars.q, fR_var4);
    

    % Compute relative residual of current solution
    if(pars.weight_scheme ~= NMP_WEIGHTS_MIL)

      if pars.sbstflg, % Use substitution (indirect) method for l1-TV
        wfs = s./vec(WF);
        v = u./vec(WF);
        relres = pcgrelres(@(x) f_IWDTWDW(x, InputDims, WF, WR, lambda), wfs, v);
      else        % Use standard (direct) method for l1-TV
        wfs = vec(WF).*s;
        relres = pcgrelres(@(x) f_IDTWD(x, InputDims, WF, WR, lambda), wfs, u);
      end

    else % NMP_WEIGHTS_MATINVLEMMA


      wfs = DxDy(S);
      wfs = wfs(:);

      if( (k>1) || isempty(pars.U0) )
        relres = pcgrelres(@(x) WrDWfDT(x, InputDims, WF, WR), wfs, z);
      else
        z = 0*wfs;
        relres = 1e-1;
      end
    end

    if( pars.adaptPCGtol )
      pcgtol = relres/pars.rrs;
    else
      pcgtol = pars.pcgtol_ini;
    end

    % Compute and display functional values
    iterdisp(S, [], mean(lambda(:)), pars.p, pars.q, U, WF, WR, k, relres, [], [], ...
             pars.sbstflg, dflag);  % FIXME: mean(lambda(:))


    % -----------------------
    % Update current solution
    % -----------------------

    switch(pars.weight_scheme)

      % ---
      case{NMP_WEIGHTS_THRESHOLD}

      if pars.sbstflg, % Use substitution (indirect) method
        [v, flg, ~, pit] = pcg(@(x) f_IWDTWDW(x, InputDims, WF, WR, lambda), ...
                                 wfs, pcgtol, pars.pcgitn,[],[],v);
        u = vec(WF).*v;

      else        % Use standard (direct) method
        [u, flg, ~, pit] = pcg(@(x) f_IDTWD(x, InputDims, WF, WR, lambda), ...
                                 wfs, pcgtol, pars.pcgitn, [], [], u);
      end


      % ---
      case{NMP_WEIGHTS_MIL}


      [z, flg , ~, pit] = pcg(@(x) WrDWfDT(x, InputDims, WF, WR), ...
                                wfs, pcgtol, pars.pcgitn,[],[],z);

      tmp = WF(:).*DxTDyT(z, InputDims);
      u = s - ( tmp.*(abs(tmp)<=1) + sign(tmp).*(abs(tmp)>1) );


%          u = s - WF(:).*DxTDyT(z, InputDims);

    end % _END_ SWITCH

    U = reshape(u, InputDims);


    % Compute and display functional values and CG status
    % FIXME: mean(lambda(:)) !?
    is = iterdisp(S, [], mean(lambda(:)), pars.p, pars.q, U, WF, WR, ...
                  [], relres, flg, pit, pars.sbstflg, dflag); % FIXME: mean(lambda(:))
    itstat = [itstat; [k is]];

    if(pit == 0) 
        break;
    end

  end


return;

% ===========================================================================
% ===========================================================================


function[U, itstat] = irntv_standard_deconv(S, KC, lambda, pars, InputDims, dflag, fHandles)


  % Start clock
  tic;

  % Initialise iteration status vector
  itstat = [];

  nmpdef;

  fF         = fHandles{1};
  fF_var3    = fHandles{2};
  fR         = fHandles{3};
  fR_var4    = fHandles{4};
  f_IDTD     = fHandles{5};
  f_IDTWD    = fHandles{6};
  f_IWDTWDW  = fHandles{7};
  f_KTKDTD   = fHandles{8};
  f_KTWKDTWD = fHandles{9};

  % KC is cell array of K and K^T function handles
  K = KC{1};
  KT = KC{2};
 

  %-----------------------------------------------------------------------
  %-----------------------------------------------------------------------

  %-----------------------------------------------------------------------

  
  if isempty(pars.U0),
    % Construct initial solution
    kts = vec(KT(S));
    [u, flg, rlr, pit] = pcg(@(x) f_KTKDTD(x, InputDims, K, KT, lambda), ...
                             kts, pars.pcgtol_ini, pars.pcgitn);
    U = reshape(u, InputDims);

  else
    U = pars.U0;
    u = U(:);
  end

  % Iterate
  for k = 1:pars.loops,


    % Set weight matrices
    WR = fR(U, InputDims, pars.q, fR_var4);


    % Compute relative residual of current solution
    if(pars.variant == NMP_TV_STANDARD) 
      ku_s = (K(U)-S);

    else % If we enter here, it means that size(lambda) = size(Input)

      ku_s = (K(U)-S).*lambda;
    end

    WF = fF(ku_s, pars.p, fF_var3);
    wfs = vec(KT(WF.*S));

    if(pars.variant == NMP_TV_STANDARD) 
      relres = pcgrelres(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, lambda), wfs, u);
    else
      relres = pcgrelres(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, 1.0), wfs, u);
    end


    % tolerance for PCG
    pcgtol = relres/pars.rrs;


    % Compute and display functional values
    % FIXME: mean(lambda(:)) !?
    iterdisp(S, K, mean(lambda(:)), pars.p, pars.q, U, WF, WR, k, relres, ...
             [], [], [], dflag);

    % Update current solution
    if(pars.variant == NMP_TV_STANDARD) 
      [u, flg, rlr, pit] = pcg(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, lambda), ...
                             wfs, pcgtol, pars.pcgitn,[],[],u);
    else
      [u, flg, rlr, pit] = pcg(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, 1), ...
                             wfs, pcgtol, pars.pcgitn,[],[],u);
    end


    U = reshape(u, InputDims);

    % Compute and display functional values and CG status
    % FIXME: mean(lambda(:)) !?
    is = iterdisp(S, K, mean(lambda(:)), pars.p, pars.q, U, WF, WR, ...
                  [], relres, flg, pit, [], dflag);
    itstat = [itstat; is];

  end



return;

% ===========================================================================
% ===========================================================================


function [U, itstat] = irntv_l1l2(S, KC, lambda, pars, InputDims, fF, fF_var3, fR, fR_var4, dflag)

nmpdef;


% Start clock
tic;

% Initialise iteration status vector
itstat = [];

% 1D representation of input image
s = S(:); 


%

if(length(InputDims) == 2)

  lambda_l1 = lambda(:,:,1);
  lambda_l2 = lambda(:,:,2);

  % FIXME: code more efficiently
  mask_l1   = lambda(:,:,1) > 0;
  scalar_l2 = sum(lambda_l2(:)) / ( sum(1-mask_l1(:)) );

else

  for d = 1:InputDims(3)
      lambda_l1(:,:,d) = lambda(:,:,d);
      lambda_l2(:,:,d) = lambda(:,:,3+d);
      mask_l1(:,:,d)   = lambda(:,:,d) > 0;
  end

  % FIXME: code more efficiently
  [k1,k2] = find(lambda_l2(:,:,1)>0);
  scalar_l2 = lambda_l2(k1(1), k2(1), 1);

end

  lambda = 1;


%
switch(pars.variant)	% check variant

  case{NMP_TV_LamdaAdapt}
    if isempty(KC)
      f_IDTD = @(var1, var2, var3) IDTD_lambdaMtx_Fid(var1, var2, var3);
      f_IDTWD = @(var1, var2, var3, var4, var5) IDTWD_lambdaMtx_Fid(var1, var2, var3, var4, var5);
      f_IWDTWDW = @(var1, var2, var3, var4, var5) IWDTWDW_lambdaMtx_Fid(var1, var2, var3, var4, var5);
    else
%  
%        f_KTKDTD = @(var1, var2, var3, var4, var5) KTKDTD(var1, var2, var3, var4, var5); 
%        f_KTWKDTWD = @(var1, var2, var3, var4, var5, var6, var7) KTWKDTWD(var1, var2, var3, var4, var5, var6, var7);
%  
      disp('TV_L1L2 Deconvolution case: not implemented');
    end
  
  otherwise
    error('For l1l2TV only the NMP_TV_LamdaAdapt is implemented');


end

if isempty(KC),  % Denoising problem


  if isempty(pars.U0),
    % Construct initial solution                        
    [u, flg, rlr, pit] = pcg(@(x) f_IDTD(x, InputDims, pars.lambda_ini ), s, ...
                             pars.pcgtol_ini, pars.pcgitn, [], [], s);
    U = reshape(u, InputDims);
  else
    U = pars.U0;
    u = U(:);
  end

  pars.sbstflg = 1;

  % Iterate
  for k = 1:pars.loops,

if 1
    if( rem(k,2) == 0) 
      S0 = adaptMedian_colfilt(U);

%        lambda_l1  =  lambda_l1 + 0.85*(S0 == 3) + ...
%                      1.1*(S0 == 5) + 1.2*(S0 == 7) + 1.4*(S0 == 9);
      lambda_l1 = 1.2*lambda_l1 +  1.5*(S0 > 0);
      lambda_l2 = scalar_l2*(S0 == 0)*(1.05^(k-1));
      mask_l1   = lambda_l1 > 0;
    end
else

    lambda_l1 = lambda_l1*1.5;
end
    % Set weight matrices
    if( (k==1) && (sum(abs(U(:)-S(:))) < 0.05/(InputDims(1)*InputDims(2)) ) )
      WF_L1 = mask_l1.*( (2)*(0.1)^( (2-1)/2 )*ones(InputDims) )./( lambda_l1 + (1-mask_l1) );
      WF_L2 = (1-mask_l1)./( lambda_l2 + mask_l1 );
    else
      WF_L1 = mask_l1.*fF(U-S, 1, fF_var3)./( lambda_l1 + (1-mask_l1) );
      WF_L2 = (1-mask_l1)./( lambda_l2 + mask_l1 );
    end

    WR = fR(U, InputDims, pars.q, fR_var4);
    

    % Compute relative residual of current solution
    if pars.sbstflg, % Use substitution (indirect) method

        WF = 1./sqrt(WF_L1 + WF_L2);
                                                  % WFN2.*DxT(WR.*Dx(WFN2.*V))
        wfs = s./vec(WF);
        v = u./vec(WF);
        relres = pcgrelres(@(x) f_IWDTWDW(x, InputDims, WF, WR, lambda), wfs, v);
    else        % Use standard (direct) method

        WF = (WF_L1 + WF_L2);

        wfs = vec(WF).*s;
        relres = pcgrelres(@(x) f_IDTWD(x, InputDims, WF, WR, lambda), wfs, u);
    end
    pcgtol = (0.95^(k-1))*relres/pars.rrs;

    % Compute and display functional values
    iterdisp(S, [], mean(lambda(:)), pars.p, pars.q, U, WF, WR, k, relres, [], [], ...
             pars.sbstflg, dflag);  % FIXME: mean(lambda(:))

    % Update current solution
    if pars.sbstflg, % Use substitution (indirect) method
        [v, flg, rlr, pit] = pcg(@(x) f_IWDTWDW(x, InputDims, WF, WR, lambda), ...
                                 wfs, pcgtol, pars.pcgitn,[],[],v);
        V = reshape(v, InputDims);
        u = vec(WF).*v;
    else        % Use standard (direct) method
        [u, flg, rlr, pit] = pcg(@(x) f_IDTWD(x, InputDims, WF, WR, lambda), ...
                                 wfs, pcgtol, pars.pcgitn, [], [], u);
    end

    U = reshape(u, InputDims);

    % Compute and display functional values and CG status
    % FIXME: mean(lambda(:)) !?
    is = iterdisp(S, [], mean(lambda(:)), pars.p, pars.q, U, WF, WR, ...
                  [], relres, flg, pit, pars.sbstflg, dflag); % FIXME: mean(lambda(:))
    itstat = [itstat; [k is]];

  end

else % General inverse problem with linear operator K


  % KC is cell array of K and K^T function handles
  K = KC{1};
  KT = KC{2};
 

  %-----------------------------------------------------------------------
  %-----------------------------------------------------------------------

  %-----------------------------------------------------------------------

  
  if isempty(pars.U0),
    % Construct initial solution
    kts = vec(KT(S));
    [u, flg, rlr, pit] = pcg(@(x) f_KTKDTD(x, InputDims, K, KT, lambda), ...
                             kts, pars.pcgtol_ini, pars.pcgitn);
    U = reshape(u, InputDims);

  else
    U = pars.U0;
    u = U(:);
  end

  % Iterate
  for k = 1:pars.loops,


    % Set weight matrices

    WR = fR(U, InputDims, pars.q, fR_var4);


    % Compute relative residual of current solution

    if(pars.variant == NMP_TV_STANDARD) 

      ku_s = (K(U)-S);
    else
      ku_s = (K(U)-S).*lambda;
    end

      WF = fF(ku_s, pars.p, fF_var3);
      wfs = vec(KT(WF.*S));

%===============================================================================
%===============================================================================
if 0
    ku_s = (K(U)-S);

    if(pars.variant == NMP_TV_STANDARD) 

      WF = fF(ku_s, pars.p, fF_var3);
      wfs = vec(KT(WF.*S));

    else % Adapt

      if( (k==1) && (sum(abs(ku_s(:))) < 0.05/(InputDims(1)*InputDims(2)) ) )
        WF = fF(lambda.*1e-4, pars.p, fF_var3);
      else
        WF = fF(lambda.*ku_s, pars.p, fF_var3);
      end
      wfs = vec(KT(lambda.*WF.*S));

    end
end
%===============================================================================
%===============================================================================



    if(pars.variant == NMP_TV_STANDARD) 
      relres = pcgrelres(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, lambda), wfs, u);
    else
      relres = pcgrelres(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, 1.0), wfs, u);
    end

    pcgtol = relres/pars.rrs;

    % Compute and display functional values
    % FIXME: mean(lambda(:)) !?
    iterdisp(S, K, mean(lambda(:)), pars.p, pars.q, U, WF, WR, k, relres, ...
             [], [], [], dflag);

    % Update current solution
    if(pars.variant == NMP_TV_STANDARD) 
      [u, flg, rlr, pit] = pcg(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, lambda), ...
                             wfs, pcgtol, pars.pcgitn,[],[],u);
    else
      [u, flg, rlr, pit] = pcg(@(x) f_KTWKDTWD(x, InputDims, K, KT, WF, WR, 1), ...
                             wfs, pcgtol, pars.pcgitn,[],[],u);
    end


    U = reshape(u, InputDims);

    % Compute and display functional values and CG status
    % FIXME: mean(lambda(:)) !?
    is = iterdisp(S, K, mean(lambda(:)), pars.p, pars.q, U, WF, WR, ...
                  [], relres, flg, pit, [], dflag);
    itstat = [itstat; is];

  end

end


% Final part of status display
iterdispclose(dflag);

return


% ===========================================================================
% ===========================================================================


function [U, itstat] = irntv_nqp(S, KC, lambda, pars, InputDims, fF, fF_var3, fR, fR_var4, dflag)

% Start clock
tic;

% Initialise iteration status vector
itstat = [];

% 1D representation of input image
s = S(:); 

if isempty(KC),  % Denoising problem

  %-----------------------------------------------------------------------
  %-----------------------------------------------------------------------

  if isempty(pars.U0),
    % Initial solution for IRN-NQP 

    % Construct initial solution
%      kts = vec(S);
%      [u, flg, rlr, pit] = pcg(@(x) IDTD(x, InputDims, pars.lambda_ini), ...
%                               kts, pars.pcgtol_ini, pars.pcgitn);
%      U = reshape(u, InputDims);

    U = S;
  else
    U = pars.U0;
  end


  % Iterate
  for k = 1:pars.loops,

    % Set weight matrices
    WR = fR(U, InputDims, pars.q, fR_var4);
    

    if( (k==1) && (pars.p == 2) )
      beta  = S;
      beta2 = beta.*beta;
      WF = [];
    end

    if( pars.p ~= 2 )

      %% WF = fF(U-S, pars.p, fF_var3); % Set weight matrices
      if k>=2
        WF = fF(U-S, pars.p, fF_var3); % Set weight matrices
      else
        WF=1;
      end

      beta  = WF.*S;
      beta2 = beta.*beta;
    end


  % Auto-adapt NQP tolerance (ala Inexact Newton)
  if(pars.gamma_NQP <= 0)

    local_eps = [];
    NQPeps   = [];

    for NQPit=1:pars.loops_NQP,

      xi_pos = fApos(U, WR, WF, lambda, InputDims);
      xi_neg = fAneg(U, WR, lambda, InputDims);

      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      % U = (U.*num)./den;
      U = min( (U.*num)./den, pars.vmax_NQP);

    end

  else % Auto-adapt NQP tolerance (ala Inexact Newton)

    NQPit=0;

    xi_pos = fApos(U, WR, WF, lambda, InputDims);
    xi_neg = fAneg(U, WR, lambda, InputDims);

    local_eps = pars.gamma_NQP*( norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:)) )^pars.alpha_NQP;
    NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));

    while( NQPeps > local_eps )

      NQPit = NQPit + 1;

      if(NQPit > pars.loops_NQP ) 
         NQPit=NQPit-1; 
         break; 
      end

      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      U = min( (U.*num)./den, pars.vmax_NQP);

      xi_pos = fApos(U, WR, WF, lambda, InputDims);
      xi_neg = fAneg(U, WR, lambda, InputDims);

      NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));

    end

  end % Auto-adapt NQP tolerance (ala Inexact Newton)

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
if 0
    % Compute relative residual of current solution
    if pars.sbstflg, % Use substitution (indirect) method
        wfs = s./vec(WF);
        v = u./vec(WF);
        relres = pcgrelres(@(x) IWDTWDW(x, InputDims, WF, WR, lambda), wfs, v);
    else        % Use standard (direct) method
        wfs = vec(WF).*s;
        relres = pcgrelres(@(x) IDTWD(x, InputDims, WF, WR, lambda), wfs, u);
    end
    pcgtol = relres/pars.rrs;

    % Compute and display functional values
    iterdisp(S, [], lambda, pars.p, pars.q, U, WF, WR, k, relres, [], [], ...
             pars.sbstflg, dflag);

    % Update current solution
    if pars.sbstflg, % Use substitution (indirect) method
        [v, flg, rlr, pit] = pcg(@(x) IWDTWDW(x, InputDims, WF, WR, lambda), ...
                                 wfs, pcgtol, pars.pcgitn,[],[],v);
        V = reshape(v, InputDims);
        u = vec(WF).*v;
    else        % Use standard (direct) method
        [u, flg, rlr, pit] = pcg(@(x) IDTWD(x, InputDims, WF, WR, lambda), ...
                                 wfs, pcgtol, pars.pcgitn, [], [], u);
    end

    U = reshape(u, InputDims);
end
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------

    % Compute and display functional values and CG status
    is = iterdisp_nqp(S, [], lambda, pars.p, pars.q, U, WF, WR, k, NQPit, local_eps, NQPeps, dflag);
    itstat = [itstat; [k is]];



  end

else % General inverse problem with linear operator K

  % KC is cell array of K and K^T function handles
  K = KC{1};
  KT = KC{2};
   
  if isempty(pars.U0),

    % Construct initial solution
    kts = vec(KT(S));
    [u, flg, rlr, pit] = pcg(@(x) KTKDTD(x, InputDims, K, KT, pars.lambda_ini), ...
                             kts, pars.pcgtol_ini, pars.pcgitn);
    U = reshape(u, InputDims);

  else

    U = pars.U0;
    if( (min(U(:)) == 0) ) 
        U = U+1; 
    end %

    disp('The NQP solver need values diferent from zero to consider');
    disp('them as part of the solution');
    disp('Your initial solution has zeros... we have added "+1"');
    disp(' ');
    disp('Consider calling inrtv with no initial solution and/or')
    disp('irntv(S+1, KC, lambda, pars)-1 ... note the "+1" and "-1"');
    disp(' ');
  end


  % Iterate
  for k = 1:pars.loops,


    % Set weight matrices
    WR = fR(U, InputDims, pars.q, fR_var4);

    if( pars.p == 2) 
      beta  = KT(S);
      beta2 = beta.*beta;
      WF = [];
    else
      WF = fF(K(U)-S, pars.p, fF_var3); % Set weight matrices
      beta  = KT(WF.*S);
      beta2 = beta.*beta;
    end


  % Auto-adapt NQP tolerance (ala Inexact Newton)
  if(pars.gamma_NQP <= 0)

    local_eps = [];
    NQP_eps   = [];

    for NQPit=1:pars.loops_NQP,

      xi_pos = fApos_K1(U, WR, WF, K, KT, lambda, InputDims);
      xi_neg = fAneg(U, WR, lambda, InputDims);


      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      % U = (U.*num)./den;
      U = min( (U.*num)./den, pars.vmax_NQP);

    end

  else % Auto-adapt NQP tolerance (ala Inexact Newton)

    NQPit=0;
    xi_pos = fApos_K1(U, WR, WF, K, KT, lambda, InputDims);
    xi_neg = fAneg(U, WR, lambda, InputDims);

    local_eps = pars.gamma_NQP*( norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:)) )^pars.alpha_NQP;

    NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));
    while( NQPeps > local_eps)

      NQPit = NQPit + 1;

      if(NQPit > pars.loops_NQP ) 
         NQPit=NQPit-1; 
         break; 
      end

      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      U = min( (U.*num)./den, pars.vmax_NQP);

      xi_pos = fApos_K1(U, WR, WF, K, KT, lambda, InputDims);
      xi_neg = fAneg(U, WR, lambda, InputDims);

      NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));
    end

  end % Auto-adapt NQP tolerance (ala Inexact Newton)

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
if 0
    % Compute relative residual of current solution
    wfs = vec(KT(WF.*S));
    relres = pcgrelres(@(x) KTWKDTWD(x, InputDims, K, KT, WF, WR, lambda), ...
                       wfs, u);
    pcgtol = relres/pars.rrs;

    % Compute and display functional values
    iterdisp(S, K, lambda, pars.p, pars.q, U, WF, WR, k, relres, ...
             [], [], [], dflag);

    % Update current solution
    [u, flg, rlr, pit] = pcg(@(x) KTWKDTWD(x, InputDims, K, KT, WF, WR, lambda), ...
                             wfs, pcgtol, pars.pcgitn,[],[],u);
    U = reshape(u, InputDims);

end
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------


    % Compute and display functional values and NQP solver status
    is = iterdisp_nqp(S, [], lambda, pars.p, pars.q, U, WF, WR, k, NQPit, local_eps, NQPeps, dflag);


  end

end


% Final part of status display
iterdispclose(dflag);

return



% ===========================================================================
% ===========================================================================

function [U, itstat] = irntv_poisson(S, KC, lambda, pars, InputDims, fR, fR_var4, dflag)

%  ----------

nmpdef;

%  ----------

% Start clock
tic;

% Initialise iteration status vector
itstat = [];

% 1D representation of input image
s = S(:); 



if isempty(KC),  % Denoising problem

  %-----------------------------------------------------------------------
  %-----------------------------------------------------------------------

  if isempty(pars.U0),
    % Initial solution for IRN-NQP 
    U = S;
  else
    U = pars.U0;
  end

    

  % Iterate
  for k = 1:pars.loops,

    % Set weight matrices
    WR = fR(U, InputDims, pars.q, fR_var4);

    % Construct matrix H
    switch(pars.variant)

      case{NMP_TV_LamdaAdapt}
        H = lambda.*(U>pars.epsF_cutoff).*S./( (U>pars.epsF_cutoff).*(U.*U) + (U<=pars.epsF_cutoff) );
        beta = -lambda.*( 1 - 2*(U>pars.epsF_cutoff).*S./( (U>pars.epsF_cutoff).*(U) + (U<=pars.epsF_cutoff) ) );
        lambdaPoisson = 1.0;

      otherwise
        H = (U>pars.epsF_cutoff).*S./( (U>pars.epsF_cutoff).*(U.*U) + (U<=pars.epsF_cutoff) );
        beta = -( 1 - 2*(U>pars.epsF_cutoff).*S./( (U>pars.epsF_cutoff).*(U) + (U<=pars.epsF_cutoff) ) );
        lambdaPoisson = lambda;

    end

    WF = H;

    beta2 = beta.*beta;


  % Auto-adapt NQP tolerance (ala Inexact Newton)
  if(pars.gamma_NQP <= 0)

    local_eps = [];
    NQP_eps   = [];

    for NQPit=1:pars.loops_NQP,

      xi_pos = fApos(U, WR, WF, lambdaPoisson, InputDims);
      xi_neg = fAneg(U, WR, lambdaPoisson, InputDims);

      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      % U = (U.*num)./den;
      U = min( (U.*num)./den, pars.vmax_NQP);

    end

  else % Auto-adapt NQP tolerance (ala Inexact Newton)

    NQPit=0;


    xi_pos = fApos(U, WR, WF, lambdaPoisson, InputDims);
    xi_neg = fAneg(U, WR, lambdaPoisson, InputDims);

    local_eps = pars.gamma_NQP*( norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:)) )^pars.alpha_NQP;
    NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));

    while( NQPeps > local_eps)

      NQPit = NQPit + 1;

      if(NQPit > pars.loops_NQP ) 
         NQPit=NQPit-1; 
         break; 
      end

      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      U = min( (U.*num)./den, pars.vmax_NQP);

      xi_pos = fApos(U, WR, WF, lambdaPoisson, InputDims);
      xi_neg = fAneg(U, WR, lambdaPoisson, InputDims);

      NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));

    end

  end % Auto-adapt NQP tolerance (ala Inexact Newton)

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------

    % Compute and display functional values and CG status
    is = iterdisp_poisson(S, [], lambdaPoisson, pars.p, pars.q, U, WF, WR, k, NQPit, local_eps, NQPeps, dflag);
    itstat = [itstat; [k is]];



  end

else % General inverse problem with linear operator K

  % KC is cell array of K and K^T function handles
  K = KC{1};
  KT = KC{2};
   
  if isempty(pars.U0),
    U = S;
  else
    U = pars.U0;
  end


  % Iterate
  for k = 1:pars.loops,


    % Set weight matrices
    WR = fR(U, InputDims, pars.q, fR_var4);


      KU = K(U);
      % H is implicit
      
    switch(pars.variant)

      case{NMP_TV_LamdaAdapt}
        WF = lambda.*(KU>pars.epsF_cutoff).*S./( (KU>pars.epsF_cutoff).*(KU.*KU) + (KU<=pars.epsF_cutoff) );
        beta = -KT( lambda.*( 1 - 2*(KU>pars.epsF_cutoff).*S./( (KU>pars.epsF_cutoff).*(KU) + (KU<=pars.epsF_cutoff) ) ) );
        lambdaPoisson = 1.0;

      otherwise
        WF = (KU>pars.epsF_cutoff).*S./( (KU>pars.epsF_cutoff).*(KU.*KU) + (KU<=pars.epsF_cutoff) );
        beta = -KT( ( 1 - 2*(KU>pars.epsF_cutoff).*S./( (KU>pars.epsF_cutoff).*(KU) + (KU<=pars.epsF_cutoff) ) ) );
        lambdaPoisson = lambda;

    end

      beta2 = beta.*beta;


  % Auto-adapt NQP tolerance (ala Inexact Newton)
  if(pars.gamma_NQP <= 0)

    local_eps = [];
    NQP_eps   = [];

    for NQPit=1:pars.loops_NQP,

      xi_pos = fApos_K1(U, WR, WF, K, KT, lambdaPoisson, InputDims);
      xi_neg = fAneg(U, WR, lambdaPoisson, InputDims);


      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      % U = (U.*num)./den;
      U = min( (U.*num)./den, pars.vmax_NQP);

    end

  else % Auto-adapt NQP tolerance (ala Inexact Newton)


    NQPit=0;
    xi_pos = fApos_K1(U, WR, WF, K, KT, lambdaPoisson, InputDims);
    xi_neg = fAneg(U, WR, lambdaPoisson, InputDims);

    local_eps = pars.gamma_NQP*( norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:)) )^pars.alpha_NQP;

    NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));
    while( NQPeps > local_eps)

      NQPit = NQPit + 1;

      if(NQPit > pars.loops_NQP ) 
         NQPit=NQPit-1; 
         break; 
      end

      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      U = min( (U.*num)./den, pars.vmax_NQP);

      xi_pos = fApos_K1(U, WR, WF, K, KT, lambdaPoisson, InputDims);
      xi_neg = fAneg(U, WR, lambdaPoisson, InputDims);

      NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));
    end

  end % Auto-adapt NQP tolerance (ala Inexact Newton)

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------


    % Compute and display functional values and NQP solver status
    is = iterdisp_nqp(S, [], lambdaPoisson, pars.p, pars.q, U, WF, WR, k, NQPit, local_eps, NQPeps, dflag);


  end

end


% Final part of status display
iterdispclose(dflag);

return



% ===========================================================================
% ===========================================================================



function [U, itstat] = irntv_speckle(S, KC, lambda, pars, InputDims, fR, fR_var4, dflag)

nmpdef;
% Start clock
tic;

% Initialise iteration status vector
itstat = [];

% 1D representation of input image
s = S(:); 

if isempty(KC),  % Denoising problem

  %-----------------------------------------------------------------------
  %-----------------------------------------------------------------------

  if isempty(pars.U0),
    % Initial solution for IRN-NQP 
    switch(pars.speckle_variant)

      case{NMP_TV_SPECKLE_AA}
      U = S;

      case{NMP_TV_SPECKLE_SO}
      U = log(S);

    end % _END_ SWITCH

  else

    U = pars.U0;

  end % _END_ IF(isempty(pars.U0))

    

  % Iterate
  for k = 1:pars.loops,


    % Set weight matrices
    WR = fR(U, InputDims, pars.q, fR_var4);


      % Model NMP_TV_SPECKLE_AA}
      UInv = (U>pars.epsF_cutoff)./( (U>pars.epsF_cutoff).*(U) + (U<=pars.epsF_cutoff) );
      UInv2 = (U>pars.epsF_cutoff)./( (U>pars.epsF_cutoff).*(U.*U) + (U<=pars.epsF_cutoff) );

      H = UInv2.*( 2*S.*UInv - 1);

      H = H.*(H>pars.epsF_cutoff) + pars.epsF_cutoff*(H<=0);
      WF = H;

      beta = -( 2*UInv - 3*S.*UInv2 );

      beta2 = beta.*beta;


    % Auto-adapt NQP tolerance (ala Inexact Newton)
    if(pars.gamma_NQP <= 0)

      local_eps = [];
      NQP_eps   = [];

      for NQPit=1:pars.loops_NQP,

        xi_pos = fApos(U, WR, WF, lambda, InputDims);
        xi_neg = fAneg(U, WR, lambda, InputDims);

        num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

        den = 2*xi_pos + (xi_pos <= 0);

        % U = (U.*num)./den;
        U = min( (U.*num)./den, pars.vmax_NQP);

      end % _END_ FOR(NQPit)

    else % Auto-adapt NQP tolerance (ala Inexact Newton)

      NQPit=0;


      xi_pos = fApos(U, WR, WF, lambda, InputDims);
      xi_neg = fAneg(U, WR, lambda, InputDims);

      local_eps = pars.gamma_NQP*( norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:)) )^pars.alpha_NQP;
      NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));

      while( NQPeps > local_eps)

        NQPit = NQPit + 1;

        if(NQPit > pars.loops_NQP ) 
          NQPit=NQPit-1; 
          break; 
        end

        num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

        den = 2*xi_pos + (xi_pos <= 0);

        U = min( (U.*num)./den, pars.vmax_NQP);

        xi_pos = fApos(U, WR, WF, lambda, InputDims);
        xi_neg = fAneg(U, WR, lambda, InputDims);

        NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));

      end % _END_ WHILE


    end % Auto-adapt NQP tolerance (ala Inexact Newton)


% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------

    % Compute and display functional values and CG status
    is = iterdisp_poisson(S, [], lambda, pars.p, pars.q, U, WF, WR, k, NQPit, local_eps, NQPeps, dflag);
    itstat = [itstat; [k is]];

  end

%    if(pars.speckle_variant == NMP_TV_SPECKLE_SO) U = exp(U); end


else % General inverse problem with linear operator K

  % KC is cell array of K and K^T function handles
  K = KC{1};
  KT = KC{2};
   
  if isempty(pars.U0),
    U = S;
  else
    U = pars.U0;
  end


  % Iterate
  for k = 1:pars.loops,


    % Set weight matrices
    WR = fR(U, InputDims, pars.q, fR_var4);

    % ------------------------

      KU = K(U);

      KUInv2 = (KU>pars.epsF_cutoff)./( (KU>pars.epsF_cutoff).*(KU.*KU) + (KU<=pars.epsF_cutoff) );
      KUInv3 = (KU>pars.epsF_cutoff)./( (KU>pars.epsF_cutoff).*(KU.*KU.*KU) + (KU<=pars.epsF_cutoff) );

      H = (2*S - KU).*KUInv3;
      H = H.*(H>pars.epsF_cutoff) + pars.epsF_cutoff*(H<=0);
      WF = H;

      beta = -KT( (2*KU - 3*S).*KUInv2 );

    % ------------------------

    beta2 = beta.*beta;


  % Auto-adapt NQP tolerance (ala Inexact Newton)
  if(pars.gamma_NQP <= 0)

    local_eps = [];
    NQP_eps   = [];

    for NQPit=1:pars.loops_NQP,

      xi_pos = fApos_K1(U, WR, WF, K, KT, lambda, InputDims);
      xi_neg = fAneg(U, WR, lambda, InputDims);


      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      % U = (U.*num)./den;
      U = min( (U.*num)./den, pars.vmax_NQP);

    end

  else % Auto-adapt NQP tolerance (ala Inexact Newton)


    NQPit=0;
    xi_pos = fApos_K1(U, WR, WF, K, KT, lambda, InputDims);
    xi_neg = fAneg(U, WR, lambda, InputDims);

    local_eps = pars.gamma_NQP*( norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:)) )^pars.alpha_NQP;

    NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));
    while( NQPeps > local_eps)

      NQPit = NQPit + 1;

      if(NQPit > pars.loops_NQP ) 
         NQPit=NQPit-1; 
         break; 
      end

      num = beta + sqrt( beta2 + 4*xi_pos.*xi_neg );

      den = 2*xi_pos + (xi_pos <= 0);

      U = min( (U.*num)./den, pars.vmax_NQP);

      xi_pos = fApos_K1(U, WR, WF, K, KT, lambda, InputDims);
      xi_neg = fAneg(U, WR, lambda, InputDims);

      NQPeps = norm(xi_pos(:)-xi_neg(:) - beta(:))/norm(beta(:));
    end

  end % Auto-adapt NQP tolerance (ala Inexact Newton)

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------


    % Compute and display functional values and NQP solver status
    is = iterdisp_nqp(S, [], lambda, pars.p, pars.q, U, WF, WR, k, NQPit, local_eps, NQPeps, dflag);


  end

end


% Final part of status display
iterdispclose(dflag);

return




% ===========================================================================
% ===========================================================================



function[x] = fApos(u, WR, WF, lambda, UDims)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)

    x(:,:,1) = fDyTpos(WR.*fDypos(u(:,:,1))) + fDxTpos(WR.*fDxpos(u(:,:,1))) + ...
               fDyTneg(WR.*fDyneg(u(:,:,1))) + fDxTneg(WR.*fDxneg(u(:,:,1)));

    for k=2:UDims(3)
      x(:,:,k) = fDyTpos(WR.*fDypos(u(:,:,k))) + fDxTpos(WR.*fDxpos(u(:,:,k))) + ...
                 fDyTneg(WR.*fDyneg(u(:,:,k))) + fDxTneg(WR.*fDxneg(u(:,:,k)));
    end

  else

    x = fDyTpos(WR.*fDypos(u)) + fDxTpos(WR.*fDxpos(u)) + ...
        fDyTneg(WR.*fDyneg(u)) + fDxTneg(WR.*fDxneg(u));

  end

  if( isempty(WF) )
      x = u + lambda.*x;
  else
      x = WF.*u + lambda.*x;
  end


return

%========================================================
%========================================================

function[x] = fApos_LambdaMtx(u, WR, WF, lambda, UDims)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)

    x(:,:,1) = fDyTpos(WR.*fDypos(u(:,:,1))) + fDxTpos(WR.*fDxpos(u(:,:,1))) + ...
               fDyTneg(WR.*fDyneg(u(:,:,1))) + fDxTneg(WR.*fDxneg(u(:,:,1)));

    for k=2:UDims(3)
      x(:,:,k) = fDyTpos(WR.*fDypos(u(:,:,k))) + fDxTpos(WR.*fDxpos(u(:,:,k))) + ...
                 fDyTneg(WR.*fDyneg(u(:,:,k))) + fDxTneg(WR.*fDxneg(u(:,:,k)));
    end

  else

    x = fDyTpos(WR.*fDypos(u)) + fDxTpos(WR.*fDxpos(u)) + ...
        fDyTneg(WR.*fDyneg(u)) + fDxTneg(WR.*fDxneg(u));

  end

  if( isempty(WF) )
      x = u.*lambda + x;
  else
      x = WF.*u.*lambda + x;
  end

return


%========================================================
%========================================================


function[x] = fTVpos(u, WR, lambda, UDims)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)

    x(:,:,1) = fDyTpos(WR.*fDypos(u(:,:,1))) + fDxTpos(WR.*fDxpos(u(:,:,1))) + ...
               fDyTneg(WR.*fDyneg(u(:,:,1))) + fDxTneg(WR.*fDxneg(u(:,:,1)));

    for k=2:UDims(3)
      x(:,:,k) = fDyTpos(WR.*fDypos(u(:,:,k))) + fDxTpos(WR.*fDxpos(u(:,:,k))) + ...
                 fDyTneg(WR.*fDyneg(u(:,:,k))) + fDxTneg(WR.*fDxneg(u(:,:,k)));
    end

  else

    x = fDyTpos(WR.*fDypos(u)) + fDxTpos(WR.*fDxpos(u)) + ...
        fDyTneg(WR.*fDyneg(u)) + fDxTneg(WR.*fDxneg(u));

  end

  x = lambda.*x;

return

%========================================================
%========================================================

function[x] = fAneg(u, WR, lambda, UDims)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)

    x(:,:,1) = fDyTneg(WR.*fDypos(u(:,:,1))) + fDxTneg(WR.*fDxpos(u(:,:,1))) + ...
               fDyTpos(WR.*fDyneg(u(:,:,1))) + fDxTpos(WR.*fDxneg(u(:,:,1)));

    for k=2:UDims(3)

      x(:,:,k) = fDyTneg(WR.*fDypos(u(:,:,k))) + fDxTneg(WR.*fDxpos(u(:,:,k))) + ...
                 fDyTpos(WR.*fDyneg(u(:,:,k))) + fDxTpos(WR.*fDxneg(u(:,:,k)));
    end

  else

    x = fDyTneg(WR.*fDypos(u)) + fDxTneg(WR.*fDxpos(u)) + ...
        fDyTpos(WR.*fDyneg(u)) + fDxTpos(WR.*fDxneg(u));

  end

  x = lambda.*x;


return

%========================================================
%========================================================

function[x] = fApos_K1(u, WR, WF, K, KT, lambda, UDims)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)

    x(:,:,1) = fDyTpos(WR.*fDypos(u(:,:,1))) + fDxTpos(WR.*fDxpos(u(:,:,1))) + ...
               fDyTneg(WR.*fDyneg(u(:,:,1))) + fDxTneg(WR.*fDxneg(u(:,:,1)));

    for k=2:UDims(3)
      x(:,:,k) = fDyTpos(WR.*fDypos(u(:,:,k))) + fDxTpos(WR.*fDxpos(u(:,:,k))) + ...
                 fDyTneg(WR.*fDyneg(u(:,:,k))) + fDxTneg(WR.*fDxneg(u(:,:,k)));
    end

  else

    x = fDyTpos(WR.*fDypos(u)) + fDxTpos(WR.*fDxpos(u)) + ...
        fDyTneg(WR.*fDyneg(u)) + fDxTneg(WR.*fDxneg(u));

  end


  if( isempty(WF) )
    y = KT(K(u));
  else
    y = KT(WF.*K(u));
  end

  x = y + lambda*x;

return

%========================================================
%========================================================


function[x] = fApos_Kpoisson(u, WR, WF, K, KT, lambda, UDims)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)

    x(:,:,1) = fDyTpos(WR.*fDypos(u(:,:,1))) + fDxTpos(WR.*fDxpos(u(:,:,1))) + ...
               fDyTneg(WR.*fDyneg(u(:,:,1))) + fDxTneg(WR.*fDxneg(u(:,:,1)));

    for k=2:UDims(3)
      x(:,:,k) = fDyTpos(WR.*fDypos(u(:,:,k))) + fDxTpos(WR.*fDxpos(u(:,:,k))) + ...
                 fDyTneg(WR.*fDyneg(u(:,:,k))) + fDxTneg(WR.*fDxneg(u(:,:,k)));
    end

  else

    x = fDyTpos(WR.*fDypos(u)) + fDxTpos(WR.*fDxpos(u)) + ...
        fDyTneg(WR.*fDyneg(u)) + fDxTneg(WR.*fDxneg(u));

  end


  if( isempty(WF) )
    y = KT(K(u));
  else
    y = KT(WF.*K(u));
  end

  x = y + lambda*x;

return
%========================================================
%========================================================

function[x] = fDypos(u)

  [Nrows, Ncols] = size(u);

  x = zeros(size(u));

  x(1:Nrows-1,:) = u(2:Nrows,:);
  x(Nrows,:)     = 0.5*u(Nrows-2,:) + 1.5*u(Nrows,:);

return

%========================================================
%========================================================


function[x] = fDyTpos(u)

  [Nrows, Ncols] = size(u);

  x = zeros(size(u));

  x(2:Nrows,:) = u(1:Nrows-1,:);
  x(Nrows-2,:) = x(Nrows-2,:) + 0.5*u(Nrows,:);
  x(Nrows,:) = x(Nrows,:) + 1.5*u(Nrows,:);

return

%========================================================
%========================================================

function[x] = fDxpos(u)

  [Nrows Ncols] = size(u);

  x = zeros(size(u));

  x(:,1:Ncols-1) = u(:,2:Ncols);
  x(:,Ncols)     = 0.5*u(:,Ncols-2) + 1.5*u(:,Ncols);

return

%========================================================
%========================================================

function[x] = fDxTpos(u)

  [Nrows Ncols] = size(u);

  x = zeros(size(u));

  x(:,2:Ncols) = u(:,1:Ncols-1);
  x(:,Ncols-2) = x(:,Ncols-2) + 0.5*u(:,Ncols);
  x(:,Ncols) = x(:,Ncols) + 1.5*u(:,Ncols);

return

%========================================================
%========================================================

function[x] = fDyneg(u)

  [Nrows Ncols] = size(u);

  x = u;
  x(Nrows,:)     = 2*u(Nrows-1,:);

return

%========================================================
%========================================================

function[x] = fDyTneg(u)

  [Nrows Ncols] = size(u);

  x = u;
  x(Nrows-1,:)  = x(Nrows-1,:) + 2*u(Nrows,:);
  x(Nrows,:)    = 0;

return

%========================================================
%========================================================

function[x] = fDxneg(u)

  [Nrows Ncols] = size(u);

  x = u;        
  x(:,Ncols)     = 2*u(:,Ncols-1);

return

%========================================================
%========================================================

function[x] = fDxTneg(u)

  [Nrows Ncols] = size(u);

  x = u;
  x(:,Ncols-1,:)  = x(:,Ncols-1) + 2*u(:,Ncols);
  x(:,Ncols,:)    = 0;

return

%========================================================
%========================================================


% Vectorise image v
function[u] = vec(v)

  u = v(:);
  
return

%========================================================
%========================================================

% Consider operator Dy to be the y derivate of a vectorised
% image. Apply this operator to the unvectorised image and return
% a gradient image.
function[u] = Dy(v)

  u = [diff(v); zeros(1,size(v,2))];
    
return

%========================================================
%========================================================

% Consider operator Dy to be the y derivate of a vectorised
% image. Apply the transpose of this operator to the unvectorised
% image and return a gradient image.
function[u] = DyT(v)

  u0 = -v(1,:);
  u1 = -diff(v);
  u2 = v(end-1,:);
  u = [u0; u1(1:(end-1),:); u2];
    
return

%========================================================
%========================================================

% Consider operator Dx to be the x derivate of a vectorised
% image. Apply this operator to the unvectorised image and return
% a gradient image.
function[u] = Dx(v)

  u = [diff(v,1,2) zeros(size(v,1),1)];
    
return

%========================================================
%========================================================

% Consider operator Dx to be the x derivate of a vectorised
% image. Apply the transpose of this operator to the unvectorised
% image and return a gradient image.
function[u] = DxT(v)

  u = DyT(v')';

return

%========================================================
%========================================================


function[u] = DxDy(v)

  u = [Dx(v), Dy(v)];

return

%========================================================
%========================================================
% This function assumes that input is a vectorized matrix; 
% moreover, the input matrix is formed by two matrices: 
% [v1, v2] where size(v1) = size(v2)
function[u] = DxTDyT(v, vDims)

  Ncols = vDims(2);
  N     = 2*Ncols;

  if( length(vDims) == 3 )
    V = reshape(v, [vDims(1), N, vDims(3)]);
  else
    V = reshape(v, [vDims(1), N]);
  end


  u = DxT(V(:,1:Ncols)) + DyT(V(:,Ncols+1:N));

  u = u(:);

return

%========================================================
%========================================================


% Compute scalar function f_F (Fidelity weights)
function[y] = fF_MatInvLemma(x, p, epsf)

  if(p ~= 2)
    absx = abs(x);
    y =  (p/2).*( ( absx ).^(2-p) );    % CHECK 2/p!!
  else
    y = 1;
  end

return

%========================================================
%========================================================


% Compute scalar function f_F (Fidelity weights)
function[y] = fF_fixed(x, p, epsf)

  absx = abs(x);
  mask = (absx >= epsf);
  y =  (2/p).*( ( (mask.*absx) + (1-mask)*epsf ).^(p-2) );


return

%========================================================
%========================================================

function y = fF_fixed_substitute(x, p, epsf)

  absx = abs(x);
  mask = (absx >= epsf);
  y =  sqrt(p/2).*( ( (mask.*absx) + (1-mask)*epsf ).^((2-p)/2) );


return

%========================================================
%========================================================

function y = fF_adapt(x, p, percentile)

  absx = abs(x);

  max_absx = max( absx(:) );
  min_absx = min( absx(:) );

  %find histogram between max and min.
  [h bin] = hist(absx(:), min_absx: (max_absx-min_absx)/999 :max_absx);
  hacc = cumsum(h);

  [dummy pos] = max( hacc/hacc(end) > percentile ); % 0 0 .. 0 1 1 .. 1
                                                    %         ^
                                                    %         | percentile
  epsf = bin(pos+1);

  mask = (absx >= epsf);
  y =  (2/p).*( ( (mask.*absx) + (1-mask)*epsf ).^(p-2) );

return

%========================================================
%========================================================

function y = fF_adapt_substitute(x, p, percentile)

  absx = abs(x);

  max_absx = max( absx(:) );
  min_absx = min( absx(:) );

  %find histogram between max and min.
  [h bin] = hist(absx(:), min_absx: (max_absx-min_absx)/999 :max_absx);
  hacc = cumsum(h);

  [dummy pos] = max( hacc/hacc(end) > percentile ); % 0 0 .. 0 1 1 .. 1
                                                    %         ^
                                                    %         | percentile
  epsf = bin(pos+1);

  mask = (absx >= epsf);
  y = sqrt(p/2).*( ( (mask.*absx) + (1-mask)*epsf ).^((2-p)/2) );

return

%========================================================
%========================================================

% Compute scalar function f_R (Regularization weights, add epsilon)
function y = fR_addEPS(U, UDims, q, epsr)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    absx = Dx(U(:,:,1)).^2 + Dy(U(:,:,1)).^2 + epsr;
    for k=2:UDims(3)
      absx = absx + Dx(U(:,:,k)).^2 + Dy(U(:,:,k)).^2;
    end
  else
    absx = Dx(U).^2 + Dy(U).^2 + epsr;
  end


  y =  (2/q).*( ( absx  ).^((q-2)/2) );


return

%========================================================
%========================================================

% Compute scalar function f_R (Regularization weights)
function y = fR_MatInvLemma(U, UDims, q, lambda)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    absx = Dx(U(:,:,1)).^2 + Dy(U(:,:,1)).^2;
    for k=2:UDims(3)
      absx = absx + Dx(U(:,:,k)).^2 + Dy(U(:,:,k)).^2;
    end
  else
    absx = Dx(U).^2 + Dy(U).^2 ;
  end


  y =  (q/2).*( ( absx ).^((2-q)/2) )./lambda;


return

%========================================================
%========================================================

% Compute scalar function f_R (Regularization weights)
function y = fR_fixed(U, UDims, q, epsr)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    absx = Dx(U(:,:,1)).^2 + Dy(U(:,:,1)).^2;
    for k=2:UDims(3)
      absx = absx + Dx(U(:,:,k)).^2 + Dy(U(:,:,k)).^2;
    end
  else
    absx = Dx(U).^2 + Dy(U).^2 ;
  end

  mask = (absx >= epsr);

  y =  (2/q).*( ( (mask.*absx) + (1-mask)*epsr ).^((q-2)/2) );


return

%========================================================
%========================================================

% Compute scalar function f_R
function y = fR_adapt(U, UDims, q, percentile)

  if( length(UDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    absx = Dx(U(:,:,1)).^2 + Dy(U(:,:,1)).^2;
    for k=2:UDims(3)
      absx = absx + Dx(U(:,:,k)).^2 + Dy(U(:,:,k)).^2;
    end
  else
    absx = Dx(U).^2 + Dy(U).^2 ;
  end

  absx = abs( absx );

  max_absx = max( absx(:) );
  min_absx = min( absx(:) );


  [h bin] = hist(absx(:), min_absx: (max_absx-min_absx)/999 :max_absx); 
  hacc = cumsum(h);

  [dummy pos] = max( hacc/hacc(end) > percentile ); % 0 0 .. 0 1 1 .. 1
                                                    %         ^
                                                    %         | percentile
  
  if( bin(pos) > 0 ) epsr = bin(pos);
  else epsr = bin(pos+1);
  end
  
  mask = (absx >= epsr);
  y =  (2/q).*( ( (mask.*absx) + (1-mask)*epsr ).^((q-2)/2) );


return


%========================================================
%========================================================

function u = IDDT(v, vDims, lambda)

  Ncols = vDims(2);
  N     = 2*Ncols;

  if( length(vDims) == 3 )
    V = reshape(v, [vDims(1), N, vDims(3)]);
  else
    V = reshape(v, [vDims(1), N]);
  end

  U = zeros(size(V));


  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    disp('fixme... code this!');
  else

    U(:,1:Ncols) = V(:,1:Ncols)./lambda + Dx( DxT(V(:,1:Ncols)) ) + Dx( DyT(V(:,Ncols+1:N)) );
    U(:,Ncols+1:N) = V(:,Ncols+1:N)./lambda + Dy( DyT(V(:,Ncols+1:N)) ) + Dy( DxT(V(:,1:Ncols)) );

    u = U(:);

  end

return

  

%========================================================
%========================================================

function u = WDDT(v, vDims, W)

  Ncols = vDims(2);
  N     = 2*Ncols;

  if( length(vDims) == 3 )
    V = reshape(v, [vDims(1), N, vDims(3)]);
  else
    V = reshape(v, [vDims(1), N]);
  end

  U = zeros(size(V));


  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    disp('fixme... code this!');
  else

    U(:,1:Ncols) = W.*V(:,1:Ncols) + Dx( DxT(V(:,1:Ncols)) ) + Dx( DyT(V(:,Ncols+1:N)) );
    U(:,Ncols+1:N) = W.*V(:,Ncols+1:N) + Dy( DyT(V(:,Ncols+1:N)) ) + Dy( DxT(V(:,1:Ncols)) );

    u = U(:);

  end

return

  

%========================================================
%========================================================

function u = WrDWfDT(v, vDims, WF, WR)

  Ncols = vDims(2);
  N     = 2*Ncols;

  if( length(vDims) == 3 )
    V = reshape(v, [vDims(1), N, vDims(3)]);
  else
    V = reshape(v, [vDims(1), N]);
  end

  U = zeros(size(V));


  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    disp('fixme... code this!');
  else

    U(:,1:Ncols)   = WR.*V(:,1:Ncols)   + Dx( WF.*DxT(V(:,1:Ncols)) ) + Dx( WF.*DyT(V(:,Ncols+1:N)) );
    U(:,Ncols+1:N) = WR.*V(:,Ncols+1:N) + Dy( WF.*DxT(V(:,1:Ncols)) ) + Dy( WF.*DyT(V(:,Ncols+1:N)) );

    u = U(:);

  end

return

  

%========================================================
%========================================================




% Compute I + lambda*D_x^T*D_x + lambda*D_y^T*D_y for vectorised
% image (since called by pcg function)

function u = IDTD(v, vDims, lambda)
  
  V = reshape(v, vDims);

  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = V(:,:,k) + lambda*DxT(Dx(V(:,:,k))) + lambda*DyT(Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else

    U = V + lambda*DxT(Dx(V)) + lambda*DyT(Dy(V));
    u = U(:);

  end

return

%========================================================
%========================================================

function u = IDTD_lambdaMtx(v, vDims, lambda)
  
  V = reshape(v, vDims);

  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = V(:,:,k) + lambda(:,:,k).*DxT(Dx(V(:,:,k))) + lambda(:,:,k).*DyT(Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else

    U = V + lambda.*DxT(Dx(V)) + lambda.*DyT(Dy(V));
    u = U(:);

  end

return


%========================================================
%========================================================


function u = IDTD_lambdaMtx_Fid(v, vDims, lambda)
  
  V = reshape(v, vDims);

  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
%        U = lambda(:,:,k).*V(:,:,k) + DxT(Dx(V(:,:,k))) + DyT(Dy(V(:,:,k)));
      U = V(:,:,k) + DxT(Dx(V(:,:,k))) + DyT(Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else

    U = lambda.*V + DxT(Dx(V)) + DyT(Dy(V));
    u = U(:);

  end

return


%========================================================
%========================================================

function u = IDTWD_lambdaMtx_Fid(v, vDims, WF, WR, lambda)

  V = reshape(v, vDims);

  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
%        U = WF(:,:,k).*lambda(:,:,k).*V(:,:,k) + DxT(WR.*Dx(V(:,:,k))) + ...
%            DyT(WR.*Dy(V(:,:,k)));

      U = WF(:,:,k).*V(:,:,k) + DxT(WR.*Dx(V(:,:,k))) + ...
          DyT(WR.*Dy(V(:,:,k)));

      u(1+(k-1)*N:k*N) = U(:);
    end

  else
    U = WF.*lambda.*V + DxT(WR.*Dx(V)) + DyT(WR.*Dy(V));
    u = U(:);
  end  

return


%========================================================
%========================================================

function u = IWDTWDW_lambdaMtx_Fid(v, vDims, WFN2, WR, lambda)
  
  V = reshape(v, vDims);

  if( length(vDims) == 3 )

    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
%        U = lambda(:,:,k).*V(:,:,k) + ...
      U = V(:,:,k) + ...
          WFN2(:,:,k).*DxT(WR.*Dx(WFN2(:,:,k).*V(:,:,k))) + ...
          WFN2(:,:,k).*DyT(WR.*Dy(WFN2(:,:,k).*V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);

    end

  else
    U = lambda.*V + WFN2.*DxT(WR.*Dx(WFN2.*V)) + ...
        WFN2.*DyT(WR.*Dy(WFN2.*V));
    u = U(:);
  end

return


%========================================================
%========================================================

% Compute KT*K + lambda*D_x^T*D_x + lambda*D_y^T*D_y for vectorised
% image (since called by pcg function)
function u = KTKDTD(v, vDims, K, KT, lambda)
 
  V = reshape(v, vDims);

  KTK_V = KT(K(V));

  if( length(vDims) == 3 )

    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = KTK_V(:,:,k) + lambda*DxT(Dx(V(:,:,k))) + lambda*DyT(Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else
%      U = KT(K(V)) + lambda*DxT(Dx(V)) + lambda*DyT(Dy(V));
    U = KTK_V + lambda*DxT(Dx(V)) + lambda*DyT(Dy(V));
    u = U(:);
  end


return

%========================================================
%========================================================

function u = KTKDTD_lambdaMtx(v, vDims, K, KT, lambda)
 
  V = reshape(v, vDims);

  KTK_V = KT(lambda.*K(V));

  if( length(vDims) == 3 )

    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = KTK_V(:,:,k) + DxT(Dx(V(:,:,k))) + DyT(Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else
    U = KTK_V + DxT(Dx(V)) + DyT(Dy(V));
    u = U(:);
  end


return

%========================================================
%========================================================

% Compute I + lambda*D_x^T*W_R*D_x + lambda*D_y^T*W_R*D_y for vectorised image
% (since called by pcg function) 
function u = IDTWD(v, vDims, WF, WR, lambda)

  V = reshape(v, vDims);

  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = WF(:,:,k).*V(:,:,k) + lambda*DxT(WR.*Dx(V(:,:,k))) + ...
          lambda*DyT(WR.*Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else
    U = WF.*V + lambda*DxT(WR.*Dx(V)) + lambda*DyT(WR.*Dy(V));
    u = U(:);
  end  

return

%========================================================
%========================================================

function u = IDTWD_lambdaMtx(v, vDims, WF, WR, lambda)

  V = reshape(v, vDims);

  if( length(vDims) == 3 ) % only two possible values: 
                           % 3 (vector) or 2 (scalar)
    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = WF(:,:,k).*V(:,:,k) + lambda(:,:,k).*DxT(WR.*Dx(V(:,:,k))) + ...
          lambda(:,:,k).*DyT(WR.*Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else
    U = WF.*V + lambda.*DxT(WR.*Dx(V)) + lambda.*DyT(WR.*Dy(V));
    u = U(:);
  end  

return


%========================================================
%========================================================

% Compute I + lambda*W_F^(-1/2)*D_x^T*W_R*D_x*W_F^(-1/2) +
% lambda*W_F^(-1/2)*D_y^T*W_R*D_y*W_F^(-1/2) for vectorised image
% (since called by pcg function) 
function u = IWDTWDW(v, vDims, WFN2, WR, lambda)
  
  V = reshape(v, vDims);
%    WFN2 = WF.^(-1/2);

  if( length(vDims) == 3 )

    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)

      U = V(:,:,k) + ...
          lambda*WFN2(:,:,k).*DxT(WR.*Dx(WFN2(:,:,k).*V(:,:,k))) + ...
          lambda*WFN2(:,:,k).*DyT(WR.*Dy(WFN2(:,:,k).*V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);

    end

  else
    U = V + lambda*WFN2.*DxT(WR.*Dx(WFN2.*V)) + ...
        lambda*WFN2.*DyT(WR.*Dy(WFN2.*V));
    u = U(:);
  end

return

%========================================================
%========================================================

% Compute I + lambda*W_F^(-1/2)*D_x^T*W_R*D_x*W_F^(-1/2) +
% lambda*W_F^(-1/2)*D_y^T*W_R*D_y*W_F^(-1/2) for vectorised image
% (since called by pcg function) 
function u = IWDTWDW_lambdaMtx(v, vDims, WFN2, WR, lambda)
  
  V = reshape(v, vDims);

  if( length(vDims) == 3 )

    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = V(:,:,k) + ...
          lambda(:,:,k).*WFN2(:,:,k).*DxT(WR.*Dx(WFN2(:,:,k).*V(:,:,k))) + ...
          lambda(:,:,k).*WFN2(:,:,k).*DyT(WR.*Dy(WFN2(:,:,k).*V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);

    end

  else
    U = V + lambda.*WFN2.*DxT(WR.*Dx(WFN2.*V)) + ...
        lambda.*WFN2.*DyT(WR.*Dy(WFN2.*V));
    u = U(:);
  end

return


%========================================================
%========================================================

% Compute KT*WF*K + lambda*D_x^T*W_R*D_x + lambda*D_y^T*W_R*D_y for
% vectorised image (since called by pcg function) 
function u = KTWKDTWD(v, vDims, K, KT, WF, WR, lambda)
  
  V = reshape(v, vDims);

  KTWK_V = KT(WF.*K(V));

  if( length(vDims) == 3 )

    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = KTWK_V(:,:,k) + lambda*DxT(WR.*Dx(V(:,:,k))) + ...
          lambda*DyT(WR.*Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else
%      U = KT(WF.*K(V)) + lambda*DxT(WR.*Dx(V)) + lambda*DyT(WR.*Dy(V));
    U = KTWK_V + lambda*DxT(WR.*Dx(V)) + lambda*DyT(WR.*Dy(V));
    u = U(:);
  end

return

%========================================================
%========================================================

% Compute KT*WF*K + lambda*D_x^T*W_R*D_x + lambda*D_y^T*W_R*D_y for
% vectorised image (since called by pcg function) 
function u = KTWKDTWD_lambdaMtx(v, vDims, K, KT, WF, WR, lambda)
  
  V = reshape(v, vDims);

  KTWK_V = KT(lambda.*WF.*K(V));

  if( length(vDims) == 3 )

    N = vDims(1)*vDims(2);
    u = zeros( size(v) );

    for k=1:vDims(3)
      U = KTWK_V(:,:,k) + DxT(WR.*Dx(V(:,:,k))) + ...
          DyT(WR.*Dy(V(:,:,k)));
      u(1+(k-1)*N:k*N) = U(:);
    end

  else
    U = KTWK_V + DxT(WR.*Dx(V)) + DyT(WR.*Dy(V));
    u = U(:);
  end

return

%========================================================
%========================================================

% Compute relative residual as used by pcg function
function y = pcgrelres(A,b,x)
  y = norm(A(x) - b)/norm(b);
return

%========================================================
%========================================================

% Compute generalised TV functional value
function [fnc, df, reg] = gentvfnctnl(S, K, lambda, p, q, U)

  sDims = size(S);
  s = S(:);
  u = U(:);

  if isempty(K),
    df = sum(abs(u-s).^p)/p;
  else
    df = sum(vec(abs(K(U)-S).^p))/p;
  end


  if(length(sDims) == 3) 
    dxu = Dx(U(:,:,1));
    dyu = Dy(U(:,:,1));
    tmp = dxu.^2 + dyu.^2;
    for k=2:sDims(3)
      dxu = Dx(U(:,:,k));
      dyu = Dy(U(:,:,k));
      tmp = tmp + dxu.^2 + dyu.^2;
    end
    reg = sum(abs(sqrt(tmp(:))).^q)/q;
  else
    dxu = Dx(U);
    dyu = Dy(U);
    reg = sum(abs(sqrt(vec(dxu.^2 + dyu.^2))).^q)/q;
  end
  fnc = df + lambda*reg;

return
  
%========================================================
%========================================================

% Compute weighted approximation to generalised TV functional value
function [fnc, df, reg] = gentvwtfnctnl(S, K, lambda, p, q, U, WF, WR)
  
  sDims = size(S);
  s = S(:);
  u = U(:);

  if isempty(K),
    df = sum(abs(vec(WF.^(1/2)).*(u-s)).^2)/2;
  else
    df = sum(abs(vec(WF.^(1/2)).*vec(K(U)-S)).^2)/2;
  end

  if(length(sDims) == 3) 
    dxu = Dx(U(:,:,1));
    dyu = Dy(U(:,:,1));
    tmp = vec((WR.^(1/2)).*dxu(:,:,1)).^2 + vec((WR.^(1/2)).*dyu(:,:,1)).^2;
    for k=2:sDims(3)
      dxu = Dx(U(:,:,k));
      dyu = Dy(U(:,:,k));
      tmp = tmp + vec((WR.^(1/2)).*dxu(:,:,1)).^2 + ...
            vec((WR.^(1/2)).*dyu(:,:,1)).^2;
    end
    reg = sum(abs( tmp ))/2;
  else
    dxu = Dx(U);
    dyu = Dy(U);
    reg = sum(abs([vec((WR.^(1/2)).*dxu); vec((WR.^(1/2)).*dyu)]).^2)/2;
  end

  fnc = df + lambda(:).*reg; % FIXME: lambda --> matrix
  
return

%========================================================
%========================================================

% Compute weighted approximation to generalised TV functional value
% for denoising with substitution (indirect) method
function [fnc, df, reg] = gentvdnwtfnctnl(S, lambda, p, q, V, WFN2, WR)

 
  sDims = size(S);
  s = S(:);
  v = V(:);

  if( length(sDims) == 3 )      % color/vectorial image

    regV = ( (WR.^(1/2)).*Dx(WFN2(:,:,1).*V(:,:,1)) ).^2 + ...
    ( (WR.^(1/2)).*Dy(WFN2(:,:,1).*V(:,:,1)) ).^2;

    for k=2:sDims(3)
      regV = regV + ( (WR.^(1/2)).*Dx(WFN2(:,:,k).*V(:,:,k)) ).^2 + ...
      ( (WR.^(1/2)).*Dy(WFN2(:,:,k).*V(:,:,k)) ).^2;
    end

  else

    regV = ( (WR.^(1/2)).*Dx(WFN2.*V) ).^2 + ( (WR.^(1/2)).*Dy(WFN2.*V) ).^2;
  
  end

  reg = sum(regV(:));
  df = sum(abs(v - (s./vec(WFN2))).^2)/2;
  fnc = df + 0.5*lambda(:).*reg;   % FIXME: lambda --> Matrix

return


% Open iteration status display
function iterdispopen(dspflag)
  if dspflag,
    disp('Itn  Fnc        DFid       Reg        WFnc       WDFid      WReg       RelRes    CGIt  Flg Time');
    disp('--------------------------------------------------------------------------------------------------');
  end
return

%========================================================
%========================================================

% Compute and display iteration status
function is = iterdisp(S, K, lambda, p, q, U, WF, WR, k, relres, pcgflg, pcgit, sbstflg, dspflag)

  t = toc;
  [tnp, dfnp, rnp] = gentvfnctnl(S, K, lambda, p, q, U);
  if isempty(K) && sbstflg,
    V = U./WF;
    [tnw, dfnw, rnw] = gentvdnwtfnctnl(S, lambda, p, q, V, WF, WR);
  else
    [tnw, dfnw, rnw] = gentvwtfnctnl(S, K, lambda, p, q, U, WF, WR);
  end
  if dspflag,
    kstr = '   ';
    if ~isempty(k),
      kstr = sprintf('%3d', k);
    end
    pistr = '     ';
    if ~isempty(pcgit),
      pistr = sprintf('%5d', pcgit);
    end
    pfstr = '  ';
    if ~isempty(pcgflg),
      pfstr = sprintf('%2d', pcgflg);
    end
    disp(sprintf('%s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %s %s %8.1e', ...
                 kstr, tnp, dfnp, rnp, tnw, dfnw, rnw, relres, pistr, pfstr, t));
  end
  is = [tnp dfnp rnp tnw dfnw rnw relres pcgit pcgflg t];

return

%========================================================
%========================================================

% Open iteration status display
function iterdispopen_nqp(dspflag)
  if dspflag,
    disp('Itn  Fnc        DFid       Reg         NQPIt      NQPeps (target)         Time');
    disp('--------------------------------------------------------------------------------');
  end
return

%========================================================
%========================================================

function is = iterdisp_nqp(S, K, lambda, p, q, U, WF, WR, k, NQPit, local_eps, NQPeps, dspflag)

  t = toc;

  if( isempty(WF) ) WF = 1; end

  [tnp, dfnp, rnp] = gentvfnctnl(S, K, lambda, p, q, U);

  if dspflag,
    kstr = '   ';
    if ~isempty(k),
      kstr = sprintf('%3d', k);
    end

    nqpitstr = '    ';
    if ~isempty(NQPit),
      nqpitstr = sprintf('%4d', NQPit);
    end

    nqpeps_str = '     ';
    if ~isempty(NQPeps),
      nqpeps_str = sprintf('%10.3e', NQPeps);
    end

    localeps_str = '     ';
    if ~isempty(local_eps),
      localeps_str = sprintf('%10.3e', local_eps);
    end

    disp(sprintf('%s %10.3e %10.3e %10.3e  %s    %s (%s) %8.1e', ...
                 kstr, tnp, dfnp, rnp, nqpitstr, nqpeps_str, localeps_str, t));
  end
  is = [tnp dfnp rnp NQPit NQPeps local_eps t];

return

%========================================================
%========================================================

function is = iterdisp_poisson(S, K, lambda, p, q, U, WF, WR, k, NQPit, local_eps, NQPeps, dspflag)

  t = toc;

  if( isempty(WF) ) WF = 1; end

%    [tnp, dfnp, rnp] = gentvfnctnl(S, K, lambda, p, q, U);

  tnp = 0;
  dfnp = 0;
  rnp = 0;
  if( isempty(K) ) KU = U; else KU = K(U); end;

  sDims = size(S);

  if(length(sDims) == 3) 
    vKU = KU(:,:,1);
    vS  = S(:,:,1);
    F = sum( vKU(:) - vS(:).*log( (vKU(:)>0).*vKU(:) + (vKU(:)<=0) ) );
    for k=2:sDims(3)
      vKU = KU(:,:,k);
      vS  = S(:,:,k);
      F = F + sum( vKU(:) - vS(:).*log( (vKU(:)>0).*vKU(:) + (vKU(:)<=0) ) );
    end
  else
    F = sum( KU(:) - S(:).*log( (KU(:)>0).*KU(:) + (KU(:)<=0) ) );
  end

%------------------
  if(length(sDims) == 3) 
    dxu = Dx(U(:,:,1));
    dyu = Dy(U(:,:,1));
    tmp = dxu.^2 + dyu.^2;
    for k=2:sDims(3)
      dxu = Dx(U(:,:,k));
      dyu = Dy(U(:,:,k));
      tmp = tmp + dxu.^2 + dyu.^2;
    end
    R = sum(abs(sqrt(tmp(:))).^q)/q;
  else
    dxu = Dx(U);
    dyu = Dy(U);
    R = sum(abs(sqrt(vec(dxu.^2 + dyu.^2))).^q)/q;
  end

  tnp  = F + R;
  dfnp = F;
  rnp  = R;

%------------------

  if dspflag,
    kstr = '   ';
    if ~isempty(k),
      kstr = sprintf('%3d', k);
    end

    nqpitstr = '    ';
    if ~isempty(NQPit),
      nqpitstr = sprintf('%4d', NQPit);
    end

    nqpeps_str = '     ';
    if ~isempty(NQPeps),
      nqpeps_str = sprintf('%10.3e', NQPeps);
    end

    localeps_str = '     ';
    if ~isempty(local_eps),
      localeps_str = sprintf('%10.3e', local_eps);
    end

    disp(sprintf('%s %10.3e %10.3e %10.3e  %s    %s (%s) %8.1e', ...
                 kstr, tnp, dfnp, rnp, nqpitstr, nqpeps_str, localeps_str, t));
  end
  is = [tnp dfnp rnp NQPit NQPeps local_eps t];

return

%========================================================
%========================================================

% Close iteration status display
function iterdispclose(dspflag)
  if dspflag,
    disp('--------------------------------------------------------------------------------------------------');
  end
return


%========================================================
%========================================================


