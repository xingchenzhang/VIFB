function [U, optinf] = sbl1tv2(S, FwdOp, lambda, opt)

% sbl1tv2 -- Split-Bregman method for l1-TV in 2D
%            The problem
%              min_u |A u-s|_1 + lambda*|sqrt((D_x u)^2 + (D_y u)^2)|_1
%            is solved with additional parameters alpha and beta
%            for the l2 terms (resulting from the split) on the
%            data fidelity and regularisation terms
%            respectively. A further parameter gamma controls
%            another split for a bound constraint on the result.
%          
% Usage:
%       [U, optinf] = sbl1tv2(S, FwdOp, lambda, opt);
%
% Input:
%       S           Input image
%       FwdOp       Cell array {A, AT} of function handles for forward 
%                   linear operator and its transpose. These operators 
%                   act on images rather than vectors, but the transpose
%                   should be considered to be in the sense of the 
%                   operator on a vectorised image. Set to [] for a 
%                   denoising problem.
%       lambda      Regularisation term weighting factor
%       opt         Structure providing algorithm options (see code
%                   for details)
% Output:
%       U           Output image
%       optinf      Optimisation progress information

  
  if nargin < 4,
    opt = [];
  end
  opt = defaultopts(opt);
  if ~isfield(opt,'alpha'), opt.alpha = 1; end
  if ~isfield(opt,'beta'), opt.beta = lambda; end
  if ~isfield(opt,'gamma'), opt.gamma = opt.alpha/100; end
  alpha = opt.alpha;
  beta = opt.beta;
  gamma = opt.gamma;

  % Check argument FwdOp
  if isempty(FwdOp),
    [Nr, Nc] = size(S);
    A = @(x) x;
    AT = @(x) x;
  elseif iscell(FwdOp) && isa(FwdOp{1}, 'function_handle') && ...
                          isa(FwdOp{2}, 'function_handle'),
    A = FwdOp{1};
    AT = FwdOp{2};
    [Nr, Nc] = size(AT(S));
  end

  optinf = [];
  optinf.itstat = [];
  optinf.opt = opt;
  J0 = [];
  ffc = Inf;

  hstr1 = 'Itn  Fnc       DFid      TV(u)       B(v)    B(x)    B(y)';
  hstr2 = '    B(z)     RelR    CGIt F';
  nsepc = 85;
  sfmst = '%3d %9.2e %9.2e %9.2e    %7.1e %7.1e %7.1e %7.1e  %7.1e %4d %1d';
  
  if opt.Verbose && opt.MaxMainIter > 0,
    disp([hstr1 hstr2]);
    disp(char('-' * ones(1,nsepc)));
  end
  
  U = zeros(Nr,Nc);
  V = zeros(size(S));
  X = Dx(U);
  Y = Dy(U);
  Z = U;
  Bv = zeros(size(V));
  Bx = zeros(size(X));
  By = zeros(size(Y));
  Bz = zeros(size(Z));

  tstart = tic;

  k = 1;
  while k <= opt.MaxMainIter & ffc > opt.FracStopTol,

    ulhs = @(x) ATADTD(x,Nr,Nc,A,AT,alpha,beta,gamma);
    urhs = vec(alpha*AT(V + S - Bv) + beta*(DxT(X - Bx) + DyT(Y - By)) ...
           + gamma*(Z - Bz));
    [u, flg, rlr, pit] = pcg(ulhs, urhs, opt.CGTol, opt.MaxCGIter);
    U = reshape(u,Nr,Nc);
    AU = A(U);
 
    V = shrink1(AU - S + Bv, 1/alpha);
    
    [X, Y] = shrink2(Dx(U) + Bx, Dy(U) +  By, lambda/beta);

    Z = U + Bz;
    Z(Z > opt.BoundRange(2)) = opt.BoundRange(2);
    Z(Z < opt.BoundRange(1)) = opt.BoundRange(1);

    Bv = Bv + AU - S - V;
    Bx = Bx + Dx(U) - X;
    By = By + Dy(U) - Y;
    Bz = Bz + U - Z;
    
    Jdf = sum(abs(vec(AU-S)));
    Jtv = sum(vec(abs(sqrt(Dx(U).^2 + Dy(U).^2))));
    J1 = Jdf + lambda*Jtv;
    Jv = 0.5*sum(vec(V - (AU - S)).^2);
    Jx = 0.5*sum(vec(X - Dx(U)).^2);
    Jy = 0.5*sum(vec(Y - Dy(U)).^2);
    Jz = 0.5*sum(vec(Z - U).^2);

    optinf.itstat = [optinf.itstat; ...
                     [J1 Jdf Jtv Jv Jx Jy Jz rlr pit flg]];

    if opt.Verbose,
      disp(sprintf(sfmst, k, J1, Jdf, Jtv, Jv, Jx, Jy, Jz, rlr, pit, flg));
    end
    
    if ~isempty(J0),
      ffc = abs((J0-J1)/J1);
    end
    J0 = J1;
    k = k + 1;
    
  end

  optinf.runtime = toc(tstart);

  if opt.Verbose && opt.MaxMainIter > 0,
    disp(char('-' * ones(1,nsepc)));
  end
  
  optinf.V = V;
  optinf.X = X;
  optinf.Y = Y;
  optinf.Z = Z;
  optinf.Bv = Bv;
  optinf.Bx = Bx;
  optinf.By = By;
  optinf.Bz = Bz;

return


function u = shrink1(v, lambda)
    
  u = sign(v).*max(0, abs(v) - lambda);
  
return 


function [u, v] = shrink2(s, t, lambda)

  a = sqrt(s.^2 + t.^2);
  b = max(0, a - lambda);
  
  b(a == 0) = 0;
  a(a == 0) = 1;
  
  u = b.*s./a;
  v = b.*t./a;

return 


% Vectorise image v
function u = vec(v)

  u = v(:);
  
return


% Consider operator Dy to be the y derivate of a vectorised
% image. Apply this operator to the unvectorised image and return
% a gradient image.
function u = Dy(v)

  u = [diff(v); zeros(1,size(v,2))];
    
return


% Consider operator Dy to be the y derivate of a vectorised
% image. Apply the transpose of this operator to the unvectorised
% image and return a gradient image.
function u = DyT(v)

  u0 = -v(1,:);
  u1 = -diff(v);
  u2 = v(end-1,:);
  u = [u0; u1(1:(end-1),:); u2];
    
return


% Consider operator Dx to be the x derivate of a vectorised
% image. Apply this operator to the unvectorised image and return
% a gradient image.
function u = Dx(v)

  u = [diff(v,1,2) zeros(size(v,1),1)];
    
return


% Consider operator Dx to be the x derivate of a vectorised
% image. Apply the transpose of this operator to the unvectorised
% image and return a gradient image.
function u = DxT(v)

  u = DyT(v')';

return


function u = ATADTD(v, Nr, Nc, A, AT, alpha, beta, gamma)

  V = reshape(v, Nr, Nc);
  U = alpha*AT(A(V)) + beta*DxT(Dx(V)) + beta*DyT(Dy(V)) + gamma*V;
  u = U(:);

return


function opt = defaultopts(opt)

  if ~isfield(opt,'Verbose'),
    opt.Verbose = 0;
  end
  if ~isfield(opt,'MaxMainIter'),
    opt.MaxMainIter = 20;
  end
  if ~isfield(opt,'FracStopTol'),
    opt.FracStopTol = 1e-3;
  end
  if ~isfield(opt,'MaxCGIter'),
    opt.MaxCGIter = 100000;
  end
  if ~isfield(opt,'CGTol'),
    opt.CGTol = 1e-5;
  end
  if ~isfield(opt,'BoundRange'),
    opt.BoundRange = [-Inf Inf];
    opt.gamma = 0;
  end

return
