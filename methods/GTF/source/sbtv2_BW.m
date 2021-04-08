function [U, optinf] = sbtv2(S, KC, lambda, gamma, nit, U0)

  if nargin < 6,
    U0 = [];
  end
  if nargin < 5,
    nit = 10;
  end
  if nargin < 4,
    gamma = lambda; 
  end

  dflag = 1;
  if nargout > 1,
    %dflag = 0;
  end

  % Check argument KC
  if isempty(KC),
    [Nr, Nc] = size(S);
    A = @(x) x;
    AT = @(x) x;
  else
    if iscell(KC) && isa(KC{1}, 'function_handle') && ...
            isa(KC{2}, 'function_handle'),
      A = KC{1};
      AT = KC{2};
      [Nr, Nc] = size(AT(S));
    else
      warning('argument KC must be empty or a cell array {K,KT} of function handles');
      U = [];
      return;
    end
  end

  s = S(:);
  if isempty(U0),
    U = zeros(Nr,Nc);
  else
    U = U0;
  end
  u = U(:);

  %imagesc(U); drawnow;
  
  tstart = tic;

  dx = Dx(U); dy = Dy(U);
  %[dx, dy] = shrink(Dx(U), Dy(U), lambda/gamma);
  bx = zeros(size(dx));
  by = zeros(size(dy));


  its = [];
  
  ats = AT(S);

  for k = 1:nit,

    rhs = vec(ats + gamma*(DxT(dx - bx) + DyT(dy - by)));
    [u,flg,rlr,pit] = pcg(@(x) ATADTD(x,Nr,Nc,A,AT,gamma), rhs, ...
                          1e-4, 1000, [], [], u);
    U = reshape(u,Nr,Nc);
    [dx, dy] = shrink(Dx(U) + bx, Dy(U) +  by, lambda/gamma);

    bx = bx + Dx(U) - dx;
    by = by + Dy(U) - dy;

    df = 0.5*norm(vec(A(U)-S))^2;
    rg = sum(vec(abs(sqrt(Dx(U).^2 + Dy(U).^2))));
    spx = 0.5*norm(vec(dx - Dx(U)))^2;
    spy = 0.5*norm(vec(dy - Dy(U)))^2;

    t = toc(tstart);
    its = [its; [k t df+lambda*rg df rg spx spy ]];
    
    if(dflag)
      disp(sprintf('%3d  %.3e %.3e %.3e   %.3e %.3e   %7.1e %4d %1d', ...
                 k, df + lambda*rg, df, rg, spx, spy, rlr, pit, flg)); 
    end

  end

  optinf.dx = dx;
  optinf.dy = dy;
  optinf.bx = bx;
  optinf.by = by;
  optinf.itstat = its;
  
return


function [u, v] = shrink(s, t, lambda)

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


function u = ATADTD(v, Nr, Nc, A, AT, gamma)

  V = reshape(v, Nr, Nc);
  U = AT(A(V)) + gamma*DxT(Dx(V)) + gamma*DyT(Dy(V));
  u = U(:);
 
return
