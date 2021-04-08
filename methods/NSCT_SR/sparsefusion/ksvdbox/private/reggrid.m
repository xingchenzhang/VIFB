function [varargout] = reggrid(sz,num,mode)
%REGGRID Regular sampling grid.
%  [I1,I2,...,Ip] = REGGRID([N1 N2 ... Np], NUM) returns the indices
%  of a regular uniform sampling grid over a p-dimensional matrix with
%  dimensions N1xN2x...xNp. NUM is the minimal number of required samples,
%  and it is ensured that the actual number of samples, given by
%  length(I1)xlength(I2)x...xlength(Ip), is at least as large as NUM.
%
%  [I1,I2,...,Ip] = REGGRID([N1 N2 ... Np], NUM,'MODE') specifies the
%  method for distributing the samples along each dimension. Valid modes
%  include 'eqdist' (the default mode) and 'eqnum'. 'eqdist' indicates an
%  equal distance between the samples in each dimension, while 'eqnum'
%  indicates an equal number of samples in each dimension.
%
%  Notes about MODE:
%
%    1. The 'eqnum' mode will generally fail when the p-th root of NUM
%    (i.e. NUM^(1/p)) is larger than min([N1 N2 ... Np]). Thus 'eqdist' is
%    the more useful choice for sampling an arbitrary number of samples
%    from the matrix (up to the total number of matrix entries).
%  
%    2. In both modes, the equality (of the distance between samples, or
%    the number of samples in each dimension) is only approximate. This is
%    because REGGRID attempts to maintain the appropriate equality while at
%    the same time find a sampling pattern where the total number of
%    samples is as close as possible to NUM. In general, the larger {Ni}
%    and NUM are, the tighter the equality.
%
%  Example: Sample a set of blocks uniformly from a 2D image.
%
%    n = 512; blocknum = 20000; blocksize = [8 8];
%    im = rand(n,n);
%    [i1,i2] = reggrid(size(im)-blocksize+1, blocknum);
%    blocks = sampgrid(im, blocksize, i1, i2);
%
%  See also SAMPGRID.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  November 2007

dim = length(sz);

if (nargin<3)
  mode = 'eqdist';
end

if (any(sz<1))
  error(['Invalid matrix size : [' num2str(sz) ']']);
end

if (num > prod(sz))
  warning(['Invalid number of samples, returning maximum number of samples.']);
elseif (num <= 0)
  if (num < 0)
    warning('Invalid number of samples, assuming 0 samples.');
  end
  for i = 1:length(sz)
    varargout{i} = [];
  end
  return;
end


if (strcmp(mode,'eqdist'))
  
  % approximate distance between samples: total volume divided by number of
  % samples gives the average volume per sample. then, taking the p-th root
  % gives the average distance between samples
  d = (prod(sz)/num)^(1/dim);
  
  % compute the initial guess for number of samples in each dimension.
  % then, while total number of samples is too large, decrese the number of
  % samples by one in the dimension where the samples are the most crowded.
  % finally, do the opposite process until just passing num, so the final
  % number of samples is the closest to num from above.
  
  n = min(max(round(sz/d),1),sz);   % set n so that it saturates at 1 and sz
  
  active_dims = find(n>1);    % dimensions where the sample num can be reduced
  while(prod(n)>num && ~isempty(active_dims))
    [y,id] = min((sz(active_dims)-1)./n(active_dims));
    n(active_dims(id)) = n(active_dims(id))-1;
    if (n(active_dims(id)) < 2)
      active_dims = find(n>1);
    end
  end

  active_dims = find(n<sz);    % dimensions where the sample num can be increased
  while(prod(n)<num && ~isempty(active_dims))
    [y,id] = max((sz(active_dims)-1)./n(active_dims));
    n(active_dims(id)) = n(active_dims(id))+1;
    if (n(active_dims(id)) >= sz(active_dims(id)))
      active_dims = find(n<sz);
    end
  end

  for i = 1:dim
    varargout{i} = round((1:n(i))/n(i)*sz(i));
    varargout{i} = varargout{i} - floor((varargout{i}(1)-1)/2);
  end
  
elseif (strcmp(mode,'eqnum'))
  
  % same idea as above
  n = min(max( ones(size(sz)) * round(num^(1/dim)) ,1),sz);

  active_dims = find(n>1);
  while(prod(n)>num && ~isempty(active_dims))
    [y,id] = min((sz(active_dims)-1)./n(active_dims));
    n(active_dims(id)) = n(active_dims(id))-1;
    if (n(active_dims(id)) < 2)
      active_dims = find(n>1);
    end
  end
  
  active_dims = find(n<sz);
  while(prod(n)<num && ~isempty(active_dims))
    [y,id] = max((sz(active_dims)-1)./n(active_dims));
    n(active_dims(id)) = n(active_dims(id))+1;
    if (n(active_dims(id)) >= sz(active_dims(id)))
      active_dims = find(n<sz);
    end
  end
  
  for i = 1:dim
    varargout{i} = round((1:n(i))/n(i)*sz(i));
    varargout{i} = varargout{i} - floor((varargout{i}(1)-1)/2);
  end
else
  error('Invalid sampling mode');
end

