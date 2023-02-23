function y = sampgrid(x,blocksize,varargin)
%SAMPGRID Sample a multi-dimensional matrix on a regular grid.
%  Y = SAMPGRID(X,BLOCKSIZE,I1,I2,...,Ip) extracts block samples of size
%  BLOCKSIZE from the p-dimensional matrix X, arranging the samples as the
%  column vectors of the matrix Y. The locations of the (1,1,..,1)-th
%  elements of each block are given in the index vectors I1,I2,..Ip. The
%  total number of samples taken is length(I1)xlength(I2)x...xlength(Ip).
%  BLOCKSIZE should either be a p-element vector of the form [N1,N2,...Np],
%  or a scalar N which is shorthand for the square block size [N N ... N].
%
%  Example: Sample a set of blocks uniformly from a 2D image.
%
%    n = 512; blocknum = 20000; blocksize = [8 8];
%    im = rand(n,n);
%    [i1,i2] = reggrid(size(im)-blocksize+1, blocknum);
%    blocks = sampgrid(im, blocksize, i1, i2);
%
%  See also REGGRID.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  November 2007


p = ndims(x);
if (p==2 && any(size(x)==1) && length(blocksize)==1)
  p = 1;
end

if (numel(blocksize)==1)
  blocksize = ones(1,p)*blocksize;
end

n = zeros(1,p);
for i = 1:p
  n(i) = length(varargin{i});
end

nsamps = prod(n);

% create y of the same class as x
y = zeros(prod(blocksize),nsamps,class(x));

% ids() contains the index of the current block in I1..Ip
ids = ones(p,1);

% block_ids contains the indices of the current block in X
block_ids = cell(p,1);
for j = 1:p
  block_ids{j} = varargin{j}(1) : varargin{j}(1)+blocksize(j)-1;
end

for k = 1:nsamps
  block = x(block_ids{:});
  y(:,k) = block(:);
  
  % increment ids() and block_ids{}
  if (k<nsamps)
    j = 1;
    while (ids(j) == n(j))
      ids(j) = 1;
      block_ids{j} = varargin{j}(1) : varargin{j}(1)+blocksize(j)-1;
      j = j+1;
    end
    ids(j) = ids(j)+1;
    block_ids{j} = varargin{j}(ids(j)) : varargin{j}(ids(j))+blocksize(j)-1;
  end
end

