function [y,nz] = ompdenoise3(params,msgdelta)
%OMPDENOISE3 OMP denoising of 3-D signals.
%  OMPDENOISE3 denoises a 3-dimensional signal using OMP denoising. The
%  function syntax is identical to OMPDENOISE, but it runs significantly
%  faster on 3-D signals. OMPDENOISE3 requires somewhat more memory than
%  OMPDENOISE (approximately the size of the input signal), so if memory is
%  limited, OMPDENOISE can be used instead.
%
%  See also OMPDENOISE.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% parse input arguments %
addpath('ompbox');

x = params.x;
D = params.dict;
blocksize = params.blocksize;


% blocksize %
if (numel(blocksize)==1)
  blocksize = ones(1,3)*blocksize;
end


% maxval %
if (isfield(params,'maxval'))
  maxval = params.maxval;
else
  maxval = 1;
end


% gain %
if (isfield(params,'gain'))
  gain = params.gain;
else
  gain = 1.15;
end


% maxatoms %
if (isfield(params,'maxatoms'))
  maxatoms = params.maxatoms;
else
  maxatoms = floor(prod(blocksize)/2);
end


% stepsize %
if (isfield(params,'stepsize'))
  stepsize = params.stepsize;
  if (numel(stepsize)==1)
    stepsize = ones(1,3)*stepsize;
  end
else
  stepsize = ones(1,3);
end
if (any(stepsize<1))
  error('Invalid step size.');
end


% noise mode %
if (isfield(params,'noisemode'))
  switch lower(params.noisemode)
    case 'psnr'
      sigma = maxval / 10^(params.psnr/20);
    case 'sigma'
      sigma = params.sigma;
    otherwise
      error('Invalid noise mode specified');
  end
elseif (isfield(params,'sigma'))
  sigma = params.sigma;
elseif (isfield(params,'psnr'))
  sigma = maxval / 10^(params.psnr/20);
else
  error('Noise strength not specified');
end


% lambda %
if (isfield(params,'lambda'))
  lambda = params.lambda;
else
  lambda = maxval/(10*sigma);
end


% msgdelta %
if (nargin <2)
  msgdelta = 5;
end
if (msgdelta<=0)
  msgdelta = -1;
end


epsilon = sqrt(prod(blocksize)) * sigma * gain;   % target error for omp


MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end


% denoise the signal %

if (memusage >= MEM_NORMAL)
  G = D'*D;
end

nz = 0;  % count non-zeros in block representations

% the denoised signal
y = zeros(size(x));

% indices of the current block
ids1 = 1:blocksize(1);
ids2 = 1:blocksize(2);
ids3 = 1:blocksize(3);

% the current batch of blocks: contains as columns the blocks from the
% sub-volume x(:,ids2,ids3)
blocks = zeros(prod(blocksize), floor((size(x,1)-blocksize(1))./stepsize(1)) + 1);
for j = 1:blocksize(3)
  allblocks = im2col(x(:,ids2,ids3(j)),blocksize(1:2),'sliding');
  blocks((j-1)*blocksize(1)*blocksize(2)+1 : j*blocksize(1)*blocksize(2), :) = allblocks(:,1:stepsize(1):end);
end

% remove DC
[blocks, dc] = remove_dc(blocks,'columns');

% denoise the blocks
if (memusage == MEM_LOW)
  gamma = omp2(D,blocks,[],epsilon,'maxatoms',maxatoms);
else
  gamma = omp2(D'*blocks,sum(blocks.*blocks),G,epsilon,'maxatoms',maxatoms);
end
nz = nz + nnz(gamma);
cleanblocks = add_dc(D*gamma, dc, 'columns');

k = 1;   % the index of the current block in the blocks matrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  sum the cleaned blocks into y
%  every time the current batch of signals is exhausted, extract the next
%  batch of signals and denoise them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% lastids contains the indices of the last block in each dimension
lastids = stepsize .* floor((size(x)-blocksize)./stepsize) + 1;

blocknum = prod(floor((size(x)-blocksize)./stepsize) + 1);
tid = timerinit('ompdenoise', blocknum); lastmsgcheck = 0;
for i = 1 : blocknum
  
  % current denoised block
  y(ids1,ids2,ids3) = y(ids1,ids2,ids3) + reshape(cleanblocks(:,k),blocksize);

  % increment block ids, and compute next block of signals if necessary
  if (i<blocknum)
    
    if (ids1(1) == lastids(1))
      ids1 = 1:blocksize(1);
      if (ids2(1) == lastids(2))
        ids2 = 1:blocksize(2);
        ids3 = ids3+stepsize(3);
      else
        ids2 = ids2+stepsize(2);
      end
      
      % current batch exhausted, compute next batch
      for j = 1:blocksize(3)
        allblocks = im2col(x(:,ids2,ids3(j)),blocksize(1:2),'sliding');
        blocks((j-1)*blocksize(1)*blocksize(2)+1 : j*blocksize(1)*blocksize(2),:) = allblocks(:,1:stepsize(1):end);
      end
      [blocks, dc] = remove_dc(blocks,'columns');
      if (memusage == MEM_LOW)
        gamma = omp2(D,blocks,[],epsilon,'maxatoms',maxatoms);
      else
        gamma = omp2(D'*blocks,sum(blocks.*blocks),G,epsilon,'maxatoms',maxatoms);
      end
      nz = nz + nnz(gamma);
      cleanblocks = add_dc(D*gamma, dc, 'columns');

      k = 1;
      
    else
      
      % just increment block ids and DC value
      ids1 = ids1+stepsize(1);
      k = k+1;
      
    end
  end
   
  % display status
  if (msgdelta>0 && i-lastmsgcheck>=30)
    lastmsgcheck=i;
    timereta(tid, i, msgdelta);
  end
  
end

if (msgdelta>0)
  timereta(tid, i);
end

nz = nz/blocknum;  % average number of non-zeros

% average the denoised and noisy signals
cnt = countcover(size(x),blocksize,stepsize);
y = (y+lambda*x)./(cnt + lambda);

