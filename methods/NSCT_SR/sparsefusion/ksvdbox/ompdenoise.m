function [y,nz] = ompdenoise(params,msgdelta)
%OMPDENOISE OMP denoising.
%  Y = OMPDENOISE(PARAMS) denoises a (possibly multi-dimensional) signal
%  using OMP denoising. OMPDENOISE first denoises each block of the noisy
%  signal individually, then combines the denoised blocks, averaging
%  overlapping blocks. A single block is denoised by removing its mean,
%  applying OMP, and restoring its mean.
%
%  Y = OMPDENOISE(PARAMS,MSGDELTA) specifies the frequency of message
%  printing during the process. MSGDELTA should be a positive number
%  representing the interval in seconds between messages. A zero or
%  negative value cancels these messages. Default is MSGDELTA=5.
%
%  [Y,NZ] = OMPDENOISE(...) also returns the average number of non-zero
%  coefficients in the representations of the denoised blocks.
%
%
%  Required fields in PARAMS:
%  --------------------------
%
%    'x' - Noisy signal.
%      The signal to denoise (can be multi-dimensional). Should be of type
%      double, and (for PSNR computations) with values within [0,1] (to
%      specify a different range, see parameter 'maxval' below).
%
%    'blocksize' - Size of block.
%      Indicates the size of the blocks to operate on. Should be either an
%      array of the form [N1 N2 ... Np], where p is the number of
%      dimensions of x, or simply a scalar N, representing the square block
%      [N N ... N]. See parameter 'stepsize' below to specify the amount of
%      overlap between the blocks.
%
%    'dict' - The dictionary.
%      Specifies the dictionary for denoising each block in x. Should be of
%      size NxL where N=N1*N2*...*Np and L is the number of atoms.
%
%    'psnr' / 'sigma' - Noise power.
%      Specifies the noise power in dB (psnr) or the noise standard
%      deviation (sigma), used to determine the target error for
%      sparse-coding each block. If both fields are present, sigma is used
%      unless the field 'noisemode' is specified (below). When specifying
%      the noise power in psnr, make sure to set the 'maxval' parameter
%      as well (below) if the signal values are not within [0,1].
%
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'stepsize' -  Interval between neighboring blocks.
%      Specifies the interval (in pixels/voxels) between neighboring blocks
%      to denoise. By default, all overlapping blocks are denoised and
%      averaged. This can be changed by specifying an alternate stepsize,
%      as an array of the form [S1 S2 ... Sp] (where p is the number of
%      dimensions of x). This sets the distance between neighboring blocks
%      to be Si in the i-th dimension. Stepsize can also be a scalar S,
%      corresponding to the step size [S S ... S]. Each Si must be >= 1,
%      and, to ensure coverage of the entire noisy signal, size(x,i)-Ni
%      should be a multiple of Si for all i. The default stepsize is 1.
%
%    'maxval' - Maximal intensity value.
%      Specifies the range of the signal values. Used to determine the
%      noise standard deviation when the noise power is specified in psnr.
%      By default, the signal values are assumed to be within [0,1]. When
%      'maxval' is specified, this range changes to [0,maxval].
%
%    'memusage' - Memory usage.
%      This parameter controls memory usage of the function. 'memusage'
%      should be one of the strings 'high', 'normal' (default) or 'low'.
%      When set to 'low', the matrix G=D'*D is not precomputed before
%      calling OMP, which reduces memory consumption when memory resources
%      are low. This, however, will dramatically increase runtime and is
%      not recommended for normal use. The 'high' and 'normal' modes are
%      equivalent settings in this function.
%
%
%  Optional fields in PARAMS - advanced:
%  -------------------------------------
%
%    'noisemode' - Noise power mode.
%      Specifies whether the 'psnr' or 'sigma' fields should be used to
%      determine the noise power. This is useful when both fields are
%      present in PARAMS. 'noisemode' should be one of the string 'psnr' or
%      'sigma'. If it is not present, and both fields are specified,
%      'sigma' is used.
%
%    'gain' - Noise gain. 
%      A positive value (usually near 1) controlling the target error for
%      sparse-coding each block. When gain=1, the target error is precisely
%      the value derived from the psnr/sigma fields. when gain is different
%      from 1, the target error is multiplied by this value. The default
%      value is gain = 1.15.
%
%    'lambda' - Weight of the input signal.
%      For each output sample, covered by K overlapping blocks, its final
%      value is obtained as a weighted average of the K denoised values and
%      the original noisy value. In the averaging, each denoised value is
%      attributed a weight of 1, and the noisy value is attributed a weight
%      of lambda. In this way, samples with less confidence (covered by
%      less blocks) are regularized by giving larger weight to the original
%      sample. Default value for lambda is 0.1*(maxval/sigma), where sigma
%      is the standard deviation of the noise. Specifying lambda=0 means
%      that the noisy signal is not averaged into the result.
%
%    'maxatoms' - Maximal number of atoms.
%      This parameter can be used to specify a hard limit on the number of
%      atoms used to sparse-code each block (see parameter 'maxatoms' in
%      OMP2 for more details). Default value is prod(blocksize)/2, i.e.
%      half the number of samples in a block.
%
%
%   Summary of all fields in PARAMS:
%   --------------------------------
%
%   Required:
%     'x'                   signal to denoise
%     'blocksize'           size of block to process
%     'dict'                dictionary to denoise each block
%     'psnr' / 'sigma'      noise power in dB / standard deviation
%
%   Optional (default values in parenthesis):
%     'stepsize'            distance between neighboring blocks (1)
%     'maxval'              maximal intensity value (1)
%     'memusage'            'low, 'normal' or 'high' ('normal')
%     'noisemode'           'psnr' or 'sigma' ('sigma')
%     'gain'                noise gain (1.15)
%     'lambda'              weight of input signal (0.1*maxval/sigma)
%     'maxatoms'            max # of atoms per block (prod(blocksize)/2)
%
%
%  References:
%  [1] M. Elad and M. Aharon, "Image Denoising via Sparse and Redundant
%      representations over Learned Dictionaries", the IEEE Trans. on Image
%      Processing, Vol. 15, no. 12, pp. 3736-3745, December 2006.
%
%  See also OMPDENOISE1, OMPDENOISE2, OMPDENOISE3, KSVDDENOISE, OMP, OMP2.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% parse input arguments %
addpath('ompbox');

x = params.x;
blocksize = params.blocksize;
D = params.dict;

p = ndims(x);
if (p==2 && any(size(x)==1) && length(blocksize)==1)
  p = 1;
  x = x(:);
end


% blocksize %
if (numel(blocksize)==1)
  blocksize = ones(1,p)*blocksize;
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
    stepsize = ones(1,p)*stepsize;
  end
else
  stepsize = ones(1,p);
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


% adaptation for 1-D signals %
if (p==1)
  blocksize = [blocksize 1];
  stepsize = [stepsize 1];
end


% denoise the signal %

if (memusage >= MEM_NORMAL)
  G = D'*D;
end

y = zeros(size(x));

% ids{} contains the indices of the current block
ids = cell(p,1);
for j = 1:p
  ids{j} = 1:blocksize(j);
end

% lastids contains the indices of the last block in each dimension
lastids = stepsize .* floor((size(x)-blocksize)./stepsize) + 1;

blocknum = prod(floor((size(x)-blocksize)./stepsize) + 1);
tid = timerinit('ompdenoise', blocknum); lastmsgcheck = 0;

nz = 0;  % count non-zeros in block representations
for i = 1 : blocknum
 
  block = x(ids{:});
  block = block(:);
  dc = mean(block);
  block = block-dc;
  if (memusage == MEM_LOW)
    gamma = omp2(D,block,[],epsilon,'maxatoms',maxatoms);
  else
    gamma = omp2(D'*block,sum(block.*block),G,epsilon,'maxatoms',maxatoms);
  end
  nz = nz + nnz(gamma);
  y(ids{:}) = y(ids{:}) + reshape(D*gamma,blocksize) + dc;

  % increment block ids
  if (i<blocknum)
    j = 1;
    while (ids{j}(1) == lastids(j))
      ids{j} = 1:blocksize(j);
      j = j+1;
    end
    ids{j} = ids{j}+stepsize(j);
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
y = reshape(y,size(params.x));
