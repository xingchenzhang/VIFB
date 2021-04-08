function [y,D,nz] = ksvddenoise(params,msgdelta)
%KSVDDENOISE K-SVD denoising.
%  [Y,D] = KSVDDENOISE(PARAMS) denoises the specified (possibly
%  multi-dimensional) signal using K-SVD denoising. Y is the denoised
%  signal and D is the trained dictionary produced by K-SVD.
%
%  [Y,D] = KSVDDENOISE(PARAMS,MSGDELTA) specifies the frequency of message
%  printing during the process. MSGDELTA should be a positive number
%  representing the interval in seconds between messages. A zero or
%  negative value cancels all messages. Default is MSGDELTA=5.
%
%  [Y,D,NZ] = KSVDDENOISE(...) also returns the average number of non-zero
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
%    'dictsize' - Size of dictionary to train.
%      Specifies the number of dictionary atoms to train by K-SVD.
%
%    'psnr' / 'sigma' - Noise power.
%      Specifies the noise power in dB (psnr) or the noise standard
%      deviation (sigma), used to determine the target error for
%      sparse-coding each block. If both fields are present, sigma is used
%      unless the field 'noisemode' is specified (below). When specifying
%      the noise power in psnr, make sure to set the 'maxval' parameter
%      as well (below) if the signal values are not within [0,1].
%
%    'trainnum' - Number of training blocks.
%      Specifies the number of training blocks to extract from the noisy
%      signal for K-SVD training.
%
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'initdict' - Initial dictionary.
%      Specifies the initial dictionary for the K-SVD training. Should be
%      either a matrix of size NxL where N=(N1*N2*...*Np), the string
%      'odct' to specify the overcomplete DCT dictionary, or the string
%      'data' to initialize using random signal blocks. When a matrix is
%      specified for 'initdict', L must be >= dictsize, and in this case
%      the dictionary is initialized using the first dictsize columns from
%      initdict. By default, initdict='odct'.
%
%    'stepsize' -  Interval between neighboring blocks.
%      Specifies the interval (in pixels/voxels) between neighboring blocks
%      to denoise in the OMP denoising step. By default, all overlapping
%      blocks are denoised and averaged. This can be changed by specifying
%      an alternate stepsize, as an array of the form [S1 S2 ... Sp] (where
%      p is the number of dimensions of x). This sets the distance between
%      neighboring blocks to be Si in the i-th dimension. Stepsize can also
%      be a scalar S, corresponding to the step size [S S ... S]. Each Si
%      must be >= 1, and, to ensure coverage of the entire noisy signal,
%      size(x,i)-Ni should be a multiple of Si for all i. The default
%      stepsize is 1.
%
%    'iternum' - Number of K-SVD iterations.
%      Specifies the number of K-SVD training iterations to perform. If not
%      specified, the default is 10.
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
%      When 'memusage' is specified, both KSVD and OMPDENOISE are invoked
%      using their corresponding memusage setting. Additionallly, when set
%      to 'low', the use of the functions OMPDENOISE2 and OMPDENOISE3 to
%      accelerate denoising of 2-D and 3-D signals is disabled, reverting
%      to the function OMPDENOISE. Note that the 'low' setting will
%      significantly increase runtime.
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
%      the value derived from the psnr/sigma fields. When gain is different
%      from 1, the target error is multiplied by this value. The default
%      value is gain = 1.15.
%
%    'lambda' - Weight of the input signal.
%      Specifies the relative weight attributed to the noisy input signal
%      in determining the output. The default value is 0.1*(maxval/sigma),
%      where sigma is the standard deviation of the noise. See function
%      OMPDENOISE for more information.
%
%    'maxatoms' - Maximal number of atoms.
%      This parameter can be used to specify a hard limit on the number of
%      atoms used to sparse-code each block. Default value is
%      prod(blocksize)/2, i.e. half the number of samples in a block. See
%      function OMP2 for more information.
%
%    'exact' - Exact K-SVD update.
%      Specifies whether the exact or approximate dictionary update should
%      be used in the K-SVD training. By default, the approximate
%      computation is used. However, specifying a nonzero value for 'exact'
%      causes the exact computation to be used instead. See function KSVD
%      for more information.
%
%
%   Summary of all fields in PARAMS:
%   --------------------------------
%
%   Required:
%     'x'                      signal to denoise
%     'blocksize'              size of block to process
%     'dictsize'               size of dictionary to train
%     'psnr' / 'sigma'         noise power in dB / standard deviation
%     'trainnum'               number of training signals
%
%   Optional (default values in parenthesis):
%     'initdict'               initial dictionary ('odct')
%     'stepsize'               distance between neighboring blocks (1)
%     'iternum'                number of training iterations (10)
%     'maxval'                 maximal intensity value (1)
%     'memusage'               'low, 'normal' or 'high' ('normal')
%     'noisemode'              'psnr' or 'sigma' ('sigma')
%     'gain'                   noise gain (1.15)
%     'lambda'                 weight of input signal (0.1*maxval/sigma)
%     'maxatoms'               max # of atoms per block (prod(blocksize)/2)
%     'exact'                  exact update instead of approximate (0)
%
%
%  References:
%  [1] M. Elad and M. Aharon, "Image Denoising via Sparse and Redundant
%      representations over Learned Dictionaries", the IEEE Trans. on Image
%      Processing, Vol. 15, no. 12, pp. 3736-3745, December 2006.
%
%  See also KSVD, OMPDENOISE, OMP2.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


%%%%% parse input parameters %%%%%
addpath('ompbox');

x = params.x;
blocksize = params.blocksize;
trainnum = params.trainnum;
dictsize = params.dictsize;

p = ndims(x);
if (p==2 && any(size(x)==1) && length(blocksize)==1)
  p = 1;
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
  params.maxval = maxval;
end


% gain %
if (isfield(params,'gain'))
  gain = params.gain;
else
  gain = 1.15;
  params.gain = gain;
end


% msgdelta %
if (nargin<2)
  msgdelta = 5;
end

verbose = 't';
if (msgdelta <= 0)
  verbose='';
  msgdelta = -1;
end


% initial dictionary %

if (~isfield(params,'initdict'))
  params.initdict = 'odct';
end

if (isfield(params,'initdict') && ischar(params.initdict))
  if (strcmpi(params.initdict,'odct'))
    params.initdict = odctndict(blocksize,dictsize,p);
  elseif (strcmpi(params.initdict,'data'))
    params = rmfield(params,'initdict');    % causes initialization using random examples
  else
    error('Invalid initial dictionary specified.');
  end
end

if (isfield(params,'initdict'))
  params.initdict = params.initdict(:,1:dictsize);
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

params.Edata = sqrt(prod(blocksize)) * sigma * gain;   % target error for omp
params.codemode = 'error';

params.sigma = sigma;
params.noisemode = 'sigma';


% make sure test data is not present in params
if (isfield(params,'testdata'))
  params = rmfield(params,'testdata');
end


%%%% create training data %%%

ids = cell(p,1);
if (p==1)
  ids{1} = reggrid(length(x)-blocksize+1, trainnum, 'eqdist');
else
  [ids{:}] = reggrid(size(x)-blocksize+1, trainnum, 'eqdist');
end
params.data = sampgrid(x,blocksize,ids{:});
params.data = remove_dc(params.data,'columns');


%%%%% KSVD training %%%%%

if (msgdelta>0)
  disp('KSVD training...');
end
D = ksvd(params,verbose,msgdelta);


%%%%%  denoise the signal  %%%%%

if (~isfield(params,'lambda'))
  params.lambda = maxval/(10*sigma);
end

params.dict = D;

if (msgdelta>0)
  disp('OMP denoising...');
end

% call the appropriate ompdenoise function
if (p>3  || (isfield(params,'memusage') && strcmpi(params.memusage,'low')))
  [y,nz] = ompdenoise(params,msgdelta);
elseif (p==1)
  [y,nz] = ompdenoise1(params,msgdelta);
elseif (p==2)
  [y,nz] = ompdenoise2(params,msgdelta);
else
  [y,nz] = ompdenoise3(params,msgdelta);
end

return;

