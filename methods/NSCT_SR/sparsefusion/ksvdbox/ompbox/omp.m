function gamma = omp(varargin)
%OMP Sparsity-constrained Orthogonal Matching Pursuit.
%  GAMMA = OMP(D,X,G,T) solves the optimization problem
%
%       min  |X - D*GAMMA|_2     s.t.  |GAMMA|_0 <= T
%      gamma
%
%  for each of the signals in X, using Batch Orthogonal Matching Pursuit.
%  Here, D is a dictionary with normalized columns, X is a matrix
%  containing column signals, T is the # if non-zeros in each signal
%  representation, and G is the Gramm matrix D'*D. The output GAMMA is a
%  matrix containing the sparse representations as its columns. 
%
%  GAMMA = OMP(D,X,[],T) performs the same operation, but without the
%  matrix G, using OMP-Cholesky. This call produces the same output as
%  Batch-OMP, but is significantly slower. Using this syntax is only
%  recommended when available memory is too small to store G.
%
%  GAMMA = OMP(DtX,G,T) is the fastest implementation of OMP, but also
%  requires the most memory. Here, DtX stores the projections D'*X. In this
%  case Batch-OMP is used, but without having to compute D'*X in advance,
%  which slightly improves runtime. Note that in general, the call
%
%    GAMMA = OMP(D'*X,G,T);
%
%  will be faster than the call
%
%    GAMMA = OMP(D,X,G,T);
%
%  due to optimized matrix multiplications in Matlab. However, when the
%  entire matrix D'*X cannot be stored in memory, one of the other two
%  versions can be used. Both compute D'*X for just one signal at a time,
%  and thus require much less memory.
%
%  GAMMA = OMP(...,PARAM1,VAL1,PARAM2,VAL2,...) specifies additional
%  parameters for OMP. Available parameters are:
%
%    'gammamode' - Specifies the representation mode for GAMMA. Can be
%                  either 'full' or 'sparse', corresponding to a full or
%                  sparse matrix, respectively. By default, GAMMA is
%                  returned as a sparse matrix.
%    'messages'  - Specifies whether progress messages should be displayed.
%                  When positive, this is the number of seconds between
%                  status prints. When negative, indicates that no messages
%                  should be displayed (this is the default).
%    'profile'   - Can be either 'on' or 'off'. When 'on', profiling
%                  information is displayed at the end of the funciton
%                  execution.
%
%
%  Summary of OMP versions:
%
%    version      |   speed     |   memory
%  --------------------------------------------------
%   OMP(DtX,G,T)  |  very fast  |  very large
%   OMP(D,X,G,T)  |  fast       |  moderate
%   OMP(D,X,[],T) |  very slow  |  small
%  --------------------------------------------------
%
%
%  References:
%  [1] M. Elad, R. Rubinstein, and M. Zibulevsky, "Efficient Implementation
%      of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit",
%      Technical Report - CS, Technion, April 2008.
%
%  See also OMP2.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009


% default options

sparse_gamma = 1;
msgdelta = -1;
profile = 0;


% determine number of parameters

paramnum = 1;
while (paramnum<=nargin && ~ischar(varargin{paramnum}))
  paramnum = paramnum+1;
end
paramnum = paramnum-1;


% parse options

for i = paramnum+1:2:length(varargin)
  paramname = varargin{i};
  paramval = varargin{i+1};

  switch lower(paramname)

    case 'gammamode'
      if (strcmpi(paramval,'sparse'))
        sparse_gamma = 1;
      elseif (strcmpi(paramval,'full'))
        sparse_gamma = 0;
      else
        error('Invalid GAMMA mode');
      end

    case 'messages'
      msgdelta = paramval;
      
    case 'profile'
      if (strcmpi(paramval,'on'))
        profile = 1;
      elseif (strcmpi(paramval,'off'))
        profile = 0;
      else
        error('Invalid profile mode');
      end

    otherwise
      error(['Unknown option: ' paramname]);
  end
  
end


% determine call type

if (paramnum==3)
  DtX = varargin{1};
  G = varargin{2};
  T = varargin{3};
  D = [];
  X = [];
elseif (paramnum==4)
  D = varargin{1};
  X = varargin{2};
  G = varargin{3};
  T = varargin{4};
  DtX = [];
else
  error('Invalid number of parameters');
end


% verify dictionary normalization

if (isempty(G))
  atomnorms = sum(D.*D);
else
  atomnorms = diag(G);
end
if (any(abs(atomnorms-1) > 1e-2))
  error('Dictionary columns must be normalized to unit length');
end


% omp

gamma = ompmex(D,X,DtX,G,T,sparse_gamma,msgdelta,profile);
