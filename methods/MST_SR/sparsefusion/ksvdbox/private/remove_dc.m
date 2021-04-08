function [y,dc] = remove_dc(x,columns)
%REMOVE_DC Remove DC channel from signals.
%   [Y,DC] = REMOVE_DC(X) removes the DC channel (i.e. the mean) from the
%   specified (possibly multi-dimensional) signal X. Y is the DC-free
%   signal and is the same size as X. DC is a scalar containing the mean of
%   the signal X.
%
%   [Y,DC] = REMOVE_DC(X,'columns') where X is a 2D matrix, treats the
%   columns of X as a set of 1D signals, removing the DC channel from each
%   one individually. Y is the same size as X and contains the DC-free
%   signals. DC is a row vector of length size(X,2) containing the means of
%   the signals in X.
%
%   See also ADD_DC.


if (nargin==2 && strcmpi(columns,'columns')), columns = 1;
else columns = 0;
end

if (columns)
  dc = mean(x);
  y = addtocols(x,-dc);
else
  if (ndims(x)==2)  % temporary, will remove in future
    warning('Treating 2D matrix X as a single signal and not each column individually');
  end
  dc = mean(x(:));
  y = x-dc;
end

