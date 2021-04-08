function f = ldfilter(fname)
% LDFILTER	Generate filter for the ladder structure network
%
%	f = ldfilter(fname)
%
% Input:
%	fname:  Available 'fname' are:
%           'pkvaN': length N filter from Phoong, Kim, Vaidyanathan and Ansari

switch fname
    case {'pkva12', 'pkva'}
        v = [0.6300   -0.1930    0.0972   -0.0526    0.0272   -0.0144];
        
    case {'pkva8'}
        v = [0.6302   -0.1924    0.0930   -0.0403];
        
    case {'pkva6'}
        v = [0.6261   -0.1794    0.0688];
        
    otherwise
        error('Unrecognized ladder structure filter name');
end

% Symmetric impulse response
f = [v(end:-1:1), v];