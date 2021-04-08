function y = nssfbrec( x1, x2, f1, f2, mup )
% NSSFBREC   Two-channel nonsubsampled filter bank reconstruction with periodic extension.
%   NSSFBREC reconstructs the image Y by a two-channel nonsubsampled filter bank.
%   The two inputs X1 and X2 will be convolved with two corresponding upsampled 
%   filters F1 and F2 by the sampling matrix MUP. The periodic extension has been
%   considered. There is no subsampling in the signals and hence the operation is
%   shift-invariant.
%  
%       nssfbrec( x1, x2, f1, f2, [mup] )
%
% INPUT:
%   x1:
%       a matrix, the input for the first branch.
%   x2:
%       a matrix, the input for the second branch.
%	f1:	
%		a matrix, the filter for the first branch.
%	f2:	
%		a matrix, the filter for the second branch.
%   mup:
%       an integer or an integer matrix, upsampling matrix. 
%       If it is an integer, the upsampling is separable with same rate in two
%       dimensions. If it does not exist, no upsampling.
%
% OUTPUT:
%	y:
%       a matrix, reconstructed image. 
%
% NOTE:
%   1. The size of X1 and X2 should be equal. 
%   2. There are two mex files (zconv2.c, zconv2S.c) that might need to 
%      be recompiled. This can be done by typing from the Matlab command window
%      >> mex zconv2.c  
%      The name of the generated mex-file is also zconv2, but the extension
%      depends on your operating system. For example, *.dll for Windows, and 
%      *.mexmac for Macintosh. 
%
% See also:     EFILTER2, ZCONV2, ZCONV2S.

%
% History: 
%   08/07/2004  Created by Jianping Zhou.
%               ZCONV2 was created by Jason Laska in July 2004.
%   08/11/2004  Modified by Jianping Zhou to add ZCONV2S. 
%

% Check input
if size(x1) ~= size(x2)
    error('The size of inputs for two branches shall be the same!');
end
if ~exist('mup', 'var')
    % Convolve with the filters with periodic extension.
    y1 = efilter2( x1, f1 );
    y2 = efilter2( x2, f2 );
    y = y1 + y2 ;
    return ;
end
% Use the built-in convolution function when there is no upsampling
if mup == 1 | mup == eye(2)
    % Convolve with the filters with periodic extension.
    y1 = efilter2( x1, f1 );
    y2 = efilter2( x2, f2 );
    y = y1 + y2 ;
    return ;
end

if size(mup) == [2, 2]
% Nonseparable sampling matrix case.

    % Convolve the input with upsampled filters with periodic extension.
    y1 = zconv2( x1, f1, mup );
    y2 = zconv2( x2, f2, mup );
    y = y1 + y2 ;
    
elseif size(mup) == [1, 1]
% Separable sampling matrix case.

    % Convolve the input with upsampled filters with periodic extension.
    mup = mup * eye(2) ;
    y1 = zconv2S ( x1, f1, mup );
    y2 = zconv2S ( x2, f2, mup ); 
    y = y1 + y2 ;
    
else
    error('The upsampling parameter should be an integer or two-dimensional integer matrix!');
end