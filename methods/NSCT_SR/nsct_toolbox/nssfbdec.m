function [y1, y2] = nssfbdec( x, f1, f2, mup )
% NSSFBDEC   Two-channel nonsubsampled filter bank decomposition with periodic extension.
%   NSSFBDEC convolve the image X by a two-channel nonsubsampled filter bank
%   with two analysis filters F1 and F2 upsampled by the sampling matrix MUP. 
%   There is no subsampling and hence the operation is shift-invariant.
%  
%       nssfbdec( x, f1, f2, [mup] )
%
% INPUT:
%   x:
%       a matrix, input image.
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
%	y1:
%       a matrix, the output for the first branch. 
%	y2:
%       a matrix, the output for the second branch. 
%
% NOTE:
%   There are two mex files (zconv2.c and zconv2s.c) that might need to 
%   be recompiled. This can be done by typing from the Matlab command window
%   >> mex zconv2.c  
%   The name of the generated mex-file is also zconv2, but the extension depends on
%   your operating system. For example, *.dll for Windows, and *.mexmac for Macintosh. 
%
% See also:     EFILTER2, ZCONV2, ZCONV2S.

% History: 
%   08/06/2004  Created by Jianping Zhou.
%               ZCONV2 was created by Jason Laska in July 2004. 
%   08/11/2004  Modified by Jianping Zhou to add ZCONV2S. 
%               

% Check input
if ~exist('mup', 'var')
    % Convolve with the filters with periodic extension.
    y1 = efilter2( x, f1 );
    y2 = efilter2( x, f2 );
    return ;
end
% Use the built-in convolution function when there is no upsampling
if mup == 1 | mup == eye(2)
    % Convolve with the filters with periodic extension.
    y1 = efilter2( x, f1 );
    y2 = efilter2( x, f2 );
    return ;
end

if size(mup) == [2, 2]
% Nonseparable sampling matrix
    
    % Convolve the input with upsampled filters with periodic extension.
    y1 = zconv2( x, f1, mup );
    y2 = zconv2( x, f2, mup );
    
elseif size(mup) == [1, 1]
% Separable sampling matrix

    % Convolve the input with upsampled filters with periodic extension.
    mup = mup * eye(2) ;
    y1 = zconv2S ( x, f1, mup );
    y2 = zconv2S ( x, f2, mup ); 
        
else
    error('The upsampling parameter should be an integer or two-dimensional integer matrix!');
end



