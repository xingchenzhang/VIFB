function y = nsctdec(x, levels, dfilt, pfilt )
% NSSCDEC   Nonsubsampled Contourlet Transform Decomposition
%
%	y = nsctdec(x, levels, [dfilt, pfilt] )
%
% Input:
%   x:      
%       a matrix, input image
%   levels:  
%       vector of numbers of directional filter bank decomposition levels 
%       at each pyramidal level (from coarse to fine scale).
%       If the number of level is 0, a critically sampled 2-D wavelet 
%       decomposition step is performed.
%       The support for the wavelet decomposition has not been verified!!!!!!
%   dfilt:  
%       a string, filter name for the directional decomposition step.
%       It is optional with default value 'dmaxflat7'. See dfilters.m for all
%       available filters.
%   pfilt:  
%       a string, filter name for the pyramidal decomposition step.
%       It is optional with default value 'maxflat'. See atrousfilters.m for 
%       all available filters. 
%
% Output:
%   y:  a cell vector of length length(nlevs) + 1, where except y{1} is 
%       the lowpass subband, each cell corresponds to one pyramidal
%       level and is a cell vector that contains bandpass directional
%       subbands from the DFB at that level.
%
% Index convention:
%   Suppose that nlevs = [l_J,...,l_2, l_1], and l_j >= 2.
%   Then for j = 1,...,J and k = 1,...,2^l_j 
%       y{J+2-j}{k}(n_1, n_2)
%   is a contourlet coefficient at scale 2^j, direction k, and position
%       (n_1 * 2^(j+l_j-2), n_2 * 2^j) for k <= 2^(l_j-1), 
%       (n_1 * 2^j, n_2 * 2^(j+l_j-2)) for k > 2^(l_j-1).
%   As k increases from 1 to 2^l_j, direction k rotates clockwise from
%   the angle 135 degree with uniform increment in cotan, from -1 to 1 for
%   k <= 2^(l_j-1), and then uniform decrement in tan, from 1 to -1 for 
%   k > 2^(l_j-1).
%
% See also:	ATROUSFILTERS, DFILTERS, NSCTREC, NSFBDEC, NSDFBDEC.

%  History:
%     02/17/04 Created by Jianping Zhou.
%     08/07/04 Modified by Jianping Zhou. Incorporte the fast implementation 
%              of convolution algorithm by Jason Laska. 
%     08/30/04 Modified by Arthur. L. Cunha, added more pyramid filters.
%     10/17/04 Modified by Arthur. L. Cunha, replaced periodic with symmetric extension
%     01/24/05 Modified by Jianping Zhou, changed function names and corrected a bug. 
%     10/31/05 Modified by Arthur L Cunha, corrected a bug in the pyramid decomposition (a wrong index)

% Check input
if ~isnumeric( levels )
    error('The decomposition levels shall be integers');
end
if isnumeric( levels )
    if round( levels ) ~= levels
        error('The decomposition levels shall be integers');
    end
end
        
if ~exist('dfilt', 'var')
    dfilt = 'dmaxflat7' ;
end;

if ~exist('pfilt', 'var')
    pfilt = 'maxflat' ; 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get fan filters, parallelogram filters, and pyramid filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the directional filters for the critically sampled DFB.
filters = cell(4) ;
[h1, h2] = dfilters(dfilt, 'd');
% A scale is required for the nonsubsampled case.
h1 = h1./sqrt(2) ;
h2 = h2./sqrt(2) ;

% Generate the first-level fan filters by modulations.
filters{1} = modulate2(h1, 'c');
filters{2} = modulate2(h2, 'c'); 

% Obtain the parallelogram filters from the diamond filters
[filters{3}, filters{4}] = parafilters( h1, h2 ) ;

[h1, h2, g1, g2] = atrousfilters(pfilt); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsubsampled Contourlet transform with tree structure filter banks
% Nonsubsampled pyramids make multiresolution decomposition.
% Nonsubsampled directional filter banks make directional decomposition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numer of levels:
clevels = length( levels ) ;
nIndex = clevels + 1 ;
% Initialize the output
y = cell(1, nIndex) ;

% Nonsubsampled pyramid decomposition
for i= 1 : clevels    
    
    % Nonsubsampled Contourlet transform
    % Nonsubsampled pyramids decomposition
    [xlo, xhi] = nsfbdec(x, h1, h2, i-1) ;
        
    if levels(nIndex-1) > 0        
        % Nonsubsampled DFB decomposition on the bandpass image
        xhi_dir = nsdfbdec(xhi, filters, levels(nIndex-1));        
        y{nIndex}=xhi_dir ;
    else
        % Copy the result directly:
        y{nIndex}=xhi ;
    end
    
    % Update the index for the Nonsubsampled Pyramdis
    nIndex = nIndex - 1 ;
    
    % Prepare for next iteration
    x = xlo ;
    
end

% The lowpass output
y{1}=x;
