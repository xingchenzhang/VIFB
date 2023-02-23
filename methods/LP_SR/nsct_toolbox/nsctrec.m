function x = nsctrec(y, dfilt, pfilt)
% NSCTREC   Nonsubsampled Contourlet Reconstruction
%
%	x = nsscrec(y, [dfilt, pfilt] )
%
% INPUT:
%   y:  a cell vector of length length(nlevs) + 1, where except y{1} is 
%       the lowpass subband, each cell corresponds to one pyramidal
%       level and is a cell vector that contains bandpass directional
%       subbands from the DFB at that level.
%   dfilt:  
%       a string, filter name for the directional decomposition step.
%       It is optional with default value 'dmaxflat7'. See dfilters.m for all
%       available filters.
%   pfilt:  
%       a string, filter name for the pyramidal decomposition step.
%       It is optional with default value 'maxflat'. See atrousfilters.m for 
%       all available filters. 
%
% OUTPUT:
%   x:      
%       a matrix, reconstructed image
%
% See also: ATROUSFILTERS, DFILTERS, NSCTDEC, NSFBREC, NSDFBREC

%
% HISTORY:
%     02/17/04 Created by Jianping Zhou.
%     08/07/04 Modified by Jianping Zhou. Incorporte the fast implementation 
%              of convolution algorithm by Jason Laska.
%     08/30/04 Modified by Arthur. L. Cunha, added more pyramid filters.
%     10/17/04 Modified by Arthur. L. Cunha, replaced periodic with symmetric extension
%     01/24/05 Modified by Jianping Zhou, changed function names and corrected a bug. 
%     10/31/05 Modified by Arthur L Cunha, corrected a bug in the pyramid decomposition (a wrong index)


% Check input
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
[h1, h2] = dfilters(dfilt, 'r');
% A scale is required for the nonsubsampled case.
h1 = h1./sqrt(2) ;
h2 = h2./sqrt(2) ;

% Generate the first-level fan filters by modulations.
filters{1} = modulate2(h1, 'c');
filters{2} = modulate2(h2, 'c'); 

% Obtain the parallelogram filters from the diamond filters
[filters{3}, filters{4}] = parafilters( h1, h2 ) ;

% Currently only one filter by Arthur Cunha
% It has been normalized.
[h1, h2, g1, g2] = atrousfilters(pfilt) ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsubsampled Contourlet transform with tree structure filter banks
% Nonsubsampled pyramids make multiresolution decomposition.
% Nonsubsampled directional filter banks make directional decomposition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(y) - 1;
xlo = y{1};
nIndex = n-1;
for i=1:n
    
    % Process the detail subbands
    if iscell( y{i+1} )
        % Nonsubsampled DFB reconstruction
        xhi = nsdfbrec( y{i+1}, filters ); 
    else
        % No DFB decomposition, copy directly
        xhi = y{i+1};
    end
                 
    % Nonsubsampled Pyramid reconstruction
    x = nsfbrec(xlo, xhi, g1, g2, nIndex);            
        
    % Prepare for the next level
    xlo = x ;
    nIndex = nIndex -1;     
end