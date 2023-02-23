function y = nsdfbdec( x, dfilter, clevels )
% NSDFBDEC   Nonsubsampled directional filter bank decomposition.
%   NSDFBDEC Decompose the image X by a nonsubsampled directional filter bank
%   with a binary-tree structure. It outputs the final branches, totally 2^clevels.
%   There is no subsampling and hence the operation is shift-invariant.
%  
%       nsdfbdec( x, dfilter, [clevels] )
%
% INPUT:
%   x:
%       an array, input image.
%	dfilter:	
%		a string, directional filter name.
%       a cell of matrices, including two directional filters and eight
%       parallelogram filters.
%   clevels:
%       a non-negative integer, number of decomposition levels.
%
% OUTPUT:
%	y:
%       a cell vector, output subbands. 
%
% See also:     DFILTERS, PARAFILTERS, NSSFBDEC.

%
% History: 
%   08/06/2004  Created by Jianping Zhou.

% Input check
if ~exist('clevels', 'var')
    clevels = 0 ;
    y{1} = x;
    return;
end
if (clevels ~= round(clevels)) | (clevels < 0)
    error('Number of decomposition levels must be a non-negative integer');
end
if clevels == 0
    % No decomposition, simply copy input to output
    y{1} = x;    
    return;
end
if ~ischar( dfilter )
    if iscell( dfilter )
        if length( dfilter ) ~= 4
            error('You shall provide a cell of two 2D directional filters and two groups of 2D parallelogram filters!');
        end
    else
        error('You shall provide the name of directional filter or all filters!');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get fan filters, parallelogram filters, and basic sampling matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the diamond filters, if necessary
if ischar( dfilter )
    
    % Get the directional filters for the critically sampled DFB.
    [h1, h2] = dfilters(dfilter, 'd');
    % A scale is required for the nonsubsampled case.
    h1 = h1./sqrt(2) ;
    h2 = h2./sqrt(2) ;
    
    % Generate the first-level fan filters by modulations.
    k1 = modulate2(h1, 'c');
    k2 = modulate2(h2, 'c'); 
    
    % Obtain the parallelogram filters from the diamond filters
    [f1, f2] = parafilters( h1, h2 ) ;

else
    % Copy the fan filters directly.
    k1 = dfilter{1} ;
    k2 = dfilter{2} ;
    
    % Copy the parallelogram filters directly.
    f1 = dfilter{3} ;
    f2 = dfilter{4} ;    
end 

% Quincunx sampling matrices
q1 = [1, -1; 1, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First-level Decompositions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clevels == 1 
    % No upsampling for filters at the first-level.
    [y{1}, y{2}] = nssfbdec( x, k1, k2 ) ;        
    
else %Others    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second-level Decompositions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % No upsampling at filters for the first-level.
    [x1, x2] = nssfbdec( x, k1, k2 ) ;

    % Convolution with upsampled filters
    [y{1}, y{2}] = nssfbdec( x1, k1, k2, q1 ) ;
    [y{3}, y{4}] = nssfbdec( x2, k1, k2, q1 ) ;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third and higher levels Decompositions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    % Third and higher levels decompositions
    for l = 3:clevels
        % Allocate space for the new subband outputs
        y_old = y;    
        y = cell(1, 2^l);
	
        % The first half channels:
        for k = 1:2^(l-2)
            
            % Compute the upsampling matrix by the formula (3.18) of Minh N. Do's 
            % thesis. The upsampling matrix for the channel k in a l-levels DFB is
            % M_k^{(l-1)} (refer to (3.18), pp. 53, Minh N. Do's thesis)
            
            % Compute s_{(l-1)}(k):
            slk = 2*floor( (k-1) /2 ) - 2^(l-3) + 1 ;
            % Compute the sampling matrix:
            mkl = 2*[ 2^(l-3), 0; 0, 1 ]*[1, 0; -slk, 1]; 
            i = mod(k-1, 2) + 1;
            % Decompose by the two-channel filter bank:
            [y{2*k-1}, y{2*k}] = nssfbdec( y_old{k}, f1{i}, f2{i}, mkl );
        end	
	
        % The second half channels:
        for k = 2^(l-2)+1 : 2^(l-1)
            
            % Compute the upsampling matrix by the extension of the formula (3.18) 
            % of Minh N. Do's thesis to the second half channels.
            % thesis. The upsampling matrix for the channel k in a l-levels DFB is
            % M_k^{(l-1)} (refer to notes by Jianping Zhou)
        
            % Compute s_{(l-1)}(k):
            slk = 2 * floor( ( k-2^(l-2)-1 ) / 2 ) - 2^(l-3) + 1 ;
            % Compute the sampling matrix:
            mkl = 2*[ 1, 0; 0, 2^(l-3) ]*[1, -slk; 0, 1]; 
            i = mod(k-1, 2) + 3;
            % Decompose by the two-channel filter bank:
            [y{2*k-1}, y{2*k}] = nssfbdec( y_old{k}, f1{i}, f2{i}, mkl );
        end
    end
end