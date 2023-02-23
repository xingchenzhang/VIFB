function [y1, y2] = parafilters( f1, f2 )
% PARAFILTERS   Generate four groups of parallelogram filters.
%   PARAFILTERS generates four groups of parallelogram filters from a pair of 
%   diamond filters F1 and F2 by modulation and rotation operations.
%  
%       parafilters( f1, f2 )
%
% INPUT:
%	f1:	
%		a matrix, the filter for the first branch.
%	f2:	
%		a matrix, the filter for the second branch.
%
% OUTPUT:
%	y1:
%       a cell vector of matrices, four outputs for the first branch. 
%	y2:
%       a cell vector of matrices, four outputs for the second branch. 
%
% NOTE:
%   Refer to pp. 51, Minh N. Do's thesis.
%
% History: 
%   08/06/2004  Created by Jianping Zhou.
%
% See also:     MODULATE2, RESAMPZ.

% Initialize output
y1 = cell(1, 4);
y2 = cell(1, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain fan filters from the diamond filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modulation operation
y1{1} = modulate2(f1, 'r');
y2{1} = modulate2(f2, 'r');
y1{2} = modulate2(f1, 'c');
y2{2} = modulate2(f2, 'c');

% Transpose operation
y1{3} = y1{1}';
y2{3} = y2{1}';
y1{4} = y1{2}';
y2{4} = y2{2}';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain parallelogram filters from the fan filters
% Use the rotation sampling matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:4
    % Resample the filters by corresponding rotation matrices
    y1{i} = resampz( y1{i}, i ) ;
    y2{i} = resampz( y2{i}, i ) ;
end