function displayIm = shownsct( y )
% SHOWNSSC   Show nonsubsampled Contourlet transform coefficients. 
%
%       shownsct(y)
% Input:
%	y:	a cell vector of length n+1, one for each layer of 
%		subband images from NSCT, y{1} is the lowpass image
%
% NOTE:
%       It need further improvement later!!!!

% History:
%   08/08/2003  Created by Jianping Zhou.

% Level of decomposition.
clevels = length( y ) ;

% Show the subband images.
for i=1:clevels
    figure;
    if iscell( y{i} )
        % The number of directional subbands.
        csubband = length( y{i} ) ;
        if csubband > 7
            col = 4 ;
        else
            col = 2 ;
        end
        row = csubband / col ;
        for j = 1:csubband
            subplot( row, col, j ) ;
            imshow( uint8(y{i}{j}) );
            title( sprintf('NSSC coefficients: level %d', i) );
        end
    else
        imshow ( uint8(y{i}) ) ;
        title( sprintf('Nonsubsampled Contourlet coefficients level %d', i) );
    end
end