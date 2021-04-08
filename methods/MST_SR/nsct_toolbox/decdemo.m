function coeffs = decdemo( im, option )
% DECDEMO  demonstrates nonsubsampled Contourlet decomposition and reconstruction. 
%
%   DECDEMO shows how to use the nonsubsampled contourlet toolbox to decompose
%   and reconstruct an image.  It provides a sample script that uses 
%   basic functions such as nsctdec, nsctrec, and shownsct.
%
%   It can be modified for applications such as image analysis, 
%   image retrieval and image processing.
%
% Input:
%	image:  a double or integer matrix for the input image.
%           The default is the zoneplate image.
%   option: option for the demos. The default value is 'auto'
%       'auto' ------  automtatical demo, no input
%       'user' ------  semi-automatic demo, simple interactive inputs
%       'expert' ----  mannual, complete interactive inputs. 
%                      (Not implmented in this version)
%
% Output:
%	coeffs: a cell vector for the contourlet decomposition coefficients.
%   
% See also:     NSCTDEC, NSCTREC, SHOWNSCT.

% History:
%   08/08/2004  Created by Jianping Zhou.

disp('Welcome to the nonsubsampled Contourlet decomposition demo! :)');
disp('Type help decdemo for help' ) ;
disp('You can also view decdemo.m for details.') ;
disp(' ');

% Input image
if ~exist('im', 'var')
    % Zoneplate image: good for illustrating multiscale and directional
    % decomposition
    im = imread ('zoneplate.png') ;
elseif isstr(im)
    im = imread ( im ) ;
else
    error('You shall input valid image name!');
end

% Show the input image
disp( 'Displaying the input image...');
clf;
imagesc(im, [0, 255]);
title( 'Input image' ) ;
axis image off;
colormap(gray);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image decomposition by nonsubsampled contourlet transform (NSCT).
% This is the iterated filter bank that computes the nonsubsampled
% contourlet transform.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameteters:
nlevels = [0, 1, 3] ;        % Decomposition level
pfilter = 'maxflat' ;              % Pyramidal filter
dfilter = 'dmaxflat7' ;              % Directional filter

% Nonsubsampled Contourlet decomposition
coeffs = nsctdec( double(im), nlevels, dfilter, pfilter );

% Display the coefficients
disp('Displaying the contourlet coefficients...') ;
shownsct( coeffs ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsubsampled Contourlet transform (NSCT) reconstruction.
% This is the inverse of nsctdec, i.e.
% imrec = nsctrec(coeffs, dfilter, pfilter);
% would reconstruct imrec = im
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reconstruct image
imrec = nsctrec( coeffs, dfilter, pfilter ) ;

disp('Displaying the reconstructed image...') ;
disp('It should be a perfect reconstruction' ) ;
disp(' ') ;

% Show the reconstruction image and the original image
figure;
subplot(1,2,1), imagesc( im, [0, 255] ); 
title('Original image' ) ;
colormap(gray);
axis image off;
subplot(1,2,2), imagesc( imrec, [0, 255] );
title('Reconstructed image' ) ;
colormap(gray);
axis image off;

mse = sum( sum( (imrec - double(im)).^2 ) );
mse = mse / prod(size(im));

disp( sprintf('The mean square error is: %f', mse ) );
disp(' ');