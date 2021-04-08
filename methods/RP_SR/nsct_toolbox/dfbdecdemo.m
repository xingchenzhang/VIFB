
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image decomposition by nonsubsampled contourlet transform (NSSC).
% This is the iterated filter bank that computes the nonsubsampled
% contourlet transform.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameteters:
nlevels = 4 ;        % Decomposition level
%pfilter = 'pkva' ;              % Pyramidal filter
dfilter = 'dmaxflat7'; %'cd' ;              % Directional filter

% Nonsubsampled Contourlet decomposition
coeffs = nssdfbdec( double(im), dfilter, nlevels );
disp( nlevels); disp(dfilter);


% Display the coefficients
%disp('Displaying the contourlet coefficients...') ;
%shownssc( coeffs ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonsubsampled Contourlet transform (NSSC) reconstruction.
% This is the inverse of nsscdec, i.e.
% imrec = nsscrec(coeffs, dfilter, pfilter);
% would reconstruct imrec = im
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reconstruct image
imrec = nssdfbrec( coeffs, dfilter ) ;


% Show the reconstruction image and the original image
if 0
figure;
subplot(1,2,1), imagesc( im, [0, 255] ); 
title('Original image' ) ;
colormap(gray);
axis image off;
subplot(1,2,2), imagesc( imrec, [0, 255] );
title('Reconstructed image' ) ;
colormap(gray);
axis image off;
end

mse = sum( sum( (imrec - double(im)).^2 ) );
mse = mse / prod(size(im));

disp( sprintf('The mean square error is: %f', mse ) );
disp(' ');