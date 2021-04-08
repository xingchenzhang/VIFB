disp(' ');
disp('fdct_wrapping_demo_disp.m -- Displays the curvelet coefficients of an image');
disp(' ');

disp('1. The low frequency (coarse scale) coefficients are stored at');
disp('   the center of the display.') 
disp(['2. The Cartesian concentric coronae show the coefficients at different']); 
disp('   scales; the outer coronae correspond to higher frequencies.');
      
disp(['3. There are four strips associated to each corona, corresponding to']);
disp('   the four cardinal points; these are further subdivided in angular panels.');     
disp(['4. Each panel represent coefficients at a specified scale and along']);
disp('   the orientation suggested by the position of the panel.');
disp(' ');

% fdct_wrapping_demo_disp.m -- Displays the curvelet coefficients of an image

%generate image
n = 512;
X = zeros(n,n);
ix = ((-n/2):(n/2-1))' * ones(1,n);
iy = ones(n,1) * ((-n/2):(n/2-1));
ix = ix./n; iy = iy./n;
X = ((ix.^2 + iy.^2) <= .33^2) .* exp( - 2.*ix.^2 - 2.*iy.^2 ); 
alpha = pi/2 + pi/8;
X = (cos(alpha).*ix + sin(alpha).*iy >= 0).*exp(-16.*ix.^2- 16.*iy.^2 );

%forward curvelet transform
disp('Take curvelet transform: fdct_wrapping');
tic; C = fdct_wrapping(X,0); toc;

%generate curvelet image (a complex array)
img = fdct_wrapping_dispcoef(C);

%display original image and the curvelet coefficient image
subplot(1,2,1); colormap gray; imagesc(X); axis('image'); title('original image');
subplot(1,2,2); colormap gray; imagesc(abs(img)); axis('image'); title('log of curvelet coefficients');

