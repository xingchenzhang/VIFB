disp(' ');
disp('fdct_wrapping_demo_recon.m -- Partial Curvelet reconstruction.');
disp(' ');
disp('We apply the curvelet transform to an image, select a percentage');
disp('of the largest coefficients (in modulus), and set the others');
disp('to zero. We then take the inverse curvelet transform to obtain');
disp('a partial reconstruction of the original image.');
disp(' ');

% fdct_wrapping_demo_recon.m -- Partial curvelet reconstruction.

% Set the percentage of coefficients used in the partial reconstruction 
pctg = 0.1;

% Load image
X = imread('Lena.jpg'); %load Lena; X = Lena; clear Lena;

% Forward curvelet transform
disp('Take curvelet transform: fdct_wrapping');
tic; C = fdct_wrapping(double(X),0); toc;

% Get threshold value
cfs =[];
for s=1:length(C)
  for w=1:length(C{s})
    cfs = [cfs; abs(C{s}{w}(:))];
  end
end
cfs = sort(cfs); cfs = cfs(end:-1:1);
nb = round(pctg*length(cfs));
cutoff = cfs(nb);

% Set small coefficients to zero
for s=1:length(C)
  for w=1:length(C{s})
    C{s}{w} = C{s}{w} .* (abs(C{s}{w})>cutoff);
  end
end

disp('Take inverse curvelet transform: ifdct_wrapping');
tic; Y = ifdct_wrapping(C,0); toc;

subplot(1,2,1); colormap gray; imagesc(real(X)); axis('image'); title('original image');
subplot(1,2,2); colormap gray; imagesc(real(Y)); axis('image'); title('partial reconstruction');
