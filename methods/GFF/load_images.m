% This procedure loads a sequence of images
%
% Arguments:
%   'path', refers to a directory which contains a sequence of images
%   'reduce' is an optional parameter that controls downsampling, e.g., reduce = .5
%   downsamples all images by a factor of 2.
%
% tom.mertens@gmail.com, August 2007
%

function I = load_images(path, reduce)

if ~exist('reduce')
    reduce = 1;
end

if (reduce > 1 || reduce <= 0)
    error('reduce must fulfill: 0 < reduce <= 1');
end

% find all JPEG or PPM files in directory
files = dir([path '/*.tif']);
N = length(files);
if (N == 0)
    files = dir([path '/*.jpg']);
    N = length(files);
    if (N == 0)
    files = dir([path '/*.gif']);
    N = length(files);
    if (N == 0)
    files = dir([path '/*.bmp']);
    N = length(files);
    if (N == 0)
    files = dir([path '/*.png']);
    N = length(files);
    if (N == 0)
    error('no files found');
    end
          end
         end
    end
end

% allocate memory
sz = size(imread([path '/' files(1).name]));
r = floor(sz(1)*reduce);
c = floor(sz(2)*reduce);
if length(sz)==3
    I=zeros(r,c,3,N);
else
I = zeros(r,c,N);
end
% read all files
for i = 1:N
    
    % load image
    filename = [path '/' files(i).name];
    im = imread(filename);
    if (size(im,1) ~= sz(1) || size(im,2) ~= sz(2))
        error('images must all have the same size');
    end
    
    % optional downsampling step
    if (reduce < 1)
    im = imresize(im,[r c],'bicubic');
    end
    if size(im,3)==3
    I(:,:,:,i) =im;
    else
    I(:,:,i) = im;
    end
    I=uint8(I);
end
