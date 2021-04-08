function L = latent_matrix(training_path, col, unit)
% Train the extractor of salient features 'L'
%   min_Z,L,E ||Z||_* + ||L||_* +¡¡lambda||E||_1,
%           s.t. X = XZ + LX +E.

fileFolder=fullfile(training_path);
dirOutput=dir(fullfile(fileFolder,'*'));
num = length(dirOutput);
% fileNames={dirOutput.name}'; % all names

B = [];
for i = 3:num
img = imread([training_path,dirOutput(i).name]);
if size(img,3)>1
    img = rgb2gray(img);
end
img = im2double(img);
b = im2col(img, [unit, unit], 'distinct');
if i == 1
    B = b;
else
    B = cat(2, B, b);
end
end
[s1, s2] = size(B);
B_rand = B(:,randperm(s2, col));
lambda = 0.4;
disp('Start-latent lrr');
tic
[Z,L,E] = latent_lrr(B_rand,lambda);
toc
disp('Done-latent lrr');

end




