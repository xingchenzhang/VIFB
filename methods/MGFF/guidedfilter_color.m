  function A = guidedfilter_color(I, p, r, eps)
%   GUIDEDFILTER_COLOR   O(1) time implementation of guided filter using a color image as the guidance.
%
%   - guidance image: I (should be a color (RGB) image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps

% r = 16;
% eps = 0.1^2;
% I=double(imread('igloo1.jpg'));
% p=I;
q = zeros(size(I));

q(:, :, 1) = guidedfilter(I(:, :, 1), p(:, :, 1), r, eps);
q(:, :, 2) = guidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
q(:, :, 3) = guidedfilter(I(:, :, 3), p(:, :, 3), r, eps);% Eqn. (16) in the paper;
 
A=cat(3, q(:,:,1),q(:,:,2),q(:,:,3));
% figure, imshow(uint8(A));
 end