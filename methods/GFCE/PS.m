% PS.m
% -------------------------------------------------------------------
% This function implements the quality metric "perceptual saliency (PS)"
% for input image img.
% This metric was proposed in "Fusion of infrared and visible images for
% night-vision context enhancement", Applied Optics, 55(23), 2016
% by Zhiqiang Zhou et al.
%
% -------------------------------------------------------------------
% Author: Zhiqiang Zhou

function F = PS(img)

R = 32;
wS = 40;
alpha = 2;

[gx, gy] = GradientSobel(img);
g = sqrt(gx.^2+gy.^2);
img = img(2:end-1, 2:end-1);

img1 = img;
g1 = g;

fun1 = @(block_struct) WinAvgPE(block_struct.data, R);
inf = blockproc(img1, [wS, wS], fun1);

fun2 = @(block_struct) sum((block_struct.data(:)).^alpha);
lambda = blockproc(g1, [wS, wS], fun2);

F = sum(inf(:).*lambda(:))/sum(lambda(:));
end


