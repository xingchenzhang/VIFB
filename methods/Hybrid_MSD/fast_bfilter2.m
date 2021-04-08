%    This function implements a fast 2-D bilateral filtering based
%    on the C/C++ codes provided by Sylvain Paris and Francois Sillion.
%    The method is described in their paper:
%       
%      Paris, S., & Durand, F. (2006). A fast approximation of the 
%      bilateral filter using a signal processing approach.In Computer  
%      Vision¨CECCV 2006 (pp. 568-580).Springer Berlin Heidelberg. 
%      
%    Zhiqiang Zhou, Beijing Institute of Technology
%    April, 2015

function B = fast_bfilter2(A,sigma)

% Verify that the input image exists and is valid.
if ~exist('A','var') || isempty(A)
   error('Input image A is undefined or invalid.');
end
if ~isfloat(A) || ~sum([1,3] == size(A,3)) || ...
      min(A(:)) < 0 || max(A(:)) > 1
   error(['Input image A must be a double precision ',...
          'matrix of size NxMx1 or NxMx3 on the closed ',...
          'interval [0,1].']);      
end

% Verify bilateral filter standard deviations.
if ~exist('sigma','var') || isempty(sigma) || ...
      numel(sigma) ~= 2 || sigma(1) <= 0 || sigma(2) <= 0
   sigma = [3 0.1];
end

%Vertify the existence of mexw32 or mexw64 file.
if exist('fast_blf.mexw32','file') == 0 && ...
        exist('fast_blf.mexw64','file') == 0
    mex -g fast_blf.cpp
    disp('Compile finished')
end

% Apply grayscale  bilateral filtering,use c/c++ code. 
%sigma=[sigma_s,space_r]
B = fast_blf(A,sigma(1),sigma(2));



