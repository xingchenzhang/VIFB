% -------------------------------------------------------------------
%
% Authors: Sun Li
% Date:    15/04/2013
% Last modified: 26/04/2013
% -------------------------------------------------------------------
function [img1, img2, name] = PickName(path1, path2, varargin)

    bGray = 1;
    if ~isempty(varargin),
        bGray = varargin{1};
    end
    img1 = imread(path1);
    img2 = imread(path2);
    
    if bGray,
        if size(img1, 3) ~= 1,
            img1 = rgb2gray(img1);
        end
         if size(img2, 3) ~= 1,
            img2 = rgb2gray(img2);
        end
    end
    
    img1 = double(img1);
    img2 = double(img2);
    
    [~, NAME, ~] = fileparts(path1);
%     name = NAME(1:end-1);
    name = NAME(1:end);
end