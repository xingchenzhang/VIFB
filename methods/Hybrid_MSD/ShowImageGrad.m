% ShowImageGrad.m
% -------------------------------------------------------------------
% 
% Authors: Sun Li
% Date:    22/03/2013
% Last modified: 06/05/2013
% -------------------------------------------------------------------

function ShowImageGrad(img, para1, gra, para2)

    if exist('para1', 'var') && isfield(para1, 'fig'),
        figure('NumberTitle', 'off', 'Name', para1.fig);
    else
        figure;
    end

    if max(img(:)) > 1+eps,
        img = max(img, 0);
        img = min(img, 255);
    else
        img = max(img, 0);
        img = min(img, 1);
        img = img*255;
    end
    
    % ------------ Opposite -------------
    if exist('para1', 'var') && isfield(para1, 'opp'),
        if para1.opp
            img = 255-img;
        end
    end
    % -----------------------------------
    
    imshow(uint8(floor(img)))
    
    if exist('para1', 'var') && isfield(para1, 'title'),
        title(para1.title);
    end
    
    if nargin < 3,
        return;
    end
    
    % Add the Gradient --------------------------
    hold on;
    step = 7;
    dddiv = 5;
    lar = 0;
    col = 'red';
    if nargin == 4,
        if isfield(para2, 'step'),
            step = para2.step;
        end
        if isfield(para2, 'dddiv'),
            dddiv = para2.dddiv;
        end
        if isfield(para2, 'color'),
            col = para2.color;
        end
        if isfield(para2, 'lar'),
            lar = para2.lar;
        end
    end
    
    [HH, WW] = size(img);
    [XX, YY] = meshgrid(1:step:WW, 1:step:HH);
    
    
    u=real(gra(1:step:HH, 1:step:WW)) / dddiv;    
    v=imag(gra(1:step:HH, 1:step:WW)) / dddiv;
    h=quiver(XX,YY,u,v, lar);
    set(h,'color',col);
    hold off
    
    
end