function [dx, dy] = GradientSobel(img)
    
        Sv = [-1 -2 -1;...
               0  0  0;...
               1  2  1];
        Sh = [-1  0  1;
              -2  0  2;
              -1  0  1];
        dx=conv2(img, Sh,'valid'); 
        dy=conv2(img, Sv,'valid');
end
