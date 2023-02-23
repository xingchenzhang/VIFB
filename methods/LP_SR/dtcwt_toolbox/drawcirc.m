function p = drawcirc(r,w,dx,dy,N)
% function p = drawcirc(r,w,dx,dy,N)
% Generate an image of size N*N pels, containing a circle
% radius r pels and centred at dx,dy relative
% to the centre of the image.  The edge of the circle is a cosine shaped
% edge of width w (from 10 to 90% points).

x = ones(N,1) * (([1:N] - (N+1)/2 - dx)/r);
y = (([1:N]' - (N+1)/2 - dy)/r) * ones(1,N);

p = 0.5 + 0.5 * sin(min(max((exp(-0.5 * (x .* x + y .* y)) - exp(-0.5))*(r*3/w), -pi/2), pi/2));

return


