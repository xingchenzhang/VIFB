function y = cimage5(x,xsat,i0,i1)

% function y = cimage5(x,xsat,i0,i1)
% Plot a complex matrix x as a colour image, using Matlab5 RGB mode.
% xsat, if given, defines the amplitude of x which corresponds
% to saturated colours.  (If xsat<0, xmax is used.)
% i0,i1, if given, define the grey and sat. colour intensities.
% i0 is the minimum grey intensity at the centre of the circle
% and i1 is the intensity at the circumference.
% The colours are mapped around the circle as:
%           yellow
%       ?        orange
%   green            red
%      cyan     magenta
%            blue
% 
% If an output is not requested, the image is plotted.
% Otherwise the image of map indices is returned as y. 
%
% For a plot of the colour circle with a white centre, use:
% t=ones(33,1)*[-16:16]; cimage5(t-j*t.',16,0.9,1.0)
% For a dark centre use:
% t=ones(33,1)*[-16:16]; cimage5(t-j*t.',16,0.3,1.2)

% Nick Kingsbury, Cambridge University, July 1998.

if nargin < 4, i1 = 1.0; end
if nargin < 3, i0 = 0.9; end

xmax = max(max(abs(x)));
if nargin < 2, xsat = xmax; end
if xsat<0, xsat = xmax; end

fprintf(1,'Fig %.0f: xsat = %f, xmax = %f\n', gcf, xsat, xmax);

% Scale x and limit to a max amplitude of 1.
x = x ./ xsat;
x = x ./ max(abs(x),1);

% Calculate red, green and blue intensities in cx.
ax = (i0 - 0.5*i1) * (1 - abs(x));
cx(:,:,1) = (0.5 * i1) * (real(x) + 1) + ax;
cx(:,:,2) = (0.25 * i1) * (2 - real(x) + imag(x)) + ax;
cx(:,:,3) = (0.5 * i1) * (1 - imag(x)) + ax;
cx = max(min(cx,1),0);

if nargout == 0,
% Draw the image.
  image(cx)
  axis image
else
  y = cx;
end

return


