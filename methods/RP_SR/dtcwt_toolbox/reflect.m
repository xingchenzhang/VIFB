function y = reflect(x,minx,maxx)

% function y = reflect(x,minx,maxx)
% Reflect the values in matrix x about the scalar values minx and maxx.
% Hence a vector x containing a long linearly increasing series
% is converted into a waveform which ramps linearly up and down
% between minx and maxx.
% If x contains integers and minx and maxx are (integers + 0.5),
% the ramps will have repeated max and min samples.
%
% Nick Kingsbury, Cambridge University, January 1999.

y = x;

% Reflect y in maxx.
t = find(y > maxx);
y(t) = 2*maxx - y(t);

% Reflect y in minx.
t = find(y < minx);
while ~isempty(t)
   y(t) = 2*minx - y(t);
   
   % Repeat until no more values out of range.
   t = find(y > maxx);
   if ~isempty(t), y(t) = 2*maxx - y(t); end
   t = find(y < minx);
end

return;

