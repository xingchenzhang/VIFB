function Y = coldwtfilt(X, h, offset)
% function Y = coldwtfilt(X, h, offset)
% Filter the columns of image X using filter vector h, decimating by 2.
% If length(h) is odd, each output sample is aligned with samples 1:2:(size(x,1)-1)
% if offset = 0, and with samples 2:2:size(x,1) if offset = 1.
% If length(h) is even, each output sample is aligned with the mid point
% of each pair of input samples.
% size(Y) = size(X) .* [0.5 1]
% The length of the columns in X must be even.
%
% This routine is inverted by coliwtfilt.m .
% To test these routines with biorthogonal wavelet filters use:
%
% xo=ones(7,2);x=[-xo;-1 -1;1 -1;xo];load antonini
% y0=coldwtfilt(x,h0o,0);y1=coldwtfilt(x,h1o,1);
% z=coliwtfilt(y0,2*g0o,0)+coliwtfilt(y1,2*g1o,1);max(abs(z(:)-x(:)))
%
% The result should be < 1e-12 for 'perfect reconstruction'.
%
% Nick Kingsbury, Cambridge University, May 2002

[r,c] = size(X);
if rem(r,2) > 0
   error('No of rows in X must be a multiple of 2!');
end

m = length(h);
m2 = fix(m/2);

if any(X(:))
   if rem(m,2),  % Odd length h.
      % Symmetrically extend with no repeat of end samples.
      xe = reflect([(1-m2):(r+m2)], 1, r); % Use 'reflect' so r < m2 works OK.
      t1 = [1:2:(r+m-2)] + offset; 
      t2 = [2:2:(r+m-3)] + offset;
   else          % Even length h.
      % Symmetrically extend with repeat of end samples.
      xe = reflect([(1-m2):(r+m2)], 0.5, r+0.5); % Use 'reflect' so r < m2 works OK.
      t2 = 2:2:(r+m-2); 
      t1 = 3:2:(r+m-1);
   end 
   
   % Separate h into 2 filters comprising odd and even taps of h,
   % so that convolutions with odd and even samples of X can be done
   % separately for maximum efficiency.
   ho = h(1:2:m); 
   he = h(2:2:m);
   
   % Perform filtering separately on the odd and even samples of the 
   % columns of extended matrix X(xe,:), keeping only the 'valid' output 
   % samples, so Y is half the height of X. 
   Y = conv2(X(xe(t1),:),ho(:),'valid') + conv2(X(xe(t2),:),he(:),'valid');
else
   Y = zeros(r/2,c);
end
return;
