function Y = colifilt(X, ha, hb)

% function Y = colifilt(X, ha, hb)
% Filter the columns of image X using the two filters ha and hb = reverse(ha).
% ha operates on the odd samples of X and hb on the even samples.
% Both filters should be even length, and h should be approx linear phase with
% a quarter sample advance from its mid pt (ie |h(m/2)| > |h(m/2 + 1)|).
%
%                   ext       left edge                      right edge       ext
% Level 2:        !               |               !               |               !
% +q filt on x      b       b       a       a       a       a       b       b       
% -q filt on o          a       a       b       b       b       b       a       a
% Level 1:        !               |               !               |               !
% odd filt on .    b   b   b   b   a   a   a   a   a   a   a   a   b   b   b   b   
% odd filt on .      a   a   a   a   b   b   b   b   b   b   b   b   a   a   a   a
%
% The output is interpolated by two from the input sample rate and the results
% from the two filters, Ya and Yb, are interleaved to give Y.
% Symmetric extension with repeated end samples is used on the composite X
% columns before each filter is applied.
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000
% Modified to be fast if X = 0, May 2002.

[r,c] = size(X);
if rem(r,2) > 0
   error('No of rows in X must be a multiple of 2!');
end

m = length(ha);
if m ~= length(hb)
   error('Lengths of ha and hb must be the same!');
end

if rem(m,2) > 0
   error('Lengths of ha and hb must be even!');
end
m2 = fix(m/2);

Y = zeros(r*2,c);
if ~any(X(:)), return; end

if rem(m2,2) == 0,
   
   % m/2 is even, so set up t to start on d samples.
   % Set up vector for symmetric extension of X with repeated end samples.
   xe = reflect([(1-m2):(r+m2)], 0.5, r+0.5); % Use 'reflect' so r < m2 works OK.
   
   
   t = [4:2:(r+m)];
   if sum(ha.*hb) > 0,
      ta = t; tb = t - 1;
   else
      ta = t - 1; tb = t;
   end
   
   hao = ha(1:2:m);
   hae = ha(2:2:m);
   hbo = hb(1:2:m);
   hbe = hb(2:2:m);
   
   s = 1:4:(r*2);
   
   Y(s,:)   = conv2(X(xe(tb-2),:),hae(:),'valid');
   Y(s+1,:) = conv2(X(xe(ta-2),:),hbe(:),'valid');
   Y(s+2,:) = conv2(X(xe(tb),:),hao(:),'valid');
   Y(s+3,:) = conv2(X(xe(ta),:),hbo(:),'valid');
   
   
else
   
   % m/2 is odd, so set up t to start on b samples.
   % Set up vector for symmetric extension of X with repeated end samples.
   xe = reflect([(1-m2):(r+m2)], 0.5, r+0.5); % Use 'reflect' so r < m2 works OK.
   
   t = [3:2:(r+m-1)];
   if sum(ha.*hb) > 0,
      ta = t; tb = t - 1;
   else
      ta = t - 1; tb = t;
   end
   
   hao = ha(1:2:m);
   hae = ha(2:2:m);
   hbo = hb(1:2:m);
   hbe = hb(2:2:m);
   
   s = 1:4:(r*2);
   
   Y(s,:)   = conv2(X(xe(tb),:),hao(:),'valid');
   Y(s+1,:) = conv2(X(xe(ta),:),hbo(:),'valid');
   Y(s+2,:) = conv2(X(xe(tb),:),hae(:),'valid');
   Y(s+3,:) = conv2(X(xe(ta),:),hbe(:),'valid');
   
end

return