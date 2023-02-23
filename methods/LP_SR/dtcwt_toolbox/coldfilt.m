function Y = coldfilt(X, ha, hb)

% function Y = coldfilt(X, ha, hb)
% Filter the columns of image X using the two filters ha and hb = reverse(ha).
% ha operates on the odd samples of X and hb on the even samples.
% Both filters should be even length, and h should be approx linear phase with
% a quarter sample advance from its mid pt (ie |h(m/2)| > |h(m/2 + 1)|).
%
%                   ext        top edge                     bottom edge       ext
% Level 1:        !               |               !               |               !
% odd filt on .    b   b   b   b   a   a   a   a   a   a   a   a   b   b   b   b   
% odd filt on .      a   a   a   a   b   b   b   b   b   b   b   b   a   a   a   a
% Level 2:        !               |               !               |               !
% +q filt on x      b       b       a       a       a       a       b       b       
% -q filt on o          a       a       b       b       b       b       a       a
%
% The output is decimated by two from the input sample rate and the results
% from the two filters, Ya and Yb, are interleaved to give Y.
% Symmetric extension with repeated end samples is used on the composite X
% columns before each filter is applied.
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000

[r,c] = size(X);
if rem(r,4) > 0, error('No. of rows in X must be a multiple of 4!'); end
m = length(ha);
if m ~= length(hb), error('Lengths of ha and hb must be the same!'); end
if rem(m,2) > 0, error('Lengths of ha and hb must be even!'); end
m2 = fix(m/2);

% Set up vector for symmetric extension of X with repeated end samples.
xe = reflect([(1-m):(r+m)], 0.5, r+0.5); % Use 'reflect' so d < m works OK.

% Select odd and even samples from ha and hb.
hao = ha(1:2:m);
hae = ha(2:2:m);
hbo = hb(1:2:m);
hbe = hb(2:2:m);
t = [6:4:(r+2*m-2)];

r2 = r/2;
Y = zeros(r2,c);
if sum(ha.*hb) > 0,
   s1 = 1:2:r2; 
   s2 = s1 + 1;
else
   s2 = 1:2:r2; 
   s1 = s2 + 1;
end

% Perform filtering on columns of extended matrix X(xe,:) in 4 ways. 
Y(s1,:) = conv2(X(xe(t-1),:),hao(:),'valid') + conv2(X(xe(t-3),:),hae(:),'valid');
Y(s2,:) = conv2(X(xe(t),:),hbo(:),'valid') + conv2(X(xe(t-2),:),hbe(:),'valid');

return