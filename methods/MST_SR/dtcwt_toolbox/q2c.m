function z = q2c(y)

% function z = q2c(y)
% Convert from quads in y to complex numbers in z.

sy = size(y);
t1 = 1:2:sy(1); t2 = 1:2:sy(2);
j2 = sqrt([0.5 -0.5]);

% Arrange pixels from the corners of the quads into
% 2 subimages of alternate real and imag pixels.
%  a----b
%  |    |
%  |    |
%  c----d

% Combine (a,b) and (d,c) to form two complex subimages. 
p = y(t1,t2)*j2(1) + y(t1,t2+1)*j2(2);     % p = (a + jb) / sqrt(2)
q = y(t1+1,t2+1)*j2(1) - y(t1+1,t2)*j2(2); % q = (d - jc) / sqrt(2)

% Form the 2 subbands in z.
z = cat(3,p-q,p+q);

return
