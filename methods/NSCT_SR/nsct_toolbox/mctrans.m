function h = mctrans(b,t)
% MCTRANS McClellan transformation
%   H = mctrans(B,T) produces the 2-D FIR filter H that
%   corresponds to the 1-D FIR filter B using the transform T.

% Convert the 1-D filter b to SUM_n a(n) cos(wn) form
n = (length(b)-1)/2;
b = rot90(fftshift(rot90(b,2)),2); % Inverse fftshift
a = [b(1) 2*b(2:n+1)];

inset = floor((size(t)-1)/2);

% Use Chebyshev polynomials to compute h
P0 = 1; P1 = t;
h = a(2)*P1; 
rows = inset(1)+1; cols = inset(2)+1;
h(rows,cols) = h(rows,cols)+a(1)*P0;
for i=3:n+1,
    P2 = 2*conv2(t,P1);
    rows = rows + inset(1); cols = cols + inset(2);
    P2(rows,cols) = P2(rows,cols) - P0;
    rows = inset(1) + [1:size(P1,1)];
    cols = inset(2) + [1:size(P1,2)];
    hh = h;
    h = a(i)*P2; h(rows,cols) = h(rows,cols) + hh;
    P0 = P1;
    P1 = P2;
end
h = rot90(h,2); % Rotate for use with filter2