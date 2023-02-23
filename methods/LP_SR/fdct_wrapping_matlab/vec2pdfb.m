function y = vec2pdfb(c, s)
% VEC2PDFB   Convert the vector form to the output structure of the PDFB
%
%       y = vec2pdfb(c, s)
%
% Input:
%   c:  1-D vector that contains all PDFB coefficients
%   s:  structure of PDFB output
%
% Output:
%   y:  PDFB coefficients in cell vector format that can be used in pdfbrec
%
% See also:	PDFB2VEC, PDFBREC

% Copy the coefficients from c to y according to the structure s
n = s(end, 1);      % number of pyramidal layers
y = cell(1, n);

% Variable that keep the current position
pos = prod(s(1, 3:4));

y{1} = reshape(c(1:pos), s(1, 3:4));

% Used for row index of s
ind = 1;

for l = 2:n
    % Number of directional subbands in this layer
    nd = length(find(s(:, 1) == l));

    y{l} = cell(1, nd);
    
    for d = 1:nd
        % Size of this subband
        p = s(ind + d, 3);
        q = s(ind + d, 4);
        ss = p * q;
        
        y{l}{d} = reshape(c(pos+[1:ss]), [p, q]);
        pos = pos + ss;
    end
    
    ind = ind + nd;
end