function y = modulate2(x, type, center)
% MODULATE2	2D modulation
%
%	y = modulate2(x, type, [center])
%
% With TYPE = {'r', 'c' or 'b'} for modulate along the row, or column or
% both directions.
%
% CENTER secify the origin of modulation as floor(size(x)/2)+1+center
% (default is [0, 0])

if ~exist('center', 'var')
    center = [0, 0];
end

% Size and origin
s = size(x);
o = floor(s / 2) + 1 + center;

n1 = [1:s(1)] - o(1);
n2 = [1:s(2)] - o(2);

switch lower(type(1))
    case 'r'
	m1 = (-1) .^ n1;
	y = x .* repmat(m1', [1, s(2)]);
	
    case 'c'
	m2 = (-1) .^ n2;
	y = x .* repmat(m2, [s(1), 1]);
	
    case 'b'
	m1 = (-1) .^ n1;
	m2 = (-1) .^ n2;
	m = m1' * m2;
	y = x .* m;
	
    otherwise
	error('Invalid input type');
end
