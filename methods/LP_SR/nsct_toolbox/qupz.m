function y = qupz(x, type)
% QUPZ   Quincunx Upsampling (with zero-pad and matrix extending)
%
% 	y = qup(x, [type])
%
% Input:
%	x:	input image
%	type:	[optional] 1 or 2 for selecting the quincunx matrices:
%			Q1 = [1, -1; 1, 1] or Q2 = [1, 1; -1, 1]
% Output:
%	y:	qunincunx upsampled image
%
%       This resampling operation does NOT involve periodicity, thus it
%       zero-pad and extend the matrix

if ~exist('type', 'var')
    type = '1';
end

% Quincunx downsampling using the Smith decomposition:
%	Q1 = R2 * [2, 0; 0, 1] * R3
% and,
%	Q2 = R1 * [2, 0; 0, 1] * R4
%
% See RESAMP for the definition of those resampling matrices
%
% Note that R1 * R2 = R3 * R4 = I so for example,
% upsample by R1 is the same with down sample by R2.
% Also the order of upsampling operations is in the reserved order
% with the one of matrix multiplication.

switch type
    case 1
	x1 = resampz(x, 4);
	[m, n] = size(x1);
	x2 = zeros(2*m-1, n);
	x2(1:2:end, :) = x1;	
	y = resampz(x2, 1);
		
    case 2
	x1 = resampz(x, 3);
	[m, n] = size(x1);
	x2 = zeros(2*m-1, n);
	x2(1:2:end, :) = x1;	
	y = resampz(x2, 2);
			
    otherwise
	error('Invalid argument type');
end
