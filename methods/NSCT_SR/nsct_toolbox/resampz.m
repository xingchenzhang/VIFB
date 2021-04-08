function y = resampz(x, type, shift)
% RESAMPZ   Resampling of matrix
%
%       y = resampz(x, type, [shift])
%
% Input:
%       x:      input matrix
%       type:   one of {1, 2, 3, 4} (see note)
%       shift:  [optional] amount of shift (default is 1)
%
% Output:
%       y:      resampled matrix
%
% Note:
%	The resampling matrices are:
%		R1 = [1, 1;  0, 1];
%		R2 = [1, -1; 0, 1];
%		R3 = [1, 0;  1, 1];
%		R4 = [1, 0; -1, 1];
%
%	This resampling program does NOT involve periodicity, thus it
%	zero-pad and extend the matrix
%
% See also:	RESAMP

if ~exist('shift', 'var')
    shift = 1;
end

sx = size(x);

switch type
    case {1, 2}
	y = zeros(sx(1) + abs(shift * (sx(2) - 1)), sx(2));
	
	if type == 1
	    shift1 = [0:(sx(2)-1)] * (-shift);
	else
	    shift1 = [0:(sx(2)-1)] * shift;
	end
	
	% Adjust to non-negative shift if needed
	if (shift1(end) < 0)
	    shift1 = shift1 - shift1(end);
	end
	
	for n = 1:sx(2)
	    y(shift1(n)+(1:sx(1)), n) = x(:, n);
	end    
	
	% Finally, delete zero rows if needed
	start = 1;
	finish = size(y, 1);
	
	while norm(y(start, :)) == 0,
	    start = start + 1;
	end
	
	while norm(y(finish, :)) == 0,
	    finish = finish - 1;
	end
	
	y = y(start:finish, :);
	    
    case {3, 4}
	y = zeros(sx(1), sx(2) + abs(shift * (sx(1) - 1)));
	
	if type == 3
	    shift2 = [0:(sx(1)-1)] * (-shift);
	else
	    shift2 = [0:(sx(1)-1)] * shift;
	end
	
	% Adjust to non-negative shift if needed
	if (shift2(end) < 0)
	    shift2 = shift2 - shift2(end);
	end
	
	for m = 1:sx(1)
	    y(m, shift2(m)+(1:sx(2))) = x(m, :);
	end    
    
	% Finally, delete zero columns if needed
	start = 1;
	finish = size(y, 2);
	
	while norm(y(:, start)) == 0,
	    start = start + 1;
	end
	
	while norm(y(:, finish)) == 0,
	    finish = finish - 1;
	end
	
	y = y(:, start:finish);
    
    otherwise
	error('The second input (type) must be one of {1, 2, 3, 4}');	
end