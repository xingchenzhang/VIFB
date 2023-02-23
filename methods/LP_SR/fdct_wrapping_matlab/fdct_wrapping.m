function C = fdct_wrapping(x, is_real, finest, nbscales, nbangles_coarse)

% fdct_wrapping.m - Fast Discrete Curvelet Transform via wedge wrapping - Version 1.0
%
% Inputs
%   x           M-by-N matrix
%
% Optional Inputs
%   is_real     Type of the transform
%                   0: complex-valued curvelets
%                   1: real-valued curvelets
%               [default set to 0]
%   finest      Chooses one of two possibilities for the coefficients at the
%               finest level:
%                   1: curvelets
%                   2: wavelets
%               [default set to 2]
%   nbscales    number of scales including the coarsest wavelet level
%               [default set to ceil(log2(min(M,N)) - 3)]
%   nbangles_coarse
%               number of angles at the 2nd coarsest level, minimum 8,
%               must be a multiple of 4. [default set to 16]
%
% Outputs
%   C           Cell array of curvelet coefficients.
%               C{j}{l}(k1,k2) is the coefficient at
%                   - scale j: integer, from finest to coarsest scale,
%                   - angle l: integer, starts at the top-left corner and
%                   increases clockwise,
%                   - position k1,k2: both integers, size varies with j
%                   and l.
%               If is_real is 1, there are two types of curvelets,
%               'cosine' and 'sine'. For a given scale j, the 'cosine'
%               coefficients are stored in the first two quadrants (low
%               values of l), the 'sine' coefficients in the last two
%               quadrants (high values of l).  
%
% See also ifdct_wrapping.m, fdct_wrapping_param.m
%
% By Laurent Demanet, 2004

X = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)));
[N1,N2] = size(X);
if nargin < 2, is_real = 0; end;
if nargin < 3, finest = 2; end;
if nargin < 4, nbscales = ceil(log2(min(N1,N2)) - 3); end;
if nargin < 5, nbangles_coarse = 16; end;

% Initialization: data structure
nbangles = [1, nbangles_coarse .* 2.^(ceil((nbscales-(nbscales:-1:2))/2))];
if finest == 2, nbangles(nbscales) = 1; end;
C = cell(1,nbscales);
for j = 1:nbscales
    C{j} = cell(1,nbangles(j));
end;

% Loop: pyramidal scale decomposition
M1 = N1/3;
M2 = N2/3;
if finest == 1,

    % Initialization: smooth periodic extension of high frequencies
    bigN1 = 2*floor(2*M1)+1;
    bigN2 = 2*floor(2*M2)+1;
    equiv_index_1 = 1+mod(floor(N1/2)-floor(2*M1)+(1:bigN1)-1,N1);
    equiv_index_2 = 1+mod(floor(N2/2)-floor(2*M2)+(1:bigN2)-1,N2);
    X = X(equiv_index_1,equiv_index_2);
        % Invariant: equiv_index_1(floor(2*M1)+1) == (N1 + 2 - mod(N1,2))/2
        % is the center in frequency. Same for M2, N2.
    window_length_1 = floor(2*M1) - floor(M1) - 1 - (mod(N1,3)==0);
    window_length_2 = floor(2*M2) - floor(M2) - 1 - (mod(N2,3)==0);
        % Invariant: floor(M1) + floor(2*M1) == N1 - (mod(M1,3)~=0)
        % Same for M2, N2.
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    if mod(N1,3)==0, lowpass_1 = [0, lowpass_1, 0]; end;
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    if mod(N2,3)==0, lowpass_2 = [0, lowpass_2, 0]; end;
    lowpass = lowpass_1'*lowpass_2;
    Xlow = X .* lowpass;

    scales = nbscales:-1:2;

else
    
    M1 = M1/2;
    M2 = M2/2;
    window_length_1 = floor(2*M1) - floor(M1) - 1;
    window_length_2 = floor(2*M2) - floor(M2) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    lowpass = lowpass_1'*lowpass_2;
    hipass = sqrt(1 - lowpass.^2);
    Xlow_index_1 = ((-floor(2*M1)):floor(2*M1)) + ceil((N1+1)/2);
    Xlow_index_2 = ((-floor(2*M2)):floor(2*M2)) + ceil((N2+1)/2);
    Xlow = X(Xlow_index_1, Xlow_index_2) .* lowpass;
    Xhi = X;
    Xhi(Xlow_index_1, Xlow_index_2) = Xhi(Xlow_index_1, Xlow_index_2) .* hipass;
    C{nbscales}{1} = fftshift(ifft2(ifftshift(Xhi)))*sqrt(prod(size(Xhi)));
    if is_real, C{nbscales}{1} = real(C{nbscales}{1}); end;
    
    scales = (nbscales-1):-1:2;

end;
for j = scales,

    M1 = M1/2;
    M2 = M2/2;
    window_length_1 = floor(2*M1) - floor(M1) - 1;
    window_length_2 = floor(2*M2) - floor(M2) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    lowpass = lowpass_1'*lowpass_2;
    hipass = sqrt(1 - lowpass.^2);
    Xhi = Xlow;                 % size is 2*floor(4*M1)+1 - by - 2*floor(4*M2)+1
    Xlow_index_1 = ((-floor(2*M1)):floor(2*M1)) + floor(4*M1) + 1;
    Xlow_index_2 = ((-floor(2*M2)):floor(2*M2)) + floor(4*M2) + 1;
    Xlow = Xlow(Xlow_index_1, Xlow_index_2);
    Xhi(Xlow_index_1, Xlow_index_2) = Xlow .* hipass;
    Xlow = Xlow .* lowpass;     % size is 2*floor(2*M1)+1 - by - 2*floor(2*M2)+1
    
    % Loop: angular decomposition
    l = 0;
    nbquadrants = 2 + 2*(~is_real);
    nbangles_perquad = nbangles(j)/4;
    for quadrant = 1:nbquadrants
        M_horiz = M2 * (mod(quadrant,2)==1) + M1 * (mod(quadrant,2)==0);
        M_vert = M1 * (mod(quadrant,2)==1) + M2 * (mod(quadrant,2)==0);
        if mod(nbangles_perquad,2),
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right(end:-1:1)];
        else
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right((end-1):-1:1)];
        end;
        wedge_endpoints = wedge_ticks(2:2:(end-1));         % integers
        wedge_midpoints = (wedge_endpoints(1:(end-1)) + wedge_endpoints(2:end))/2;
                % integers or half-integers
        
        % Left corner wedge
        l = l+1;
        first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles_perquad) + 1);
        length_corner_wedge = floor(4*M_vert) - floor(M_vert) + ceil(first_wedge_endpoint_vert/4);
        Y_corner = 1:length_corner_wedge;
        [XX,YY] = meshgrid(1:(2*floor(4*M_horiz)+1),Y_corner);
        width_wedge = wedge_endpoints(2) + wedge_endpoints(1) - 1;
        slope_wedge = (floor(4*M_horiz) + 1 - wedge_endpoints(1))/floor(4*M_vert);
        left_line = round(2 - wedge_endpoints(1) + slope_wedge*(Y_corner - 1));
                                                            % integers
        [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
                % Coordinates of the top-left corner of the wedge wrapped
                % around the origin. Some subtleties when the wedge is
                % even-sized because of the forthcoming 90 degrees rotation
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_data(new_row,:) = Xhi(row,admissible_cols) .* (cols > 0);
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;
        slope_wedge_right = (floor(4*M_horiz)+1 - wedge_midpoints(1))/floor(4*M_vert);
        mid_line_right = wedge_midpoints(1) + slope_wedge_right*(wrapped_YY - 1);
                % not integers in general
        coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(2) - wedge_endpoints(1)) * ...
            (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
        C2 = 1/(1/(2*(floor(4*M_horiz))/(wedge_endpoints(1) - 1) - 1) + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
        C1 = C2 / (2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) + 1;
        coord_corner = C1 + C2 * ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert))) ./ ...
            (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert))));
        wl_left = fdct_wrapping_window(coord_corner);
        [wl_right,wr_right] = fdct_wrapping_window(coord_right);
        wrapped_data = wrapped_data .* (wl_left .* wr_right);

        switch is_real
            case 0
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
            case 1
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                C{j}{l} = sqrt(2)*real(x);
                C{j}{l+nbangles(j)/2} = sqrt(2)*imag(x);
        end;
                
        % Regular wedges
        length_wedge = floor(4*M_vert) - floor(M_vert);
        Y = 1:length_wedge;
        first_row = floor(4*M_vert)+2-ceil((length_wedge+1)/2)+...
            mod(length_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        for subl = 2:(nbangles_perquad-1);
            l = l+1;
            width_wedge = wedge_endpoints(subl+1) - wedge_endpoints(subl-1) + 1;
            slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(subl))/floor(4*M_vert);
            left_line = round(wedge_endpoints(subl-1) + slope_wedge*(Y - 1));
            [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_wedge,width_wedge));
            first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
                mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
                wrapped_data(new_row,:) = Xhi(row,cols);
                wrapped_XX(new_row,:) = XX(row,cols);
                wrapped_YY(new_row,:) = YY(row,cols);             
            end;
            slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(subl-1))/floor(4*M_vert);
            mid_line_left = wedge_midpoints(subl-1) + slope_wedge_left*(wrapped_YY - 1);
            coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl) - wedge_endpoints(subl-1)) * ...
                (wrapped_XX - mid_line_left)./(floor(4*M_vert)+1 - wrapped_YY);
            slope_wedge_right = ((floor(4*M_horiz)+1) - wedge_midpoints(subl))/floor(4*M_vert);
            mid_line_right = wedge_midpoints(subl) + slope_wedge_right*(wrapped_YY - 1);
            coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl+1) - wedge_endpoints(subl)) * ...
                (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
            wl_left = fdct_wrapping_window(coord_left);
            [wl_right,wr_right] = fdct_wrapping_window(coord_right);
            wrapped_data = wrapped_data .* (wl_left .* wr_right);
            switch is_real
                case 0
                    wrapped_data = rot90(wrapped_data,-(quadrant-1));
                    C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                case 1
                    wrapped_data = rot90(wrapped_data,-(quadrant-1));
                    x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                    C{j}{l} = sqrt(2)*real(x);
                    C{j}{l+nbangles(j)/2} = sqrt(2)*imag(x);
            end;
        end;

        % Right corner wedge
        l = l+1;
        width_wedge = 4*floor(4*M_horiz) + 3 - wedge_endpoints(end) - wedge_endpoints(end-1);
        slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(end))/floor(4*M_vert);
        left_line = round(wedge_endpoints(end-1) + slope_wedge*(Y_corner - 1));
        [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_data(new_row,:) = Xhi(row,admissible_cols) .* (cols <= (2*floor(4*M_horiz)+1));
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;
        slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(end))/floor(4*M_vert);
        mid_line_left = wedge_midpoints(end) + slope_wedge_left*(wrapped_YY - 1);
        coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(end) - wedge_endpoints(end-1)) * ...
            (wrapped_XX - mid_line_left)./(floor(4*M_vert) + 1 - wrapped_YY);
        C2 = -1/(2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1 + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
        C1 = -C2 * (2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY - 1)/floor(4*M_vert)) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY - 1)/floor(4*M_vert)) - 1;
        coord_corner = C1 + C2 * (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert)))) ./ ...
            ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert)));
        wl_left = fdct_wrapping_window(coord_left);
        [wl_right,wr_right] = fdct_wrapping_window(coord_corner);

        wrapped_data = wrapped_data .* (wl_left .* wr_right);
        switch is_real
            case 0
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
            case 1
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                C{j}{l} = sqrt(2)*real(x);
                C{j}{l+nbangles(j)/2} = sqrt(2)*imag(x);
        end;

        if quadrant < nbquadrants, Xhi = rot90(Xhi); end;
    end;
end;

% Coarsest wavelet level
C{1}{1} = fftshift(ifft2(ifftshift(Xlow)))*sqrt(prod(size(Xlow)));
if is_real == 1,
    C{1}{1} = real(C{1}{1});
end;
