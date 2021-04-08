function x = ifdct_wrapping(C, is_real, M, N)

% ifdct_wrapping.m - Inverse Fast Discrete Curvelet Transform via wedge wrapping - Version 1.0
% This is in fact the adjoint, also the pseudo-inverse
%
% Inputs
%   C           Cell array containing curvelet coefficients (see
%               description in fdct_wrapping.m)
%   is_real     As used in fdct_wrapping.m
%   M, N        Size of the image to be recovered (not necessary if finest
%               = 2)
%
% Outputs
%   x           M-by-N matrix
%
% See also fdct_wrapping.m
%
% By Laurent Demanet, 2004

% Initialization
nbscales = length(C);
nbangles_coarse = length(C{2});
nbangles = [1, nbangles_coarse .* 2.^(ceil((nbscales-(nbscales:-1:2))/2))];
if length(C{end}) == 1, finest = 2; else finest = 1; end;
if finest == 2, nbangles(nbscales) = 1; end;
if nargin < 2, is_real = 0; end;
if nargin < 4,
    if finest == 1, error('Syntax: IFCT_wrapping(C,M,N) where the matrix to be recovered is M-by-N'); end;
    [N1,N2] = size(C{end}{1});
else
    N1 = M;
    N2 = N;
end;

M1 = N1/3;
M2 = N2/3;

if finest == 1;
    
    bigN1 = 2*floor(2*M1)+1;
    bigN2 = 2*floor(2*M2)+1;
    X = zeros(bigN1,bigN2);

    % Initialization: preparing the lowpass filter at finest scale
    window_length_1 = floor(2*M1) - floor(M1) - 1 - (mod(N1,3)==0);
    window_length_2 = floor(2*M2) - floor(M2) - 1 - (mod(N2,3)==0);
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    if mod(N1,3)==0, lowpass_1 = [0, lowpass_1, 0]; end;
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    if mod(N2,3)==0, lowpass_2 = [0, lowpass_2, 0]; end;
    lowpass = lowpass_1'*lowpass_2;

    scales = nbscales:-1:2;
   
else

    M1 = M1/2;
    M2 = M2/2;
    
    bigN1 = 2*floor(2*M1)+1;
    bigN2 = 2*floor(2*M2)+1;
    X = zeros(bigN1,bigN2);
    
    window_length_1 = floor(2*M1) - floor(M1) - 1;
    window_length_2 = floor(2*M2) - floor(M2) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    lowpass = lowpass_1'*lowpass_2;
    hipass_finest = sqrt(1 - lowpass.^2);
    
    scales = (nbscales-1):-1:2;
    
end;

% Loop: pyramidal reconstruction

Xj_topleft_1 = 1;
Xj_topleft_2 = 1;
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
    lowpass_next = lowpass_1'*lowpass_2;
    hipass = sqrt(1 - lowpass_next.^2);
    Xj = zeros(2*floor(4*M1)+1,2*floor(4*M2)+1);
    
    % Loop: angles
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
        
        % Left corner wedge
        
        l = l+1;
        first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles_perquad) + 1);
        length_corner_wedge = floor(4*M_vert) - floor(M_vert) + ceil(first_wedge_endpoint_vert/4);
        Y_corner = 1:length_corner_wedge;
        [XX,YY] = meshgrid(1:(2*floor(4*M_horiz)+1),Y_corner);
        width_wedge = wedge_endpoints(2) + wedge_endpoints(1) - 1;
        slope_wedge = (floor(4*M_horiz) + 1 - wedge_endpoints(1))/floor(4*M_vert);
        left_line = round(2 - wedge_endpoints(1) + slope_wedge*(Y_corner - 1));
        [wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;

        slope_wedge_right = (floor(4*M_horiz)+1 - wedge_midpoints(1))/floor(4*M_vert);
        mid_line_right = wedge_midpoints(1) + slope_wedge_right*(wrapped_YY - 1);
                                                            % not integers
                                                            % in general
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
        switch is_real
         case 0
          wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(prod(size(C{j}{l})));
          wrapped_data = rot90(wrapped_data,(quadrant-1));
         case 1
          x = C{j}{l} + sqrt(-1)*C{j}{l+nbangles(j)/2};
          wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)))/sqrt(2);
          wrapped_data = rot90(wrapped_data,(quadrant-1));
        end;
        wrapped_data = wrapped_data .* (wl_left .* wr_right);
 
        % Unwrapping data
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            Xj(row,admissible_cols) = Xj(row,admissible_cols) + wrapped_data(new_row,:);
                                % We use the following property: in an assignment
                                % A(B) = C where B and C are vectors, if
                                % some value x repeats in B, then the
                                % last occurrence of x is the one
                                % corresponding to the eventual assignment.
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
            [wrapped_XX, wrapped_YY] = deal(zeros(length_wedge,width_wedge));
            first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
                mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
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
            switch is_real
             case 0
              wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(prod(size(C{j}{l})));
              wrapped_data = rot90(wrapped_data,(quadrant-1));
             case 1
              x = C{j}{l} + sqrt(-1)*C{j}{l+nbangles(j)/2};
              wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)))/sqrt(2);
              wrapped_data = rot90(wrapped_data,(quadrant-1));
            end;
            wrapped_data = wrapped_data .* (wl_left .* wr_right);
            
            % Unwrapping data
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
                Xj(row,cols) = Xj(row,cols) + wrapped_data(new_row,:);
            end;

        end;    % for subl
        
        % Right corner wedge
        l = l+1;
        width_wedge = 4*floor(4*M_horiz) + 3 - wedge_endpoints(end) - wedge_endpoints(end-1);
        slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(end))/floor(4*M_vert);
        left_line = round(wedge_endpoints(end-1) + slope_wedge*(Y_corner - 1));
        [wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);        
        end;
        YY = Y_corner'*ones(1,width_wedge);
        slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(end))/floor(4*M_vert);
        mid_line_left = wedge_midpoints(end) + slope_wedge_left*(wrapped_YY - 1);
        coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(end) - wedge_endpoints(end-1)) * ...
            (wrapped_XX - mid_line_left)./(floor(4*M_vert)+1 - wrapped_YY);
        C2 = -1/(2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1 + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
        C1 = -C2 * (2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY-1)/floor(4*M_vert)) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY-1)/floor(4*M_vert)) - 1;
        coord_corner = C1 + C2 * (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert)))) ./ ...
            ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert)));
        wl_left = fdct_wrapping_window(coord_left);
        [wl_right,wr_right] = fdct_wrapping_window(coord_corner);
        switch is_real
         case 0
          wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(prod(size(C{j}{l})));
          wrapped_data = rot90(wrapped_data,(quadrant-1));
         case 1
          x = C{j}{l} + sqrt(-1)*C{j}{l+nbangles(j)/2};
          wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)))/sqrt(2);
          wrapped_data = rot90(wrapped_data,(quadrant-1));
        end;
        wrapped_data = wrapped_data .* (wl_left .* wr_right);
        
         % Unwrapping data
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            Xj(row,fliplr(admissible_cols)) = Xj(row,fliplr(admissible_cols)) + wrapped_data(new_row,end:-1:1);
                                % We use the following property: in an assignment
                                % A(B) = C where B and C are vectors, if
                                % some value x repeats in B, then the
                                % last occurrence of x is the one
                                % corresponding to the eventual assignment.
        end;

        Xj = rot90(Xj);
        
    end;    % for quadrant
    
    Xj = Xj .* lowpass;
    Xj_index1 = ((-floor(2*M1)):floor(2*M1)) + floor(4*M1) + 1;
    Xj_index2 = ((-floor(2*M2)):floor(2*M2)) + floor(4*M2) + 1;
    Xj(Xj_index1, Xj_index2) = Xj(Xj_index1, Xj_index2) .* hipass;
    
    loc_1 = Xj_topleft_1 + (0:(2*floor(4*M1)));
    loc_2 = Xj_topleft_2 + (0:(2*floor(4*M2)));
    X(loc_1,loc_2) = X(loc_1,loc_2) + Xj;

    % Preparing for loop reentry or exit
    
    Xj_topleft_1 = Xj_topleft_1 + floor(4*M1) - floor(2*M1);
    Xj_topleft_2 = Xj_topleft_2 + floor(4*M2) - floor(2*M2);
    
    lowpass = lowpass_next;
    
end;    % for j

if is_real
    Y = X;
    X = rot90(X,2);
    X = X + conj(Y);
end
    
% Coarsest wavelet level
M1 = M1/2;
M2 = M2/2;
Xj = fftshift(fft2(ifftshift(C{1}{1})))/sqrt(prod(size(C{1}{1})));
loc_1 = Xj_topleft_1 + (0:(2*floor(4*M1)));
loc_2 = Xj_topleft_2 + (0:(2*floor(4*M2)));
X(loc_1,loc_2) = X(loc_1,loc_2) + Xj .* lowpass;

% Finest level
M1 = N1/3;
M2 = N2/3;
if finest == 1,

    % Folding back onto N1-by-N2 matrix
    shift_1 = floor(2*M1)-floor(N1/2);
    shift_2 = floor(2*M2)-floor(N2/2);
    Y = X(:,(1:N2)+shift_2);
    Y(:,N2-shift_2+(1:shift_2)) = Y(:,N2-shift_2+(1:shift_2)) + X(:,1:shift_2);
    Y(:,1:shift_2) = Y(:,1:shift_2) + X(:,N2+shift_2+(1:shift_2));
    X = Y((1:N1)+shift_1,:);
    X(N1-shift_1+(1:shift_1),:) = X(N1-shift_1+(1:shift_1),:) + Y(1:shift_1,:);
    X(1:shift_1,:) = X(1:shift_1,:) + Y(N1+shift_1+(1:shift_1),:);
    
else
    
    % Extension to a N1-by-N2 matrix
    Y = fftshift(fft2(ifftshift(C{nbscales}{1})))/sqrt(prod(size(C{nbscales}{1})));
    X_topleft_1 = ceil((N1+1)/2) - floor(M1);
    X_topleft_2 = ceil((N2+1)/2) - floor(M2);
    loc_1 = X_topleft_1 + (0:(2*floor(M1)));
    loc_2 = X_topleft_2 + (0:(2*floor(M2)));
    Y(loc_1,loc_2) = Y(loc_1,loc_2) .* hipass_finest + X;
    X = Y;
    
end;

x = fftshift(ifft2(ifftshift(X)))*sqrt(prod(size(X)));
if is_real, x = real(x); end;
