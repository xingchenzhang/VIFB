function [X_rows, X_cols, F_rows, F_cols, N_rows, N_cols] = fdct_wrapping_param(C,M,N)

% fdct_wrapping_param.m - Gives the location of each curvelet in phase space
%
% Inputs
%   C           Cell array containing curvelet coefficients (see
%               description in fdct_wrapping.m)
%   M, N        Size of the original image (not necessary if finest
%               = 2)
%
% Outputs
%   X_rows      Cell array X_rows{j}{l}(k1,k2) giving the row index on the
%               original M-by-N image grid, of the center of the curvelet
%               indexed by j, l, k1 and k2. Could be non-integer-valued.
%   X_cols      Same for the column index.
%   F_rows      Cell array F_rows{j}{l} giving the frequency row index of
%               all the curvelets indexed by j and l, i.e., their center in
%               frequency. Could be non-integer-valued.
%   F_cols      Same for the column index. The angle (mod pi/2) can then be
%               computed as (-1)*atan(F_rows{j}{l}/F_cols{j}{l})
%   N_rows      Cell array N_rows{j}{l} giving the row size of each array
%               C{j}{l}(k1,k2) computed by fdct_wrapping.
%   N_cols      Same for the column size.
%
% See also fdct_wrapping.m, ifdct_wrapping.m
%
% By Laurent Demanet, 2004

nbscales = length(C);
nbangles_coarse = length(C{2});
nbangles = [1, nbangles_coarse .* 2.^(ceil((nbscales-(nbscales:-1:2))/2))];
if length(C{end}) == 1, finest = 2; else finest = 1; end;
if finest == 2, nbangles(nbscales) = 1; end;
if nargin < 3,
    if finest == 1, error('Syntax: fdct_wrapping_param(C,M,N) where the original image x is M-by-N'); end;
    [N1,N2] = size(C{end}{1});
else
    N1 = M;
    N2 = N;
end;

[X_rows,X_cols,F_rows,F_cols,N_rows,N_cols] = deal(cell(1,nbscales));
for j = 1:nbscales
    [X_rows{j},X_cols{j},F_rows{j},F_cols{j},N_rows{j},N_cols{j}] = deal(cell(1,nbangles(j)));
end;

M1 = N1/3;
M2 = N2/3;

if finest == 1,
    scales = nbscales:-1:2;
else
    [len, width] = size(C{nbscales}{1});
    F_rows{nbscales}{1} = 0;
    F_cols{nbscales}{1} = 0;
    xloc = 1 + N2*(0:(1/width):(1-1/width));
    X_cols{nbscales}{1} = ones(len,1) * xloc;
    yloc = 1 + N1*(0:(1/len):(1-1/len));
    X_rows{nbscales}{1} = yloc' * ones(1,width);
    N_rows{nbscales}{1} = len;
    N_cols{nbscales}{1} = width;

    M1 = M1/2;
    M2 = M2/2;
    scales = (nbscales-1):-1:2;
end;

dyadic_tick_1 = - floor(2*M1) + floor(N1/2);
dyadic_tick_2 = - floor(2*M2) + floor(N2/2);
for j = scales,

    M1 = M1/2;
    M2 = M2/2;
    dyadic_tick_1 = dyadic_tick_1 + floor(4*M1) - floor(2*M1);
    dyadic_tick_2 = dyadic_tick_2 + floor(4*M2) - floor(2*M2);
   
    % Loop: angles
    l = 0;
    nbquadrants = 4;
    nbangles_perquad = nbangles(j)/nbquadrants;
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
        
        % Left corner wedge
        
        l = l+1;
        first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles_perquad) + 1);
        length_corner_wedge = floor(4*M_vert) - floor(M_vert) + ceil(first_wedge_endpoint_vert/4);
        width_wedge = wedge_endpoints(2) + wedge_endpoints(1) - 1;
        slope_wedge = (floor(4*M_horiz) + 1 - wedge_endpoints(1))/floor(4*M_vert);

        switch quadrant
            case 1
                w1 = dyadic_tick_1;
                w2 = floor(N2/2) + 1 + slope_wedge*(w1 - floor(N1/2) - 1);
            case 2
                w2 = N2 + 1 - dyadic_tick_2;
                w1 = floor(N1/2) + 1 - slope_wedge*(w2 - floor(N2/2) - 1);
            case 3
                w1 = N1 + 1 - dyadic_tick_1;
                w2 = floor(N2/2) + 1 + slope_wedge*(w1 - floor(N1/2) - 1);
            case 4
                w2 = dyadic_tick_2;
                w1 = floor(N1/2) + 1 - slope_wedge*(w2 - floor(N2/2) - 1);
        end;
        F_rows{j}{l} = w1 - ceil((N1+1)/2);
        F_cols{j}{l} = w2 - ceil((N2+1)/2);
        
        size_wedge_horiz = width_wedge * (mod(quadrant,2)==1) + length_corner_wedge * (mod(quadrant,2)==0);
        size_wedge_vert = width_wedge * (mod(quadrant,2)==0) + length_corner_wedge * (mod(quadrant,2)==1);
        xloc = 1 + N2*(0:(1/size_wedge_horiz):(1-1/size_wedge_horiz));
        X_cols{j}{l} = ones(size_wedge_vert,1) * xloc;
        yloc = 1 + N1*(0:(1/size_wedge_vert):(1-1/size_wedge_vert));
        X_rows{j}{l} = yloc' * ones(1,size_wedge_horiz);
        N_rows{j}{l} = size_wedge_vert;
        N_cols{j}{l} = size_wedge_horiz;
        
        % Regular wedges
        length_wedge = floor(4*M_vert) - floor(M_vert);
        for subl = 2:(nbangles_perquad-1);
            l = l+1;
            width_wedge = wedge_endpoints(subl+1) - wedge_endpoints(subl-1) + 1;
            slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(subl))/floor(4*M_vert);

            switch quadrant
                case 1
                    w1 = dyadic_tick_1;
                    w2 = floor(N2/2) + 1 + slope_wedge*(w1 - floor(N1/2) - 1);
                case 2
                    w2 = N2 + 1 - dyadic_tick_2;
                    w1 = floor(N1/2) + 1 - slope_wedge*(w2 - floor(N2/2) - 1);
                case 3
                    w1 = N1 + 1 - dyadic_tick_1;
                    w2 = floor(N2/2) + 1 + slope_wedge*(w1 - floor(N1/2) - 1);
                case 4
                    w2 = dyadic_tick_2;
                    w1 = floor(N1/2) + 1 - slope_wedge*(w2 - floor(N2/2) - 1);
            end;
            F_rows{j}{l} = w1 - ceil((N1+1)/2);
            F_cols{j}{l} = w2 - ceil((N2+1)/2);
            
            size_wedge_horiz = width_wedge * (mod(quadrant,2)==1) + length_wedge * (mod(quadrant,2)==0);
            size_wedge_vert = width_wedge * (mod(quadrant,2)==0) + length_wedge * (mod(quadrant,2)==1);
            xloc = 1 + N2*(0:(1/size_wedge_horiz):(1-1/size_wedge_horiz));
            X_cols{j}{l} = ones(size_wedge_vert,1) * xloc;
            yloc = 1 + N1*(0:(1/size_wedge_vert):(1-1/size_wedge_vert));
            X_rows{j}{l} = yloc' * ones(1,size_wedge_horiz);
            N_rows{j}{l} = size_wedge_vert;
            N_cols{j}{l} = size_wedge_horiz;
            
        end;    % for subl
        
        % Right corner wedge
        l = l+1;
        width_wedge = 4*floor(4*M_horiz) + 3 - wedge_endpoints(end) - wedge_endpoints(end-1);
        slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(end))/floor(4*M_vert);

        switch quadrant
            case 1
                w1 = dyadic_tick_1;
                w2 = floor(N2/2) + 1 + slope_wedge*(w1 - floor(N1/2) - 1);
            case 2
                w2 = N2 + 1 - dyadic_tick_2;
                w1 = floor(N1/2) + 1 - slope_wedge*(w2 - floor(N2/2) - 1);
            case 3
                w1 = N1 + 1 - dyadic_tick_1;
                w2 = floor(N2/2) + 1 + slope_wedge*(w1 - floor(N1/2) - 1);
            case 4
                w2 = dyadic_tick_2;
                w1 = floor(N1/2) + 1 - slope_wedge*(w2 - floor(N2/2) - 1);
        end;
        F_rows{j}{l} = w1 - ceil((N1+1)/2);
        F_cols{j}{l} = w2 - ceil((N2+1)/2);
        
        size_wedge_horiz = width_wedge * (mod(quadrant,2)==1) + length_corner_wedge * (mod(quadrant,2)==0);
        size_wedge_vert = width_wedge * (mod(quadrant,2)==0) + length_corner_wedge * (mod(quadrant,2)==1);
        xloc = 1 + N2*(0:(1/size_wedge_horiz):(1-1/size_wedge_horiz));
        X_cols{j}{l} = ones(size_wedge_vert,1) * xloc;
        yloc = 1 + N1*(0:(1/size_wedge_vert):(1-1/size_wedge_vert));
        X_rows{j}{l} = yloc' * ones(1,size_wedge_horiz);
        N_rows{j}{l} = size_wedge_vert;
        N_cols{j}{l} = size_wedge_horiz;
        
    end;    % for quadrant
    
end;    % for j

% Coarsest wavelet level
M1 = M1/2;
M2 = M2/2;
[len, width] = size(C{1}{1});
F_rows{1}{1} = 0;
F_cols{1}{1} = 0;
xloc = 1 + N2*(0:(1/width):(1-1/width));
X_cols{1}{1} = ones(len,1) * xloc;
yloc = 1 + N1*(0:(1/len):(1-1/len));
X_rows{1}{1} = yloc' * ones(1,width);
N_rows{1}{1} = len;
N_cols{1}{1} = width;
