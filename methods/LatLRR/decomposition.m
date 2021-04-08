function [I_b, I_d, I_d_v, count_matric] = decomposition(img, L, unit, is_overlap)
% obtained the detail parts and base parts

% image
[t1,t2] = size(img);

% salient parts by project matrix
I_d = zeros(t1,t2);

% the matrices for overlapping
count_matric = zeros(t1,t2);
ones_matric = ones(t1,t2);

if is_overlap == 1
    temp_vector = [];
    for i=1:t1-unit+1
        for j=1:t2-unit+1
            temp = img(i:(unit+i-1), j:(unit+j-1));
            temp_vector(:,(i-1)*(t2-unit+1)+j) = temp(:);
            % record the overlapping number
            count_matric(i:(unit+i-1), j:(unit+j-1)) =...
                count_matric(i:(unit+i-1), j:(unit+j-1)) + ones_matric(i:(unit+i-1), j:(unit+j-1));
        end
    end
    % calculate salient features
    temp_vector = L*temp_vector;
    I_d_v = temp_vector;
    % reconstruction
    for i=1:t1-unit+1
        for j=1:t2-unit+1
            temp = temp_vector(:,(i-1)*(t2-unit+1)+j);
            I_d(i:(unit+i-1), j:(unit+j-1)) = I_d(i:(unit+i-1), j:(unit+j-1)) + reshape(temp, [unit unit]);
        end
    end
    % average operation for overlapping position
    I_d = I_d./count_matric;
else
    I_col = im2col(img, [unit, unit], 'distinct');
    salient = L*I_col;
    I_d_v = salient;

    I_d = col2im(salient,[unit, unit],[t1,t2],'distinct');
    I_d(I_d<0) = 0;
end
% base parts
I_b = img - I_d;

end