% ZCA & l1-norm operation
function output = whitening_norm(features)
    features = im2double(features);
%     output = features;
    output = zca(features);
    output = l1_norm(output);
%     output = l2_norm(output);
%     output = nuclear(output);
end

% l1-norm
function output = l1_norm(features)
    [h,w,c] = size(features);
    output = zeros(h, w);
    unit = 5;
    % padding
    pad = (unit-1)/2;
    feature_temp = zeros(h+unit-1, w+unit-1, c);
    feature_temp(1+pad:h+pad, 1+pad:w+pad, :) = features;
    for i=1+pad:h+pad
        for j=1+pad:w+pad
            temp = feature_temp(i-pad:i+pad, j-pad:j+pad, :);
            norm = sum(sum(sum(abs(temp),3)))/(unit*unit);
            output(i-pad, j-pad) = norm;
        end
    end
end

% l2-norm
function output = l2_norm(features)
    [h,w,c] = size(features);
    output = zeros(h, w);
    unit = 5;
    % padding
    pad = (unit-1)/2;
    feature_temp = zeros(h+unit-1, w+unit-1, c);
    feature_temp(1+pad:h+pad, 1+pad:w+pad, :) = features;
    for i=1+pad:h+pad
        for j=1+pad:w+pad
            temp = feature_temp(i-pad:i+pad, j-pad:j+pad, :);
            temp = temp.*temp;
            norm = sum(sum(sqrt(sum(temp,3))))/(unit*unit);
            output(i-pad, j-pad) = norm;
        end
    end
end

% nuclear-norm
function output = nuclear(features)
    [h,w,c] = size(features);
    output = zeros(h, w);
    unit = 5;
    % padding
    pad = (unit-1)/2;
    feature_temp = zeros(h+unit-1, w+unit-1, c);
    feature_temp(1+pad:h+pad, 1+pad:w+pad, :) = features;
    for i=1+pad:h+pad
        for j=1+pad:w+pad
            temp = reshape(feature_temp(i-pad:i+pad, j-pad:j+pad, :), [unit*unit c]);
            [U, S, V] = svd(temp, 'econ');
            nu_norm = sum(diag(S));
            output(i-pad, j-pad) = nu_norm;
        end
    end

end



