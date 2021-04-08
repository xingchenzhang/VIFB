% ZCA operation
function output = zca(features)
    output = zeros(size(features));
    if size(features,3)>1
        for i = 1:size(features,3)
            temp = features(:,:,i);
            if mean2(temp)~=0
                matrix_co = (temp-mean2(temp))*(temp-mean2(temp))';
                [U1,S1,V1] = svd(matrix_co, 'econ');
                S1_ = diag(1./sqrt(diag(S1)+0.00001));
                % whitening
                output(:,:,i) = (U1*S1_*U1')*(temp-mean2(temp));
            else
                output(:,:,i) = temp;
            end
        end
    else
        if mean2(features)~=0
            features_co = (features-mean2(features))*(features-mean2(features))';
            [U1,S1,V1] = svd(features_co, 'econ');
            S1_ = diag(1./sqrt(diag(S1)+0.00001));
            % whitening
            output = (U1*S1_*U1')*(features-mean2(features));
        else
            output = features;
        end
    end
end
