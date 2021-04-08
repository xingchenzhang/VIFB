function f = vector_fuison_detail_parts(s1_or, s2_or, norm)
% fousion strtegy for detail parts with different norms

[w,h,c] = size(s1_or);
unit = sqrt(w);

s1 = s1_or;
s2 = s2_or;

if norm == 'l1'
    % l1-norm
    s1_l1 = sum(abs(s1),1);
    s2_l1 = sum(abs(s2),1);
    
    s1_norm = s1_l1;
    s2_norm = s2_l1;
elseif norm == 'nu'
    % nuclear-norm
    s1_nu = zeros(1,h);
    s2_nu = zeros(1,h);
    for i=1:h
        temp1 = reshape(s1(:,i), [unit unit]);
        [U, S, V] = svd(temp1, 'econ');
        nu_norm1 = sum(diag(S));
        s1_nu(:,i)=nu_norm1;
        temp2 = reshape(s2(:,i), [unit unit]);
        [U, S, V] = svd(temp2, 'econ');
        nu_norm2 = sum(diag(S));
        s2_nu(:,i)=nu_norm2;
    end
    s1_norm = s1_nu;
    s2_norm = s2_nu;
end

w1_value = s1_norm./(s1_norm+s2_norm);
w2_value = s2_norm./(s1_norm+s2_norm);

w1 = repmat(w1_value, w, 1);
w2 = repmat(w2_value, w, 1);
    
f = w1.*s1_or+w2.*s2_or;
% figure;imshow(f);

end