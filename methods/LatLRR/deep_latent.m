function [I_b_deep,I_d_deep, show_matrix, count_all] = deep_latent(img, de_level, L, unit, is_overlap)
% decomposing the source images layer by layer
% de_level - the number of decomposition layers

[t1,t2] = size(img);
show_matrix = [];

I_b = zeros(t1,t2);
for i=1:de_level
    if i==1
        [temp_b, temp_d, temp_d_v, count_m] = decomposition(img, L, unit, is_overlap);
    else
        [temp_b, temp_d, temp_d_v, count_m] = decomposition(temp_b, L, unit, is_overlap);
    end
    I_d(:,:,i) = temp_d_v(:,:);
    count_all(:,:,i) = count_m(:,:);
    show_matrix = cat(2,show_matrix, temp_d);
end

show_matrix = cat(2,show_matrix, temp_b);
I_b(:,:) = temp_b(:,:);

I_b_deep = I_b;
I_d_deep = I_d;

end