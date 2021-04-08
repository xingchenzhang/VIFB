function F = fusion_method(img1, img2, L, unit, de_level, norm, is_overlap)
% fusion method for detail parts and base parts

% img1 = imread(path1);
% img1 = im2double(img1);
% % if size(img1,3)==3
% %     img1=rgb2gray(img1);
% % end
% img2 = imread(path2);
% img2 = im2double(img2);
% % if size(img2,3)==3
% %     img2=rgb2gray(img2);
% % end

% deep decomposition
[I_b1,I_d1_deep, show_matrix1, count1] = deep_latent(img1, de_level, L, unit, is_overlap);
[I_b2,I_d2_deep, show_matrix2, count2] = deep_latent(img2, de_level, L, unit, is_overlap);
count_all = count1;

% fuion for base parts - average strategy
I_bf = fuison_base_parts(I_b1, I_b2);
% figure;imshow(I_bf);

% fuion for salient parts - l1-norm and nuclear-norm
[t1,t2] = size(img1);
if is_overlap==1
    % for overlapping
    I_d_temp = zeros(t1,t2);
    for i=1:de_level
        temp1 = I_d1_deep(:,:,i);
        temp2 = I_d2_deep(:,:,i);
        temp_vector = fuison_detail_parts(temp1, temp2, norm);
        % reconstruction
        for ii=1:t1-unit+1
            for jj=1:t2-unit+1
                temp = temp_vector(:,(ii-1)*(t2-unit+1)+jj);
                I_d_temp(ii:(unit+ii-1), jj:(unit+jj-1)) = I_d_temp(ii:(unit+ii-1), jj:(unit+jj-1)) + reshape(temp, [unit unit]);
            end
        end
        % average operation for overlapping position
        I_df_temp(:,:,i) = I_d_temp./count_all(:,:,i);
    end
    I_df = sum(I_df_temp,3);
else
    % for distinct
    if de_level>1
        temp1 = sum(I_d1_deep,3);
        temp2 = sum(I_d2_deep,3);
    else
        temp1 = I_d1_deep;
        temp2 = I_d2_deep;
    end

    I_df_vector = fuison_detail_parts(temp1, temp2, norm, isZCA);

    [t1,t2] = size(img1);
    I_df = col2im(I_df_vector,[unit, unit],[t1,t2],'distinct');
    I_df(I_df<0) = 0;
end
% figure;imshow(I_df);

% reconstructed the fused image
F = I_bf + I_df;
% figure;imshow(F);

end



