% weighting
function [gen, weight_a, weight_b] = fusion_strategy(features_a, features_b, source_a, source_b)

[m1,n1] = size(source_a);
% resize
resize_temp1 = imresize(features_a, [m1, n1]);
resize_temp2 = imresize(features_b, [m1, n1]);
% soft-max
weight_ave_temp1 = resize_temp1./(resize_temp1+resize_temp2);
weight_ave_temp2 = resize_temp2./(resize_temp1+resize_temp2);
% reconstruction
gen = source_a.*weight_ave_temp1 + source_b.*weight_ave_temp2;
% figure;imshow(gen);

weight_a = weight_ave_temp1;
weight_b = weight_ave_temp2;
end