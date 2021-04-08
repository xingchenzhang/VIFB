function [ diff_median_mask]= diff_median_blk( In, win, thres)

p = 2*win + 1;

median_val = medfilt2(In, [p p]);

diff_median_mat = abs( median_val - In);

diff_median_mask= diff_median_mat >= ( thres* max( diff_median_mat(:)));
