function [ diff_median_mask]= diff_median( i01, window_side, thres)

i01_size= size( i01);
diff_median_mat= zeros( [ i01_size]);

for j= 1+ window_side: i01_size(1)- window_side                             % rows
    for k= 1+ window_side: i01_size(2)- window_side                         % columns
       
        i01_temp= i01( j- window_side: j+ window_side, k- window_side:...
            k+ window_side);
        center= i01_temp( window_side+ 1, window_side+ 1);
        median_val= median( i01_temp( :)); 
        diff_median_mat( j, k)= abs( median_val- center);
        
    end
end

diff_median_mask= diff_median_mat>= ( thres* max( diff_median_mat( :)));
