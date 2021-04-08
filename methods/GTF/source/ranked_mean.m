function [ ranked_mean_mask]= ranked_mean( i01, window_side, road_size, thres)

i01_size= size( i01);
ranked_mean_mat= zeros( [ i01_size]);
window_size= ( 2* window_side+ 1)^ 2;
road_mask= [ ones( road_size, 1); zeros( window_size- road_size, 1)];

for j= 1+ window_side: i01_size(1)- window_side                             % rows
    for k= 1+ window_side: i01_size(2)- window_side                         % columns
       
        i01_temp= i01( j- window_side: j+ window_side, k- window_side:...
            k+ window_side);
        center= i01_temp( window_side+ 1, window_side+ 1);
        diff= abs( i01_temp- center);
        diff( 3, 3)= inf;
        sort_v= sort( diff( :));
        sort_v( end)= 0;
        road= sum( road_mask.* sort_v)/ road_size;
        ranked_mean_mat( j, k)= road;
        
    end
end

ranked_mean_mask= ranked_mean_mat> thres* max( ranked_mean_mat( :));
