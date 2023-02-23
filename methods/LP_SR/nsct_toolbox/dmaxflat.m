function h = dmaxflat(N,d)

% returns 2-D diamond maxflat filters of order 'N' 
% the filters are nonseparable and 'd' is the (0,0) coefficient, being 1 or 0 depending on use
% by Arthur L. da Cunha, University of Illinois Urbana-Champaign
% Aug 2004

if (N > 7 | N < 1), error('N must be in {1,2,3,4,5,6,7}');end;


switch N
 case 1
     h = [0 1 0; 1 0 1; 0 1 0]/4;
     h(2,2) = d;
 case 2
     h = [0 -1 0; -1 0 10; 0 10 0]; 
     h = [h fliplr(h(:,1:end-1))];
     h = [h ; flipud(h(1:end-1,:))]/32;
     h(3,3) = d;
 case 3
     h = [0  3   0  2;
          3  0  -27 0; 
          0 -27  0 174;
          2  0  174  0];
      
     h = [h fliplr(h(:,1:end-1))];
     h = [h ; flipud(h(1:end-1,:))]/512;
     h(4,4) = d;
 case 4
     h = [0   -5     0   -3    0    ;
         -5    0    52    0    34   ; 
          0    52    0   -276  0    ;
         -3    0   -276   0   1454  ;
          0    34    0   1454  0   ]/2^12;
      
     h = [h fliplr(h(:,1:end-1))];
     h = [h ; flipud(h(1:end-1,:))];
     h(5,5) = d;
 case 5
     h = [0    35    0    20    0    18;
          35   0   -425   0   -250    0;
          0   -425   0   2500   0    1610 ;
         20    0    2500  0  -10200   0;
          0   -250    0 -10200  0   47780;
         18    0    1610  0   47780   0]/2^17;
     
     h = [h fliplr(h(:,1:end-1))];
     h = [h ; flipud(h(1:end-1,:))];
     h(6,6) = d;
 case 6
     h = [0    -63    0    -35    0    -30    0   ;
         -63    0    882    0    495    0    444  ;
          0    882    0   -5910   0   -3420   0   ;
         -35    0   -5910   0   25875   0   16460 ;
          0    495    0    25875  0   -89730  0   ;
         -30    0   -3420   0  -89730   0   389112;
          0    44     0    16460  0   389112  0   ]/2^20;
       
     h = [h fliplr(h(:,1:end-1))];
     h = [h ; flipud(h(1:end-1,:))];
     h(7,7) = d;
 case 7
     h = [0    231     0     126     0      105      0      100   ;
         231    0    -3675    0    -2009     0     -1715     0    ;
          0   -3675    0    27930    0     15435     0     13804  ;
         126    0     27930   0   -136514    0    -77910     0    ;
          0   -2009    0   -136514   0     495145    0     311780 ;
         105    0     15435   0    495145    0    -1535709   0    ;
          0   -1715    0   -77910    0    -1535709   0    6305740 ;
         100    0    13804    0    311780    0    6305740    0    ]/2^24 ;
     h = [h fliplr(h(:,1:end-1))];    
     h = [h ; flipud(h(1:end-1,:))];
     h(8,8) = d;
 end
 
     
     