% The source code is from the Internet
% The interface is modified by the authors of VIFB to integrate it into VIFB. 
%
% Reference for the metric:
% A. M. Eskicioglu and P. S. Fisher, ¡°Image quality measures and their performance,¡± IEEE Transactions on
% communications, vol. 43, no. 12, pp. 2959¨C2965, 1995.

function res = metricsSpatial_frequency(img1, img2, fused)
 
    fused=double(fused);
    [m,n]=size(fused);
    RF=0;
    CF=0;

    for fi=1:m
        for fj=2:n
            RF=RF+(fused(fi,fj)-fused(fi,fj-1)).^2;
        end
    end

    RF=RF/(m*n);

    for fj=1:n
        for fi=2:m
            CF=CF+(fused(fi,fj)-fused(fi-1,fj)).^2;
        end
    end

    CF=CF/(m*n);

    res=sqrt(RF+CF);