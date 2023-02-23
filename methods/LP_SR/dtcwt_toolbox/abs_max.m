function Y=abs_max(Y1,Y2)
 %functon to get the abs max
 %Y1 and Y2 must have the same dimention
 size_Y=size(Y1);
for i=1:size_Y(1)
    for j=1:size_Y(2)
        if abs(Y1)>=abs(Y2)
            Y=Y1;
        else
            Y=Y2;
        end
    end
end                        %更简洁的方法??