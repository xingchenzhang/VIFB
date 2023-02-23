function M=M_interp(I)
[r,c]=size(I);
M=zeros(2*r,2*c);
for i=1:2*r
    for j=1:2*c
        if (mod(i,2)==0)||(mod(j,2)==0)
            M(i,j)=0;
        else 
            M(i,j)=I(floor(i/2+1),floor(j/2+1));
        end
    end
end



