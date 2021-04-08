function out = Visual_Weight_Map(I)

img=uint8(255*I);
[count, x] = imhist(img);
Sal_Tab = zeros(256,1);
for j=0:255,
    for i=0:255,
    Sal_Tab(j+1) = Sal_Tab(j+1)+count(i+1)*abs(j-i);    
    end      
end
out=zeros(size(img));
for i=0:255,
    out(img==i)=Sal_Tab(i+1);
end 
out=mat2gray(out);

end






