function output=gaussFilter(I,sigma)

output = I.*0; 
window=double(uint8(3*sigma)*2+1);%窗口大小一半为3*sigma  
H=fspecial('gaussian', window, sigma);%fspecial('gaussian', hsize, sigma)产生滤波模板  
%为了不出现黑边，使用参数'replicate'（输入图像的外部边界通过复制内部边界的值来扩展）

for c=1:size(I,3)
    output(:,:,c)=imfilter(I(:,:,c),H,'replicate');
end

end