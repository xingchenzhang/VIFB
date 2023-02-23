%函数sobel.m: sobel算子%
function y=sobel(image,th)  
%image: 输入图像%
%th: 临界值%
[Nn Nm]=size(image);           %图片大小
h=[-1 -2 -1;0 0 0;1 2 1];      %Sobel 算子
Gx=filter2(h,image);           %得梯度向量的Gx分量  
Gy=filter2(h',image);          %得梯度向量的Gy分量
F=abs(Gx)+abs(Gy);        
Bw_Gx=double(Gx)/255;          %将Gx范围从[0, 255]缩到[0,1]  
Bw_Gy=double(Gy)/255;          %将Gy范围从[0, 255]缩到[0,1]  
Bw_F=double(F)/255;            %将F范围从[0, 255]缩到[0,1]
                               %设定临界值，在[0,1]之间 
for i=1:Nn
  for j=1:Nm 
    if Bw_F(i,j)<th     
      y(i,j)=0;
    else
      y(i,j)=1;
    end
  end
end
