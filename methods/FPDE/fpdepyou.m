
function [frth]=fpdepyou(I,T)

[x y z]=size(I);
I=double(I);



dt=0.9;
  

 
I1=I;
I2=I;

k=4;


for  t=1:T
    [Ix,Iy]=gradient(I1); 
    [Ixx,Iyt]=gradient(Ix);
    [Ixt,Iyy]=gradient(Iy);

c=1./(1+(sqrt(Ixx.^2+Iyy.^2)./k).^2);

    [div1,divt1]=gradient(c.*Ixx);
    [divt2,div2]=gradient(c.*Iyy);
    [div11,divt3]=gradient(div1);
    [divt4,div22]=gradient(div2);
    div=div11+div22;
    I2=I1-(dt.*div);
    I1=I2;
end;
frth=uint8(I1);
    