[x,y,z]=peaks(7); figure(1),  mesh(x,y,z)

[xi,yi]=meshgrid(-3:0.2:3,-3:0.2:3);

z1=interp2(x,y,z,xi,yi,'nearest');

z2=interp2(x,y,z,xi,yi,'linear');

z3=interp2(x,y,z,xi,yi,'spline');

z4=interp2(x,y,z,xi,yi,'cubic');

figure(2),  subplot(2,2,1)

 mesh(xi,yi,z1) 

title('nearest')

subplot(2,2,2) 

mesh(xi,yi,z2)

 title('linear')

subplot(2,2,3)

 mesh(xi,yi,z3) 

title('spiine')

subplot(2,2,4) 

mesh(xi,yi,z4)

 title('cubic')

