function F = fuse_MGF(I1, I2, r, eps);
% clc;
% clear all;
% close all;

% parameters 
r1=r; eps1=eps;
r2=r; eps2=eps;
r3=r; eps3=eps;
r4=r; eps4=eps;
            

I1=double(I1);
I2=double(I2); 

%Level 2
A11= guidedfilter(double(I1), double(I2), r1, eps1);
D11=double(I1)-A11;
A21=guidedfilter(double(I2), double(I1), r1, eps1);
D21=double(I2)-A21;
%Level 2
A12= guidedfilter(double(A11), double(A21), r2, eps2);
D12=double(A11)-A12;
A22= guidedfilter(double(A21), double(A11), r2, eps2);
D22=double(A21)-A22;
%Level3
A13= guidedfilter(double(A12), double(A22), r3, eps3);
D13=double(A12)-A13;
A23= guidedfilter(double(A22), double(A12), r3, eps3);
D23=double(A22)-A23;
%Level 4
A14= guidedfilter(double(A13), double(A23), r4, eps4);
D14=double(A13)-A14;
A24= guidedfilter(double(A23), double(A13), r4, eps4);
D24=double(A23)-A24;
%Level 1
S11=abs(double(I1)-double(A11));
S21=abs(double(I2)-double(A21));
%Level 2
S12=abs(double(A11)-double(A12));
S22=abs(double(A21)-double(A22));
%Level 3
S13=abs(double(A12)-double(A13));
S23=abs(double(A22)-double(A23));
% Level 4
S14=abs(double(A13)-double(A14));
S24=abs(double(A23)-double(A24));
%Level 1
w11=S11./(S11+S21);
w21=S21./(S11+S21);
%Level 2
w12=S12./(S12+S22);
w22=S22./(S12+S22);
%Level 3
w13=S13./(S13+S23);
w23=S23./(S13+S23);
%Level 4
w14=S14./(S14+S24);
w24=S24./(S14+S24);

B=1/2*(A14+A24);
 
DL1=D11.*w11+D21.*w21;
DL2=D12.*w12+D22.*w22;
DL3=D13.*w13+D23.*w23;
DL4=D14.*w14+D24.*w24;

FB=B;
FD=DL1+DL2+DL3+DL4;
F=FB+FD;
end
