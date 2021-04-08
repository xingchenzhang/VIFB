function  h11=derotated_dtcwt(n,Yh)
%FUNCTION to perform a derotated DT_CWT

% n  -> No. of levels of wavelet decomposition
% Yh -> A cell array containing the 6 complex highpass subimages for
%  each level.
% h11 ->detotated Yh
h11{n}=Yh{n};
for k=n:-1:2
    for m=1:6
        xp=imresize(Yh{k}(:,:,m),2);%
        argxp=angle(xp);
        argx=angle(Yh{k-1}(:,:,m));
        argx=argx-2.*argxp;
        absx=abs(Yh{k-1}(:,:,m));
        xa=absx.*cos(argx);
        xb=absx.*sin(argx);
        h11{k-1}(:,:,m)=complex(xa,xb);
    end
end
%figure;
%cimage5(h11{1}(:,:,4));
%figure;
%cimage5(h1{1}(:,:,4));
%figure;
%cimage5(h11{2}(:,:,4));
%figure;
%cimage5(h1{2}(:,:,4));