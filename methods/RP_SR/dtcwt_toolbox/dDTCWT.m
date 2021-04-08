
function [Y1 h11]=derotated_dtcwt(I1,n,biot,Qshift)
%     X -> 2D real matrix/Image
%
%     nlevels -> No. of levels of wavelet decomposition
%
%     biort ->  'antonini'   => Antonini 9,7 tap filters.
%               'legall'     => LeGall 5,3 tap filters.
%               'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               'near_sym_b' => Near-Symmetric 13,19 tap filters.
%
%     qshift -> 'qshift_06' => Quarter Sample Shift Orthogonal (Q-Shift) 10,10 tap filters, 
%                              (only 6,6 non-zero taps).
%               'qshift_a' =>  Q-shift 10,10 tap filters,
%                              (with 10,10 non-zero taps, unlike qshift_06).
%               'qshift_b' => Q-Shift 14,14 tap filters.
%               'qshift_c' => Q-Shift 16,16 tap filters.
%               'qshift_d' => Q-Shift 18,18 tap filters.
%               
%
%     Y1     -> The real lowpass image from the final level
%     h11     -> A cell array containing the 6 complex highpass subimages
%     for each level.
%clear all;
%clc;
[Y1,h1] = dtwavexfm2(I1,n,biot,Qshift);
%[Y2,h2] = dtwavexfm2(I2,n,biot,Qshift);

h11{n}=h1{n};
for k=n:-1:2
    for m=1:6
        xp=imresize(h1{k}(:,:,m),2);%
        argxp=angle(xp);
        argx=angle(h1{k-1}(:,:,m));
        argx=argx-2.*argxp;
        absx=abs(h1{k-1}(:,:,m));
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