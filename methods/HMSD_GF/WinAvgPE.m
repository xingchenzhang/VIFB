% Compute the per-pixel perceptual energy over each window

function F = WinAvgPE(img, R)

    ff = fft2(img);
    ffc = fftshift(ff);

    [m, n]=size(img);
    [u,v]=freqspace([m,n],'meshgrid');
%     u = u/2*m; v = v/2*n;
%     u = u/m*R; v = v/n*R;
    u = u/2*R; v = v/2*R;
    r=sqrt(u.^2+v.^2);
    fw = MannosSkarision(r);
    
    affc = abs(ffc(:));
    energ = sqrt(sum(affc.*affc));
    F = abs(ffc.*fw);

    F = sqrt(sum(F(:).*F(:)))/(m*n);
end


%%
function fc = MannosSkarision(r)     
    fc = 2.6*(0.0192+0.114*r).*exp(-(0.114*r).^1.1);
end