function im_3c = make_3c(I)
[m,n] = size(I);
im_temp = zeros(m,n,3);

im_temp(:,:,1) = I;
im_temp(:,:,2) = I;
im_temp(:,:,3) = I;

im_3c = im_temp;
end