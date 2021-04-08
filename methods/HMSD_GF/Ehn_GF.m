
function R = Ehn_GF(img)

mi = min(img(:));
ma = max(img(:));
img1 = (img-mi)/(ma-mi)*255+0;

log_img = log(img1+1);

[m,n] = size(img1);
r = floor(0.04*max(m,n));
eps = 0.1;  
% ***using fast guided filter
base = 255*fastguidedfilter_md(img1/255, img1/255, r, eps^2, max(1,r/4)); 
% base = 255*guidedfilter(img1/255, img1/255, r, eps^2);
log_base = log(base+1);

log_detail = log_img - log_base;

target_contrast = log(4);
compressfactor = target_contrast/(max(log_base(:))-min(log_base(:)));
log_absolute_scale = (1-compressfactor)*max(log_base(:));
log_output = compressfactor*log_base + log_detail + log_absolute_scale;
R = exp(log_output);
R = min(R, 255);
end


