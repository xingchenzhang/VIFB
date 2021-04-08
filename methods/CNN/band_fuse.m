 function Y = band_fuse(M1, M2, W, th)	

% Burts method 
um = 3; 
% compute salience 
S1 = conv2(es2(M1.*M1, floor(um/2)), ones(um), 'valid'); 
S2 = conv2(es2(M2.*M2, floor(um/2)), ones(um), 'valid'); 
% compute match 
MA = conv2(es2(M1.*M2, floor(um/2)), ones(um), 'valid');
MA = 2 * MA ./ (S1 + S2 + eps);

m1 = MA > th; m2 = S1 > S2; 

Y  = (~m1) .* ((m2.*M1) + ((~m2).*M2));
Y  = Y + m1 .* (W.*M1+(1-W).*M2);