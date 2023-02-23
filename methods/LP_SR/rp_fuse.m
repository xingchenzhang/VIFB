function Y = rp_fuse(M1, M2, zt, ap, mp)
%Y = fuse_rat(M1, M2, zt, ap, mp) image fusion with ratio pyramid
%
%    M1 - input image A
%    M2 - input image B
%    zt - maximum decomposition level
%    ap - coefficient selection highpass (see selc.m) 
%    mp - coefficient selection base image (see selb.m) 
%
%    Y  - fused image   

%    (Oliver Rockinger 16.08.99)

% check inputs 
[z1 s1] = size(M1);
[z2 s2] = size(M2);
if (z1 ~= z2) | (s1 ~= s2)
  error('Input images are not of same size');
end;

% define filter 
w  = [1 4 6 4 1] / 16;

% define eps
eps = 1e-6;

% cells for selected images
E = cell(1,zt);

% loop over decomposition depth -> analysis
for i1 = 1:zt 
  % calculate and store actual image size 
  [z s]  = size(M1); 
  zl(i1) = z; sl(i1)  = s;
  
  % check if image expansion necessary 
  if (floor(z/2) ~= z/2), ew(1) = 1; else, ew(1) = 0; end;
  if (floor(s/2) ~= s/2), ew(2) = 1; else, ew(2) = 0; end;

  % perform expansion if necessary
  if (any(ew))
  	M1 = adb(M1,ew);
  	M2 = adb(M2,ew);
  end;	
  
  % perform filtering 
  G1 = conv2(conv2(es2(M1,2), w, 'valid'),w', 'valid');
  G2 = conv2(conv2(es2(M2,2), w, 'valid'),w', 'valid');
 
  % decimate, undecimate and interpolate 
  M1T = conv2(conv2(es2(undec2(dec2(G1)), 2), 2*w, 'valid'),2*w', 'valid');
  M2T = conv2(conv2(es2(undec2(dec2(G2)), 2), 2*w, 'valid'),2*w', 'valid');
 
  % select coefficients and store them
  E(i1) = {selc(M1./(M1T+eps), M2./(M2T+eps), ap)};
  
  % decimate 
  M1 = dec2(G1);
  M2 = dec2(G2);
end;

% select base coefficients of last decompostion stage
M1 = selb(M1,M2,mp);

% loop over decomposition depth -> synthesis
for i1 = zt:-1:1
  % undecimate and interpolate 
  M1T = conv2(conv2(es2(undec2(M1), 2), 2*w, 'valid'), 2*w', 'valid');
  % add coefficients
  M1  = (M1T+eps) .* E{i1};
  % select valid image region 
  M1 	= M1(1:zl(i1),1:sl(i1));
end;

% copy image
Y = M1;
