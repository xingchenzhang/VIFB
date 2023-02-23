function [a,b] = fdct_wrapping_pos2idx(NX,NY,SX,SY,s,w,x,y)

% fdct_wrapping_pos2idx.m - For a fixed scale and a fixed direction, returns
%		the curvelet which is closest to a certain point on the image
%
% Inputs
%   NX,NY,SX,SY     values returned by fdct_wrapping_param
%   s               scale index
%   w               wedge (angular) index
%   x,y             position in image
%
% Outputs
%   a,b             Index of the curvelet at scale s and angle w which is nearest to (x,y)
%

  nx = NX{s}{w};  ny = NY{s}{w};
  bx = SX{s}{w}(1,1);  by = SY{s}{w}(1,1);
  hx = SX{s}{w}(2,1)-bx;  hy = SY{s}{w}(1,2)-by;
  a = 1+round((x-bx)/hx);
  b = 1+round((y-by)/hy);
  if(a<1)    a = a+nx;  end
  if(a>nx)    a = a-nx;  end
  if(b<1)    b = b+ny;  end
  if(b>ny)    b = b-ny;  end
  
