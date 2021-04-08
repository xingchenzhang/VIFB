function[u] = DyT(v)

  u0 = -v(1,:);
  u1 = -diff(v);
  u2 = v(end-1,:);
  u = [u0; u1(1:(end-1),:); u2];
    
return
