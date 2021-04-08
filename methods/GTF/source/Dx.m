function[u] = Dx(v)

  u = [diff(v,1,2) zeros(size(v,1),1)];
    
return
