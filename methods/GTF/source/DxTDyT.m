function[u] = DxTDyT(v, vDims)

  Ncols = vDims(2);
  N     = 2*Ncols;

  if( length(vDims) == 3 )
    V = reshape(v, [vDims(1), N, vDims(3)]);
  else
    V = reshape(v, [vDims(1), N]);
  end


  u = DxT(V(:,1:Ncols)) + DyT(V(:,Ncols+1:N));

%    u = u(:);

return

