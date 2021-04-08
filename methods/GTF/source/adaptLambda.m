function[lambda1 ] = adaptLambda(dI, lambda0, rho, S, val)

[Nrows Ncols depth] = size(dI);

if(depth == 1)
  [Rows, Cols] = find(S>0);
else

  [CellRows{1}, CellCols{1}] = find(S(:,:,1)>0);
  [CellRows{2}, CellCols{2}] = find(S(:,:,2)>0);
  [CellRows{3}, CellCols{3}] = find(S(:,:,3)>0);

end

thresh = ( sum(S(:)>0)/(Nrows*Ncols*depth) )*0.5*rho;
lambda1 = lambda0;

for d = 1:depth

if(depth == 3) 

  Rows = CellRows{d};
  Cols = CellCols{d};

end

N = length(Rows);

for k = 1:N

    w = S(Rows(k), Cols(k), d);

    k_min = Rows(k) - w; if(k_min <= 0) k_min = 1; end
    k_max = Rows(k) + w; if(k_max > Nrows) k_max = Nrows; end
    n_min = Cols(k) - w; if(n_min <= 0) n_min = 1; end
    n_max = Cols(k) + w; if(n_max > Ncols) n_max = Ncols; end
  
    v = dI(k_min:k_max,n_min:n_max,d);
    
    if( norm(v(:),1)/(2*w+1) > thresh ) lambda1(Rows(k), Cols(k), d) = lambda1(Rows(k), Cols(k), d) / val ;
    else lambda1(Rows(k), Cols(k), d) = lambda1(Rows(k), Cols(k), d) * val ;
    end

end

end

return
