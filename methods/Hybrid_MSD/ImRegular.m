function Iout= ImRegular(Iin)

         m = size(Iin, 1);
         n = size(Iin, 2);
         Imin = min(Iin(:));
         Iout = Iin;
         
         if Imin<0
             a = -Imin;
             for j = 1:m
                 for i = 1:n
                     if Iin(j,i) < a
                         Iout(j,i) = (Iin(j,i)-Imin)^2/(4*a);
                     end
                 end
             end
         end
         
         Imax = max(Iout(:));
         if Imax>255
             b = 510-Imax;
             for j = 1:m
                 for i = 1:n
                     if Iout(j,i) > b
                         Iout(j,i) = (Iout(j,i)-Imax)^2/(4*(b-255)) + 255;
                     end
                 end
             end
         end 
end
