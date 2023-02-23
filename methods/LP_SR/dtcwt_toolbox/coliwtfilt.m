function Y = coliwtfilt(X, h, offset)
% function Y = coliwtfilt(X, h, offset)
% Filter the columns of image X using filter vector h, interpolating by 2.
% If length(h) is odd, each input sample is aligned with output samples 
% 1:2:(size(x,1)-1) if offset = 0, and with samples 2:2:size(x,1) if offset = 1.
% If length(h) is even, each input sample is aligned with the mid point
% of each pair of output samples.
% size(Y) = size(X) .* [2 1]
%
% Nick Kingsbury, Cambridge University, May 2002

[r,c] = size(X);
m = length(h);
m2 = fix(m/2);

Y = zeros(r*2,c);

if any(X(:))
   % Adjust symmetric extension and sampling offsets according to 
   % length of h and offset (if length(h) is odd).
   if rem(m,2),  % Odd length h.
      m4 = fix((m+1)/4);
      if offset
         % Symmetrically extend with repeat of lefthand end sample.
         xe = reflect([(1-m4):(r+m4)], 0.5, r);
         if rem(m2,2)
            h1 = h(1:2:m); h2 = h(2:2:m);
            t1 = 1:(r+m2); t2 = 2:(r+m2);
         else
            h1 = h(2:2:m); h2 = h(1:2:m);
            t1 = 1:(r+m2-1); t2 = 1:(r+m2);
         end
      else
         % Symmetrically extend with repeat of righthand end sample.
         xe = reflect([(1-m4):(r+m4)], 1, r+0.5);
         if rem(m2,2)
            h1 = h(2:2:m); h2 = h(1:2:m);
            t1 = 2:(r+m2); t2 = 2:(r+m2+1);
         else
            h1 = h(1:2:m); h2 = h(2:2:m);
            t1 = 1:(r+m2); t2 = 2:(r+m2);
         end
      end
   else          % Even length h.
      m4 = fix(m/4);
      % Symmetrically extend with repeat of end samples.
      xe = reflect([(1-m4):(r+m4)], 0.5, r+0.5);
      if rem(m2,2)
         h1 = h(1:2:m); h2 = h(2:2:m);
         t1 = 1:(r+m2-1); t2 = t1;
      else
         h1 = h(2:2:m); h2 = h(1:2:m);
         t1 = 1:(r+m2-1); t2 = 2:(r+m2);
      end
   end 
   
   % Perform filtering separately on the odd and even samples of the 
   % columns of extended matrix X(xe,:), keeping only the 'valid' output 
   % samples, so Y is half the height of X. 
   s = 1:2:(r*2);
   Y(s,:) = conv2(X(xe(t1),:),h1(:),'valid');
   Y(s+1,:) = conv2(X(xe(t2),:),h2(:),'valid');
end
return;
