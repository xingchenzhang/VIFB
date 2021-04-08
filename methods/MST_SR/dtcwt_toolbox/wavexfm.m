function [Yl,Yh,Yscale] = wavexfm(X,nlevels,biort);

% Function to perform a n-level DWT decompostion on a 1-D column vector X
% (or on the columns of a matrix X).
%
% [Yl,Yh,Yscale] = wavexfm(X,nlevels,biort);
%
%     X -> real 1-D signal column vector (or matrix of vectors)
%
%     nlevels -> No. of levels of wavelet decomposition
%
%     biort ->  'antonini'   => Antonini 9,7 tap filters.
%               'legall'     => LeGall 5,3 tap filters.
%               'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               'near_sym_b' => Near-Symmetric 13,19 tap filters.
%
%     Yl     -> The lowpass subband from the final level.
%     Yh     -> A cell array containing the highpass subband for each level.
%     Yscale -> This is an OPTIONAL output argument, that is a cell array containing 
%               the lowpass coefficients at every scale.
%
% 
% Example: [Yl,Yh] = wavexfm(X,5,'near_sym_b');
% performs a 5-level DWT on the real image X using the 13,19-tap filters.
%
% Nick Kingsbury, Cambridge University, May 2002

if isstr(biort)		% Check if the biort input is a string
   biort_exist = exist([biort '.mat']);
   if biort_exist == 2,  % Check to see if the filter exists as a .mat file
      load (biort);
   else
      error('Please enter the correct name of the Biorthogonal Filter, see help WAVEXFM for details.');
   end
else
   error('Please enter the name of the Biorthogonal Filter as shown in help WAVEXFM.');
end

L = size(X);

if any(rem(L(1),2)),	 % ensure that X is an even length, thus enabling it to be extended if needs be.
   error('Size of X must be a multiple of 2');
end

%initialise
Yh=cell(nlevels,1);
if nargout == 3
   Yscale=cell(nlevels,1);   % This is only required if the user specifies a third output component.
end

Lo = X;
for level = 1:nlevels;  
   if rem(size(Lo,1),4),	% Check to see if height of Lo is divisable by 4, if not extend.
      Lo = [Lo(1,:); Lo; Lo(end,:)];
   end     
   Yh{level} = coldwtfilt(Lo,h1o,1);
   Lo = coldwtfilt(Lo,h0o,0); 
   if nargout == 3
      Yscale{level} = Lo;
   end
end   
Yl = Lo;

return
