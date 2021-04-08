function [Yl,Yh,Yscale] = wavexfm2(X,nlevels,biort);

% Function to perform a n-level DWT-2D decompostion on a 2-D matrix X.
%
% [Yl,Yh,Yscale] = dtwavexfm2(X,nlevels,biort);
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
% Example: [Yl,Yh] = wavexfm2(X,4,'near_sym_b');
% performs a 4-level 2-D DWT on the real image X using the 13,19-tap filters.
%
% Nick Kingsbury, Cambridge University, May 2002

if isstr(biort)		% Check if the biort input is a string
   biort_exist = exist([biort '.mat']);
   if biort_exist == 2,  % Check to see if the filter exists as a .mat file
      load (biort);
   else
      error('Please enter the correct name of the Biorthogonal Filter, see help WAVEXFM2 for details.');
   end
else
   error('Please enter the name of the Biorthogonal Filter as shown in help WAVEXFM2.');
end

L = size(X);

if any(rem(L,2)),	 % ensure that X is an even length, thus enabling it to be extended if needs be.
   error('Size of X must be a multiple of 2');
end

%initialise
Yh=cell(nlevels,1);
if nargout == 3
   Yscale=cell(nlevels,1);   % This is only required if the user specifies a third output component.
end

LoLo = X;
for level = 1:nlevels;  
   if rem(size(LoLo,1),4),	% Check to see if height of LoLo is divisable by 4, if not extend.
      LoLo = [LoLo(1,:); LoLo; LoLo(end,:)];
   end     
   if rem(size(LoLo,2),4),	% Check to see if height of LoLo is divisable by 4, if not extend.
      LoLo = [LoLo(:,1)  LoLo  LoLo(:,end)];
   end     
   
   % Do filters on rows.
   Lo = coldwtfilt(LoLo,h0o,0).';
   Hi = coldwtfilt(LoLo,h1o,1).';
   
   % Do filters on columns.
   LoLo = coldwtfilt(Lo,h0o,0).';	%LoLo
   Yh{level} = zeros([size(LoLo)  3]);
   Yh{level}(:,:,1) = coldwtfilt(Hi,h0o,0).';	% Horizontal
   Yh{level}(:,:,3) = coldwtfilt(Lo,h1o,1).';	% Vertical
   Yh{level}(:,:,2) = coldwtfilt(Hi,h1o,1).';	% Diagonal   
   if nargout == 3
      Yscale{level} = LoLo;
   end
end   
Yl = LoLo;

return
