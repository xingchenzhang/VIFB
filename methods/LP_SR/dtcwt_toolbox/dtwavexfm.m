function [Yl,Yh,Yscale] = dtwavexfm(X,nlevels,biort,qshift);

% Function to perform a n-level DTCWT decompostion on a 1-D column vector X
% (or on the columns of a matrix X).
%
% [Yl,Yh,Yscale] = dtwavexfm(X,nlevels,biort,qshift);
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
%     qshift -> 'qshift_06' => Quarter Sample Shift Orthogonal (Q-Shift) 10,10 tap filters, 
%                              (only 6,6 non-zero taps).
%               'qshift_a' =>  Q-shift 10,10 tap filters,
%                              (with 10,10 non-zero taps, unlike qshift_06).
%               'qshift_b' => Q-Shift 14,14 tap filters.
%               'qshift_c' => Q-Shift 16,16 tap filters.
%               'qshift_d' => Q-Shift 18,18 tap filters.
%               
%
%     Yl     -> The real lowpass subband from the final level.
%     Yh     -> A cell array containing the complex highpass subband for each level.
%     Yscale -> This is an OPTIONAL output argument, that is a cell array containing 
%               the real lowpass coefficients at every scale.
%
% 
% Example: [Yl,Yh] = dtwavexfm(X,5,'near_sym_b','qshift_b');
% performs a 5-level transform on the real image X using the 13,19-tap filters 
% for level 1 and the Q-shift 14-tap filters for levels >= 2.
%
% Nick Kingsbury and Cian Shaffrey
% Cambridge University, May 2002

if isstr(biort) & isstr(qshift)		%Check if the inputs are strings
   biort_exist = exist([biort '.mat']);
   qshift_exist = exist([qshift '.mat']);
   if biort_exist == 2 & qshift_exist == 2;  %Check to see if the filters exist as .mat files
      load (biort);
      load (qshift);
   else
      error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEXFM for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEXFM.');
end

L = size(X);

if any(rem(L(1),2)),	 % ensure that X is an even length, thus enabling it to be extended if needs be.
   error('Size of X must be a multiple of 2');
end

if nlevels == 0, return; end

%initialise
Yh=cell(nlevels,1);
if nargout == 3
   Yscale=cell(nlevels,1);   % This is only required if the user specifies a third output component.
end

j = sqrt(-1);

% Level 1.
Hi = colfilter(X, h1o);   
Lo = colfilter(X, h0o);
t = 1:2:size(Hi,1);
Yh{1} = Hi(t,:) + j*Hi(t+1,:); % Convert Hi to complex form.
if nargout == 3
   Yscale{1} = Lo;
end

if nlevels >= 2;  % Levels 2 and above.
   for level = 2:nlevels;  
      if rem(size(Lo,1),4),	% Check to see if height of Lo is divisable by 4, if not extend.
         Lo = [Lo(1,:); Lo; Lo(end,:)];
      end     
      Hi = coldfilt(Lo,h1b,h1a);
      Lo = coldfilt(Lo,h0b,h0a); 
	   t = 1:2:size(Hi,1);
   	Yh{level} = Hi(t,:) + j*Hi(t+1,:); % Convert Hi to complex form.
   	if nargout == 3
      	Yscale{level} = Lo;
   	end
   end   
end

Yl = Lo;

return
