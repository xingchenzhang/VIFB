function Z = dtwaveifm(Yl,Yh,biort,qshift,gain_mask);

% Function to perform an n-level dual-tree complex wavelet (DTCWT)
% 1-D reconstruction.
%
% Z = dtwaveifm(Yl,Yh,biort,qshift,gain_mask);
%    
%     Yl -> The real lowpass subband from the final level
%     Yh -> A cell array containing the complex highpass subband for each level.
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
%     gain_mask -> Gain to be applied to each subband. 
%                  gain_mask(l) is gain for wavelet subband at level l.
%                  If gain_mask(l) == 0, no computation is performed for band (l).
%                  Default gain_mask = ones(1,length(Yh)).
%
%     Z -> Reconstructed real signal vector (or matrix).
%
% 
% For example:  Z = dtwaveifm(Yl,Yh,'near_sym_b','qshift_b');
% performs a reconstruction from Yl,Yh using the 13,19-tap filters 
% for level 1 and the Q-shift 14-tap filters for levels >= 2.
%
% Nick Kingsbury and Cian Shaffrey
% Cambridge University, May 2002

a = length(Yh); % No of levels.
if nargin < 5, gain_mask = ones(1,a); end  % Default gain_mask.

if isstr(biort) & isstr(qshift)		%Check if the inputs are strings
   biort_exist = exist([biort '.mat']);
   qshift_exist = exist([qshift '.mat']);
   if biort_exist == 2 & qshift_exist == 2;        		%Check to see if the inputs exist as .mat files
      load (biort);
      load (qshift);
   else
      error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEIFM for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEIFM.');
end

level = a; 	% No of levels = no of rows in L.

Lo = Yl;
while level >= 2;  % Reconstruct levels 2 and above in reverse order.
   Hi = c2q1d(Yh{level}*gain_mask(level));
   Lo = colifilt(Lo, g0b, g0a) + colifilt(Hi, g1b, g1a);
   
   if size(Lo,1) ~= 2*size(Yh{level-1},1)  % If Lo is not the same length as the next Yh => t1 was extended.
      Lo = Lo(2:size(Lo,1)-1,:);     	% Therefore we have to clip Lo so it is the same height as the next Yh.
   end
   if any(size(Lo) ~= size(Yh{level-1}).*[2 1]),
      error('Yh sizes are not valid for DTWAVEIFM');
   end
   
   level = level - 1;
end

if level == 1;  % Reconstruct level 1.
   Hi = c2q1d(Yh{level}*gain_mask(level));
   Z = colfilter(Lo,g0o) + colfilter(Hi,g1o);
end

return


%==========================================================================================
%				**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function z = c2q1d(x)

% An internal function to convert a 1D Complex vector back to a real array, 
% which is twice the height of x.
[a b] = size(x);
z = zeros(a*2,b);
skip = 1:2:(a*2);
z(skip,:) = real(x);
z(skip+1,:) = imag(x);

return