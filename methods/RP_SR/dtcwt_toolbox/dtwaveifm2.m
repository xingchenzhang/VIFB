function Z = dtwaveifm2(Yl,Yh,biort,qshift,gain_mask);

% Function to perform an n-level dual-tree complex wavelet (DTCWT)
% 2-D reconstruction.
%
% Z = dtwaveifm2(Yl,Yh,biort,qshift,gain_mask);
%    
%     Yl -> The real lowpass image from the final level
%     Yh -> A cell array containing the 6 complex highpass subimages for each level.
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
%                  gain_mask(d,l) is gain for subband with direction d at level l.
%                  If gain_mask(d,l) == 0, no computation is performed for band (d,l).
%                  Default gain_mask = ones(6,length(Yh)).
%
%     Z -> Reconstructed real image matrix
%
% 
% For example:  Z = dtwaveifm2(Yl,Yh,'near_sym_b','qshift_b');
% performs a 3-level reconstruction from Yl,Yh using the 13,19-tap filters 
% for level 1 and the Q-shift 14-tap filters for levels >= 2.
%
% Nick Kingsbury and Cian Shaffrey
% Cambridge University, May 2002

a = length(Yh); % No of levels.
if nargin < 5, gain_mask = ones(6,a); end  % Default gain_mask.

if isstr(biort) & isstr(qshift)		%Check if the inputs are strings
   biort_exist = exist([biort '.mat']);
   qshift_exist = exist([qshift '.mat']);
   if biort_exist == 2 & qshift_exist == 2; 	%Check to see if the inputs exist as .mat files
      load (biort);
      load (qshift);
   else
      error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEIFM2 for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEIFM2.');
end

current_level = a;
Z = Yl;

while current_level >= 2; ;  %this ensures that for level -1 we never do the following
   lh = c2q(Yh{current_level}(:,:,[1 6]),gain_mask([1 6],current_level));
   hl = c2q(Yh{current_level}(:,:,[3 4]),gain_mask([3 4],current_level));
   hh = c2q(Yh{current_level}(:,:,[2 5]),gain_mask([2 5],current_level));
   
   % Do even Qshift filters on columns.
   y1 = colifilt(Z,g0b,g0a) + colifilt(lh,g1b,g1a);
   y2 = colifilt(hl,g0b,g0a) + colifilt(hh,g1b,g1a);
   % Do even Qshift filters on rows.
   Z = (colifilt(y1.',g0b,g0a) + colifilt(y2.',g1b,g1a)).'; 
   
   % Check size of Z and crop as required
   [row_size col_size] = size(Z);
   S = 2*size(Yh{current_level-1});
   if row_size ~= S(1)		%check to see if this result needs to be cropped for the rows
      Z = Z(2:row_size-1,:);
   end 
   if col_size ~= S(2)		%check to see if this result needs to be cropped for the cols
      Z = Z(:,2:col_size-1);
   end 
   if any(size(Z) ~= S(1:2)),
      error('Sizes of subbands are not valid for DTWAVEIFM2');
   end
   
   current_level = current_level - 1;
end

if current_level == 1;
   
   lh = c2q(Yh{current_level}(:,:,[1 6]),gain_mask([1 6],current_level));
   hl = c2q(Yh{current_level}(:,:,[3 4]),gain_mask([3 4],current_level));
   hh = c2q(Yh{current_level}(:,:,[2 5]),gain_mask([2 5],current_level));

   % Do odd top-level filters on columns.
   y1 = colfilter(Z,g0o) + colfilter(lh,g1o);
   y2 = colfilter(hl,g0o) + colfilter(hh,g1o);
   % Do odd top-level filters on rows.
   Z = (colfilter(y1.',g0o) + colfilter(y2.',g1o)).';
   
end

return

%==========================================================================================
%						**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function x = c2q(w,gain)

% function z = c2q(w,gain)
% Scale by gain and convert from complex w(:,:,1:2) to real quad-numbers in z.
%
% Arrange pixels from the real and imag parts of the 2 subbands
% into 4 separate subimages .
%  A----B     Re   Im of w(:,:,1)
%  |    |
%  |    |
%  C----D     Re   Im of w(:,:,2)

sw = size(w);
x = zeros(2*sw(1:2));

if any(w(:)) & any(gain)
   sc = sqrt(0.5) * gain;
   P = w(:,:,1)*sc(1) + w(:,:,2)*sc(2);
   Q = w(:,:,1)*sc(1) - w(:,:,2)*sc(2);
   
   t1 = 1:2:size(x,1);
   t2 = 1:2:size(x,2);
   
   % Recover each of the 4 corners of the quads.
   x(t1,t2)     = real(P);  % a = (A+C)*sc; 
   x(t1,t2+1)   = imag(P);  % b = (B+D)*sc;
   x(t1+1,t2)   = imag(Q);  % c = (B-D)*sc;
   x(t1+1,t2+1) = -real(Q); % d = (C-A)*sc;
end

return