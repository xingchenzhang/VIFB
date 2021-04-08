function Z = waveifm(Yl,Yh,biort,gain_mask);

% Function to perform an n-level dual-tree complex wavelet (DTCWT)
% 1-D reconstruction.
%
% Z = waveifm(Yl,Yh,biort,gain_mask);
%    
%     Yl -> The real lowpass subband from the final level
%     Yh -> A cell array containing the complex highpass subband for each level.
%
%     biort ->  'antonini'   => Antonini 9,7 tap filters.
%               'legall'     => LeGall 5,3 tap filters.
%               'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               'near_sym_b' => Near-Symmetric 13,19 tap filters.
%
%     gain_mask -> Gain to be applied to each subband. 
%                  gain_mask(l) is gain for wavelet subband at level l.
%                  If gain_mask(l) == 0, no computation is performed for band (l).
%                  Default gain_mask = ones(1,length(Yh)).
%
%     Z -> Reconstructed real signal vector (or matrix).
%
% 
% For example:  Z = waveifm(Yl,Yh,'near_sym_b');
% performs a reconstruction from Yl,Yh using the 13,19-tap filters. 
%
% Nick Kingsbury, Cambridge University, May 2002

nlevels = length(Yh); % No of levels.
if nargin < 4, gain_mask = ones(1,nlevels); end  % Default gain_mask.

if isstr(biort)		% Check if the biort input is a string
   biort_exist = exist([biort '.mat']);
   if biort_exist == 2,  % Check to see if the filter exists as a .mat file
      load (biort);
   else
      error('Please enter the correct name of the Biorthogonal Filter, see help WAVEIFM for details.');
   end
else
   error('Please enter the name of the Biorthogonal Filter as shown in help WAVEIFM.');
end

Lo = Yl;
for level = nlevels:-1:1,  % Reconstruct levels in reverse order.
   if size(Lo,1) ~= size(Yh{level},1)  % If Lo is not the same length as the next Yh => t1 was extended.
      Lo = Lo(2:size(Lo,1)-1,:);     	% Therefore we have to clip Lo so it is the same height as the next Yh.
   end
   if any(size(Lo) ~= size(Yh{level})),
      error('Yh sizes are not valid for WAVEIFM');
   end
   Hi = Yh{level}*gain_mask(level);
   Lo = coliwtfilt(Lo,2*g0o,0) + coliwtfilt(Hi,2*g1o,1);
end
Z = Lo;

return
