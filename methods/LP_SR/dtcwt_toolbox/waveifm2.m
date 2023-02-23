function Z = waveifm2(Yl,Yh,biort,gain_mask);

% Function to perform an n-level DWT 2-D reconstruction.
%
% Z = waveifm2(Yl,Yh,biort,gain_mask);
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
% For example:  Z = waveifm2(Yl,Yh,'near_sym_b');
% performs a reconstruction from Yl,Yh using the 13,19-tap filters. 
%
% Nick Kingsbury, Cambridge University, May 2002

nlevels = length(Yh); % No of levels.
if nargin < 4, gain_mask = ones(3,nlevels); end  % Default gain_mask.

if isstr(biort)		% Check if the biort input is a string
   biort_exist = exist([biort '.mat']);
   if biort_exist == 2,  % Check to see if the filter exists as a .mat file
      load (biort);
   else
      error('Please enter the correct name of the Biorthogonal Filter, see help WAVEIFM2 for details.');
   end
else
   error('Please enter the name of the Biorthogonal Filter as shown in help WAVEIFM2.');
end

LoLo = Yl;
for level = nlevels:-1:1,  % Reconstruct levels in reverse order.
   if size(LoLo,1) ~= size(Yh{level},1)  % If LoLo is not the same height as the next Yh => t1 was extended.
      LoLo = LoLo(2:size(LoLo,1)-1,:);     	% Therefore we have to clip LoLo so it is the same height as the next Yh.
   end
   if size(LoLo,2) ~= size(Yh{level},2)  % If LoLo is not the same width as the next Yh => t1 was extended.
      LoLo = LoLo(:,2:size(LoLo,2)-1);     	% Therefore we have to clip LoLo so it is the same width as the next Yh.
   end
   if any([size(LoLo) 3] ~= size(Yh{level})),
      error('Yh sizes are not valid for WAVEIFM2');
   end
   
   lh = Yh{level}(:,:,1) * gain_mask(1,level);
   hl = Yh{level}(:,:,3) * gain_mask(3,level);
   hh = Yh{level}(:,:,2) * gain_mask(2,level);
   
   % Do even Qshift filters on columns.
   y1 = coliwtfilt(LoLo,2*g0o,0) + coliwtfilt(lh,2*g1o,1);
   y2 = coliwtfilt(hl,2*g0o,0) + coliwtfilt(hh,2*g1o,1);
   % Do even Qshift filters on rows.
   LoLo = (coliwtfilt(y1.',2*g0o,0) + coliwtfilt(y2.',2*g1o,1)).'; 
end
Z = LoLo;

return
