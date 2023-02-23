function [Yl,Yh,Yscale] = dtwavexfm2(X,nlevels,biort,qshift);

% Function to perform a n-level DTCWT-2D decompostion on a 2D matrix X
%
% [Yl,Yh,Yscale] = dtwavexfm2(X,nlevels,biort,qshift);
%
%     X -> 2D real matrix/Image
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
%     Yl     -> The real lowpass image from the final level
%     Yh     -> A cell array containing the 6 complex highpass subimages for each level.
%     Yscale -> This is an OPTIONAL output argument, that is a cell array containing 
%               real lowpass coefficients for every scale.
%
% 
% Example: [Yl,Yh] = dtwavexfm2(X,3,'near_sym_b','qshift_b');
% performs a 3-level transform on the real image X using the 13,19-tap filters 
% for level 1 and the Q-shift 14-tap filters for levels >= 2.
%
% Nick Kingsbury and Cian Shaffrey
% Cambridge University, Sept 2001


if isstr(biort) & isstr(qshift)		%Check if the inputs are strings
   biort_exist = exist([biort '.mat']);
   qshift_exist = exist([qshift '.mat']);
   if biort_exist == 2 & qshift_exist == 2;        		%Check to see if the inputs exist as .mat files
      load (biort);
      load (qshift);
   else
      error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEXFM2 for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEXFM2.');
end 

orginal_size = size(X);

if ndims(X) >= 3;
   error(sprintf('The entered image is %dx%dx%d, please enter each image slice separately.',orginal_size(1),orginal_size(2),orginal_size(3)));
end

% The next few lines of code check to see if the image is odd in size, if so an extra ...
% row/column will be added to the bottom/right of the image
initial_row_extend = 0;  %initialise
initial_col_extend = 0;
if any(rem(orginal_size(1),2)), %if sx(1) is not divisable by 2 then we need to extend X by adding a row at the bottom
   X = [X; X(end,:)];           %Any further extension will be done in due course.
   initial_row_extend = 1;
end
if any(rem(orginal_size(2),2)), 	%if sx(2) is not divisable by 2 then we need to extend X by adding a col to the left
   X = [X X(:,end)];          %Any further extension will be done in due course.
   initial_col_extend = 1;
end
extended_size = size(X);

if nlevels == 0, return; end

%initialise
Yh=cell(nlevels,1);
if nargout == 3
   Yscale=cell(nlevels,1);   %this is only required if the user specifies a third output component.
end

S = [];
sx = size(X);
if nlevels >= 1,
   
   % Do odd top-level filters on cols.
   Lo = colfilter(X,h0o).';
   Hi = colfilter(X,h1o).';
   
   % Do odd top-level filters on rows.
   LoLo = colfilter(Lo,h0o).';			% LoLo
   Yh{1} = zeros([size(LoLo)/2  6]);
   Yh{1}(:,:,[1 6]) = q2c(colfilter(Hi,h0o).');			% Horizontal pair
   Yh{1}(:,:,[3 4]) = q2c(colfilter(Lo,h1o).');			% Vertical pair
   Yh{1}(:,:,[2 5]) = q2c(colfilter(Hi,h1o).');	      % Diagonal pair
   S = [ size(LoLo) ;S];
   if nargout == 3
      Yscale{1} = LoLo;
   end
end

if nlevels >= 2;
   for level = 2:nlevels;
      [row_size col_size] = size(LoLo);
      if any(rem(row_size,4)),		% Extend by 2 rows if no. of rows of LoLo are divisable by 4;
         LoLo = [LoLo(1,:); LoLo; LoLo(end,:)];
      end 
      if any(rem(col_size,4)),		% Extend by 2 cols if no. of cols of LoLo are divisable by 4;
         LoLo = [LoLo(:,1)  LoLo  LoLo(:,end)];
      end 
      
      % Do even Qshift filters on rows.
      Lo = coldfilt(LoLo,h0b,h0a).';
      Hi = coldfilt(LoLo,h1b,h1a).';
      
      % Do even Qshift filters on columns.
      LoLo = coldfilt(Lo,h0b,h0a).';	%LoLo
      Yh{level} = zeros([size(LoLo)/2  6]);
      Yh{level}(:,:,[1 6]) = q2c(coldfilt(Hi,h0b,h0a).');	% Horizontal
      Yh{level}(:,:,[3 4]) = q2c(coldfilt(Lo,h1b,h1a).');	% Vertical
      Yh{level}(:,:,[2 5]) = q2c(coldfilt(Hi,h1b,h1a).');	% Diagonal   
      S = [ size(LoLo) ;S];
      if nargout == 3
         Yscale{level} = LoLo;
      end
   end
end

Yl = LoLo;

if initial_row_extend == 1 & initial_col_extend == 1;
   warning(sprintf(' \r\r The image entered is now a %dx%d NOT a %dx%d \r The bottom row and rightmost column have been duplicated, prior to decomposition. \r\r ',...
      extended_size(1),extended_size(2),orginal_size(1),orginal_size(2)));
end

if initial_row_extend == 1 ;
   warning(sprintf(' \r\r The image entered is now a %dx%d NOT a %dx%d \r Row number %d has been duplicated, and added to the bottom of the image, prior to decomposition. \r\r',...
      extended_size(1),extended_size(2),orginal_size(1),orginal_size(2),orginal_size(1)));
end

if initial_col_extend == 1;
   warning(sprintf(' \r\r The image entered is now a %dx%d NOT a %dx%d \r Col number %d has been duplicated, and added to the right of the image, prior to decomposition. \r\r',...
      extended_size(1),extended_size(2),orginal_size(1),orginal_size(2),orginal_size(2)));
end
return

%==========================================================================================
%						**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function z = q2c(y)

% function z = q2c(y)
% Convert from quads in y to complex numbers in z.

sy = size(y);
t1 = 1:2:sy(1); t2 = 1:2:sy(2);
j2 = sqrt([0.5 -0.5]);

% Arrange pixels from the corners of the quads into
% 2 subimages of alternate real and imag pixels.
%  a----b
%  |    |
%  |    |
%  c----d

% Combine (a,b) and (d,c) to form two complex subimages. 
p = y(t1,t2)*j2(1) + y(t1,t2+1)*j2(2);     % p = (a + jb) / sqrt(2)
q = y(t1+1,t2+1)*j2(1) - y(t1+1,t2)*j2(2); % q = (d - jc) / sqrt(2)

% Form the 2 subbands in z.
z = cat(3,p-q,p+q);

return
