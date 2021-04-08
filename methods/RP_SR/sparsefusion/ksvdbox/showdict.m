function x = showdict(D,sz,n,m,varargin)
%SHOWDICT Display a dictionary of image patches.
%  SHOWDICT(D,SZ,N,M) displays the contents of the dictionary D, whos
%  columns are 2-D image patches (in column-major order). SZ = [SX SY] is
%  the size of the image patches. SHOWDICT displays the atoms on an N x M
%  grid. If there are more atoms in D then only the first N*M are
%  displayed.
%
%  SHOWDICT(...,'lines') separates the dictionary atoms by black lines.
%  SHOWDICT(...,'whitelines') separates the dictionary atoms by white
%  lines.
%
%  SHOWDICT(...,'linewidth',W) when used with either 'lines' or
%  'whitelines' sets the width of the lines to W pixels (default=1).
%
%  SHOWDICT(...,'highcontrast') increases the contrast of the figure by
%  normalizing the intensity values of each atom individually to the range
%  of [0,1] (the default behavior is to normalize the values of the entire
%  figure to [0,1] as one image). Note that in this way, the relative
%  intensities of the atoms are not maintained.
%
%  X = SHOWDICT(...) returns a bitmat of the dictionary image without
%  displaying the figure.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009


if (size(D,2) < n*m)
  D = [D zeros(size(D,1),n*m-size(D,2))];
end


%%%  parse input arguments  %%%

linewidth = 1;
highcontrast = 0;
drawlines = 0;
linecolor = 0;

for i = 1:length(varargin)
  if (~ischar(varargin{i}))
    continue;
  end
  switch(varargin{i})
    case 'highcontrast'
      highcontrast = 1;
    case 'lines'
      drawlines = 1;
    case 'whitelines'
      drawlines = 1;
      linecolor = 1;
    case 'linewidth'
      linewidth = varargin{i+1};
  end
end



%%%  create dictionary image  %%%


if (drawlines)
  
  D = [D ; nan(sz(1)*linewidth,size(D,2))];
  sz(2) = sz(2)+linewidth;
  x = col2im(D(:,1:n*m),sz,[n m].*sz,'distinct');
  sz = [sz(2) sz(1)];
  D = im2col(x',sz,'distinct');
  D = [D ; nan(sz(1)*linewidth,size(D,2))];
  sz(2) = sz(2)+linewidth;
  x = col2im(D(:,1:n*m),sz,[m n].*sz,'distinct');
  x = x';
  x = x(1:end-linewidth,1:end-linewidth);
  
  if (highcontrast)
    for i = 0:n-1
      for j = 0:m-1
        x(i*sz(1)+1:i*sz(1)+sz(1)-linewidth, j*sz(2)+1:j*sz(2)+sz(2)-linewidth) = ...
          imnormalize(x(i*sz(1)+1:i*sz(1)+sz(1)-linewidth, j*sz(2)+1:j*sz(2)+sz(2)-linewidth));
      end
    end
  else
    x = imnormalize(x);
  end
  
  x(isnan(x)) = linecolor;
  
else
  
  x = col2im(D(:,1:n*m),sz,[n m].*sz,'distinct');
  
  if (highcontrast)
    for i = 0:n-1
      for j = 0:m-1
        x(i*sz(1)+1:i*sz(1)+sz(1), j*sz(2)+1:j*sz(2)+sz(2)) = ...
          imnormalize(x(i*sz(1)+1:i*sz(1)+sz(1), j*sz(2)+1:j*sz(2)+sz(2)));
      end
    end
  else
    x = imnormalize(x);
  end
end


if (nargout==0)
  imshow(x);
end
