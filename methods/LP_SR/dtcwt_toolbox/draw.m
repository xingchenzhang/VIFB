function y=draw(a,map,posn);
% function y=draw(a,map,posn);
% Function to display an image in a box matched to the pixel size.
% If map is not given, a gray-scale of full peak-peak range is used.
% This version tests printon and only uses xpmap if printon = 0.
% If posn is specified, it determines the x and y position of the
% bottom left corner of the figure.

global printon

if nargin < 2, map = []; end

if length(map) > 0,
  graysc = 1;
  mina = 1;
  maxa = max(size(map));
  if maxa > 236,
    map = map(1:236,:);
    a = min(a,236);
  end
  a = round(min(maxa, max(mina, a)));
else
  mina = min(a(:));
  maxa = max(a(:));

% Generate a linear greyscale map.
  map = (([0:63]' + 0.5)/64) * ones(1,3);

  graysc = (length(map) - 1) / (maxa-mina);
  fprintf(1,'Fig %.0f: black = %f, white = %f\n', gcf, mina, maxa);
  a = round((a - mina) * graysc) + 1;
end

sz=get(0,'screensize');

xsize=sz(3);
ysize=sz(4);

[m,n]=size(a);

scale=round(min([xsize/(2*n),(ysize)/(12*m/8)]));
% scale=round(min([512/n,256/m]));  % Better for 2^n image sizes.

figure(gcf);

if scale==0

clf reset
image(a);
colormap(map);
set(gca,'position',[0.01 0.01 .98 .98]);
axis('off');
axis('image');
set(gcf,'position',get(0,'screensize'));

else

clf reset
image(a);
colormap(map);
axis('off');
axis('image')

if nargin < 3,
  pos=[xsize-scale*n+12-gcf*20 rem(gcf-1,10)*20+34 scale*n+4 scale*m+4];
else
  pos=[posn(1)  posn(2)  scale*n+4  scale*m+4];
end
set(gcf,'position',pos);

% scale
u=get(gca,'units');
pos=[2 2 scale*n scale*m];
set(gca,'units','pixels');
set(gca,'position',pos);
set(gca,'units',u);

end



