% shift_test_2D.m
%
% M-file to perform a 4-level wavelet transform on a circle using Q-shift 
% dual wavelet tree and DWT, and to compare shift invariance properties.
%
% Nick Kingsbury, Cambridge University, May 2002.

clear all
close all

% Draw a circular disc.
x = round((drawcirc(64,1,0,0,256) - 0.5) * 170);
setfig(1); 
colormap(gray(256))
image(min(max(x+128,1),256));
set(gca,'position',[0.1 0.25 .25 .5]);
axis('off');
axis('image');

% draw(xx); 
title('Input (256 x 256)','FontSize',14); drawnow
settitle('Input');
drawnow
% print -depsc circ1.eps

% Do 4 levels of CWT.
[Yl,Yh] = dtwavexfm2(x,4,'near_sym_b','qshift_b');

% Loop to reconstruct output from coefs at each level in turn.
% Starts with the finest level.
titl = ['1st';'2nd';'3rd';'4th';'Low'];

yy = zeros(size(x) .* [2 3]);
yt1 = 1:size(x,1); yt2 = 1:size(x,2);

for mlev = 1:5,
   mask = zeros(6,5);
   mask(:,mlev) = 1;
   z = dtwaveifm2(Yl*mask(1,5),Yh,'near_sym_b','qshift_b',mask);
   figure;draw(z);drawnow
   settitle([titl(mlev,:) ' level DTCWT subbands']);
   yy(yt1,yt2) = z;
   yt2 = yt2 + size(x,2)/2;
end

% disp('Press a key ...')
% pause

% Now do same with DWT.

% Do 4 levels of Real DWT using 'antonini' (9,7)-tap filters.
[Yl,Yh] = wavexfm2(x,4,'antonini');

yt1 = [1:size(x,1)] + size(x,1); yt2 = 1:size(x,2);

for mlev = 1:5,
   mask = zeros(3,5);
   mask(:,mlev) = 1;
   z = waveifm2(Yl*mask(1,5),Yh,'antonini',mask);
   figure;draw(z);drawnow
   settitle([titl(mlev,:) ' level DWT subbands']);
   yy(yt1,yt2) = z;
   yt2 = yt2 + size(x,2)/2;
end

figure;
setfig(gcf); 
colormap(gray(256))
image(min(max(yy+128,1),256));
set(gca,'position',[0.1 0.1 .8 .8]);
axis('off');
axis('image');
hold on
plot(128*[[1;1]*[1:4] [0;6]]+1,128*[[0;4]*[1 1 1 1] [2;2]]+1,'-k');
hold off

title('Components of reconstructed ''disc'' images','FontSize',14);
text(-0.01*size(yy,2),0.25*size(yy,1),'DT CWT','horiz','r');
text(0.02*size(yy,2),1.02*size(yy,1),'wavelets:','horiz','r','vert','t');
text(-0.01*size(yy,2),0.75*size(yy,1),'DWT','horiz','r');
for k=1:4, text(k*128-63,size(yy,1)*1.02,sprintf('level %d',k),'FontSize',14,'horiz','c','vert','t'); end
text(5*128+1,size(yy,1)*1.02,'level 4 scaling fn.','FontSize',14,'horiz','c','vert','t');
drawnow
% print -deps circrecq.eps

disp('Press a key to see perfect reconstruction property ...')
pause

% Accumulate the images from lowband upwards to show perfect reconstruction.
sy = size(x,2)/2;
for mlev = 4:-1:1,
   yt2 = [1:sy] + (mlev-1)*sy;
   yy(:,yt2) = yy(:,yt2) + yy(:,yt2+sy);
end

figure;
setfig(gcf); 
colormap(gray(256))
image(min(max(yy+128,1),256));
set(gca,'position',[0.1 0.1 .8 .8]);
axis('off');
axis('image');

title('Accumulated reconstructions from each level of DT CWT ','FontSize',14);
text(size(yy,2)*0.5,size(yy,1)*1.02,'Accumulated reconstructions from each level of DWT ','FontSize',14,'hor','c','vert','t');
drawnow


return


