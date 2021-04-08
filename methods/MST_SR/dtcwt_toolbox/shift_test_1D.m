% SHIFT_TEST_1D.M
% Plot how step responses vary at level 1 to 4 with 16 shifted steps
% for DT CWT and DWT.
%
% Nick Kingsbury, Cambridge University, May 2002.

% Generate 16 step functions as columns of matrix X.
n = 120;
x = cumsum([zeros(n,16);diag(ones(16,1));zeros(n,16)]) - 0.5;

% Specify DT CWT filters here.
biort = 'near_sym_b';
qshift = 'qshift_b';

% Forward DT CWT on each column of x.
[Yl,Yh] = dtwavexfm(x,4,biort,qshift);

% Inverse DT CWT, one subband at a time, using gain_mask to select subbands.
z1 = dtwaveifm(Yl*0,Yh,biort,qshift,[1 0 0 0]);
z01 = dtwaveifm(Yl*0,Yh,biort,qshift,[0 1 0 0]);
z001 = dtwaveifm(Yl*0,Yh,biort,qshift,[0 0 1 0]);
z0001 = dtwaveifm(Yl*0,Yh,biort,qshift,[0 0 0 1]);
z0000 = dtwaveifm(Yl,Yh,biort,qshift,[0 0 0 0]);

% Check for perfect reconstruction: abs(error) < 1e-12.
z = z1 + z01 + z001 + z0001 + z0000;
DTCWT_error = max(abs(z(:)-x(:)))

% Plot 16 input steps, x, and 16 step responses at each level, z1 to z0000.
setfig(1)
set(gcf,'DefaultTextFontSize',12,'Color',[1 1 1]);
set(gcf,'numbertitle','off','name',['Step responses at 16 adjacent shifts']);
subplot('position',[0.15 0.1 0.3 0.8]);
sc = 1.2; % Scale offsets so sets of curves do not overlap.
on = ones(size(x,2)-1,1)/5; 
offset = sc*cumsum([0;on;1;on;1;on;1;on;1;on;1;on]/4).'; % Offsets for each curve.
zdtcwt = [x+0.5 z1 z01 z001 z0001 z0000-0.5] - ones(size(x,1),1)*offset;
tt = -40:40;  % Limits of plot horizontally.
plot(tt,zdtcwt(n+8+tt,:),'-b');
axis off
text(0,-7*sc,'(a) Dual Tree CWT','horiz','c')
xpos = -42; % Position of text labels.
text(xpos,0*sc,'Input','horiz','r','vert','m');
text(xpos,-0.9*sc,'Wavelets','horiz','r','vert','m');
text(xpos,-1.4*sc,'Level 1','horiz','r','vert','m');
text(xpos,-2.4*sc,'Level 2','horiz','r','vert','m');
text(xpos,-3.4*sc,'Level 3','horiz','r','vert','m');
text(xpos,-4.4*sc,'Level 4','horiz','r','vert','m');
text(xpos,-5.5*sc,'Scaling fn','horiz','r','vert','m');
text(xpos,-6*sc,'Level 4','horiz','r','vert','m');
drawnow


% Now do DWT for comparison.

% Specify DWT filters here.
biort = 'antonini';

% Forward DWT on each column of x.
[Yl,Yh] = wavexfm(x,4,biort);

% Inverse DWT, one subband at a time, using gain_mask to select subbands.
z1 = waveifm(Yl*0,Yh,biort,[1 0 0 0]);
z01 = waveifm(Yl*0,Yh,biort,[0 1 0 0]);
z001 = waveifm(Yl*0,Yh,biort,[0 0 1 0]);
z0001 = waveifm(Yl*0,Yh,biort,[0 0 0 1]);
z0000 = waveifm(Yl,Yh,biort,[0 0 0 0]);

% Check for perfect reconstruction: abs(error) < 1e-12.
z = z1 + z01 + z001 + z0001 + z0000;
DWT_error = max(abs(z(:)-x(:)))

subplot('position',[0.47 0.1 0.3 0.8]);
on = ones(size(x,2)-1,1)/5;
offset = sc*cumsum([0;on;1;on;1;on;1;on;1;on;1;on]/4).';
zdwt = [x+0.5 z1 z01 z001 z0001 z0000-0.5] - ones(size(x,1),1)*offset;
tt = -40:40;
plot(tt,zdwt(n+8+tt,:),'-r');
axis off
% xlabel('(b) Real DWT')
text(0,-7*sc,'(b) Real DWT','horiz','c')

return

