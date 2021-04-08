disp(' ');
disp('fdct_wrapping_wavedemo.m - Curvelets propagating in vacuum');
disp(' ');
disp('A sum of two curvelets is taken as initial condition for the 2-D acoustic');
disp('initial-value problem in a uniform medium with speed of sound c = 1. The');
disp('domain is the periodized unit square [0,1]^2. The interplay between');
disp('pressure and velocity is chosen so that each curvelet is polarized along');
disp('its co-direction and does not split initially. The equations are run for');
disp('2.5 seconds. The computation is made in the Fourier domain.');
disp(' ');
disp('Notice the long-time transverse dispersion of each wave packet (in the');
disp('direction perpendicular to the direction of propagation).');
disp(' ');
disp('By Laurent Demanet and Lexing Ying, 2004')
disp(' ');
disp('An animation should be running in figure(1) ...');

% fdct_wrapping_wavedemo.m - Curvelets propagating in vacuum

N = 256;                    % Grid size, per dimension. Domain is [0,1]^2
h = 1/N;

% Frequencies after fftshift, -pi*N <= wx, wy <= pi*N - 2pi
omega = ((-N/2):(N/2 - 1))*2*pi;
[WX,WY] = meshgrid(omega, omega);
absW = sqrt(WX.^2+WY.^2+(1e-12));

% Initial condition
p0 = zeros(N);
C = fdct_wrapping(p0,0,2,5,32);
C{4}{21}(7,45) = 1;
C{3}{4}(3,3) = 1;
p0 = ifdct_wrapping(C);
clear C;

% Acoustic polarization of each curvelet along its co-direction
p0hat = fftshift(fft2(p0));
u0hat = WX./absW .* p0hat;
v0hat = WY./absW .* p0hat;
fplus0 = (WX./absW .* u0hat + WY./absW .* v0hat + p0hat)/sqrt(2);
fminus0 = (-WX./absW .* u0hat - WY./absW .* v0hat + p0hat)/sqrt(2);

dt = 0.05;
times = 0:dt:2.5;
%M = avifile('fdct_wrapping_demo_wave.avi');
figure(1); clf;
colormap(1-gray);

p0 = real(ifft2(ifftshift(p0hat)));
clow = min(min(p0));
chigh = max(max(p0));
clim = [clow, chigh];

for T = times,
  fplus = fplus0 .* exp(i*T*absW);
  fminus = fminus0 .* exp(-i*T*absW);
  phat = (fplus + fminus)/sqrt(2);
  p = real(ifft2(ifftshift(phat)));
  imagesc(p, clim); axis('image');
  %F = getframe(gca);
  %M = addframe(M,F);
  pause(0.1);
end;
%M = close(M);
