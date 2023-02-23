function ompdemo
%OMPDEMO Demonstration of the OMP toolbox.
%  OMPDEMO generates a random sparse mixture of cosines and spikes, adds
%  noise, and applies OMP to recover the original signal.
%
%  To run the demo, type OMPDEMO from the Matlab prompt.
%
%  See also OMPSPEEDTEST.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009


disp(' ');
disp('  **********  OMP Demo  **********');
disp(' ');
disp('  This demo generates a random mixture of cosines and spikes, adds noise,');
disp('  and uses OMP to recover the mixture and the original signal. The graphs');
disp('  show the original, noisy and recovered signal with the corresponding SNR');
disp('  values. The true and recovered coefficients are shown in the bar graphs.');
disp(' ');


% generate DCT-spike dictionary %

n = 256;

D = zeros(n);
D(:,1) = 1/sqrt(n);
for i = 2:n
  v = cos((0:n-1)*pi*(i-1)/n)';
  v = v-mean(v);
  D(:,i) = v/norm(v);
end

D = [D eye(n)];


% generate random sparse mixture of cosines and spikes %

g = sparse(2*n,1);

% sinusoid coefs are random within +/-[0.5,1.5]
lowfreq = 20;
p = randperm(lowfreq);
g(p(1:2)) = (rand(2,1)+0.5).*randsig(2,1);           % two random low frequencies
p = randperm(n-lowfreq);
g(lowfreq+p(1:2)) = (rand(2,1)+0.5).*randsig(2,1);   % two random high frequencies

% spike coefs are random within +/-[0.25,0.75]
p = randperm(n);
g(p(1:3)+n) = (rand(3,1)/2+0.25).*randsig(3,1);      % three random spikes

x = D*g;


% add gaussian noise %

r = randn(size(x));
r = r/norm(r)*norm(x)/4;
y = x + r;


% perform omp %

gamma = omp(D'*y, D'*D, nnz(g));
err = x-D*gamma;


% show results %

v=[1 n -max(abs([x;y]))*1.1 max(abs([x;y]))*1.1];
figure; plot(x); axis(v); title('Original signal');
figure; plot(y); axis(v); title(sprintf('Noisy signal, SNR=%.1fdB', 10*log10((x'*x)/(r'*r))));
figure; plot(D*gamma); axis(v); title(sprintf('Reconstructed signal, SNR=%.1fdB', 10*log10((x'*x)/(err'*err))));

v = [1 2*n -max(abs([g;gamma]))*1.1 max(abs([g;gamma]))*1.1];
figure; bar(full(g)); axis(v); title('True signal decomposition');
figure; bar(full(gamma)); axis(v); title('Decomposition recovered by OMP');

return;


% random matrix with +/-1
function y = randsig(varargin)
y = round(rand(varargin{:}))*2 - 1;
return;
