function y = im2gray( x)

% xmin=min(x(:));
% xmax=max(x(:));
% y = (x-xmin)/(xmax-xmin);

x(x<0)=0;
x(x>1)=1;
y=x;