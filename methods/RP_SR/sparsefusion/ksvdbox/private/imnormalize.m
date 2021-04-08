function y = imnormalize(x)
%IMNORMALIZE Normalize image values.
%  Y = IMNORMALIZE(X) linearly transforms the intensity values of the image
%  X to tightly cover the range [0,1]. If X has more than one channel, the
%  channels are handled as one and normalized using the same transform.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  May 2004


maxval = max(x(:));
minval = min(x(:));

y = (x-minval) / (maxval-minval);
