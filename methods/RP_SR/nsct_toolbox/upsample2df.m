function ho=upsample2df(h, power);
% upsample filter by 2^power;

[m,n]=size(h);
ho   = zeros(2^power * m,2^power * n);
ho(1:2^power:end,1:2^power:end)=h;
