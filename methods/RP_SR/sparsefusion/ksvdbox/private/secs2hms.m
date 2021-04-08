function [h,m,s] = secs2hms(t)
%SECS2HMS Convert seconds to hours, minutes and seconds.
%  [H,M,S] = SECS2HMS(T) converts the specified number of seconds T to
%  hours, minutes and seconds. H and M are whole numbers, and S is real.
%
%  Example: Estimate the remaining time of a loop
%
%    n = 10; tic;
%    for i = 1:n
%      pause(1);
%      [h,m,s] = secs2hms( (n-i)*toc/i );
%      printf('estimated remaining time: %02d:%02d:%05.2f',h,m,s);
%    end


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2008


s = t;
h = fix(s/3600);
s = rem(s,3600);
m = fix(s/60);
s = rem(s,60);
