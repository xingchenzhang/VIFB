function varargout = timereta(tid,iter,delay)
%TIMERETA Estimated remaining time.
%   S = TIMERETA(TID,ITER) returns the estimated remaining time (in
%   seconds) for the process associated with timer TID, assuming the
%   process has completed ITER iterations. Note: the function will exit
%   with an error if the timer TID does not have an associated number of
%   iterations (see function TIMERINIT).
%
%   [H,M,S] = TIMERETA(TID,ITER) returns the estimated remaining time in
%   hours, minutes and seconds.
%
%   TIMERETA(TID,ITER), with no output arguments, prints the estimated
%   remaining time to the screen. The time is displayed in the format
%
%     TIMERNAME: iteration ITER / ITERNUM, estimated remaining time: HH:MM:SS.SS
%
%   If the timer has no assigned name, the display format changes to
%
%     Iteration ITER / ITERNUM, estimated remaining time: HH:MM:SS.SS
%
%   TIMERETA(TID,ITER,DELAY) only displays the remaining time if the
%   time elapsed since the previous printout is at least DELAY seconds.
%
%   See also TIMERINIT, TIMERCLEAR.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  June 2008


global utiltbx_timer_start_times 
global utiltbx_timer_iternums
global utiltbx_timer_lastiter
global utiltbx_time_lastdisp
global utiltbx_timer_name


if (tid<1 || tid>length(utiltbx_timer_iternums))
  error('Unknown timer id');
end

if (utiltbx_timer_iternums(tid) < 0)
  error('Specified timer does not have an associated number of iterations');
end

% update last reported iteration number
utiltbx_timer_lastiter(tid) = iter;

% compute elapsed time
starttime = utiltbx_timer_start_times(tid,:);
currtime = clock;
timediff = etime(currtime, starttime);

% total iteration number
iternum = utiltbx_timer_iternums(tid);

% compute eta
timeremain = (iternum-iter)*timediff/iter;

% return eta in seconds
if (nargout==1)
  varargout{1} = timeremain;
  
% return eta in hms
elseif (nargout==3)
  [varargout{1}, varargout{2}, varargout{3}] = secs2hms(timeremain);
  
  
% print eta
elseif (nargout==0)
  
  % check last display time
  lastdisptime = utiltbx_time_lastdisp(tid,:);
  if (nargin>2 && etime(currtime,lastdisptime) < delay)
    return;
  end
  
  % update last display time
  utiltbx_time_lastdisp(tid,:) = currtime;

  % display timer
  [hrs,mins,secs] = secs2hms(timeremain);
  if (isempty(utiltbx_timer_name{tid}))
    printf('Iteration %d / %d, estimated remaining time: %02d:%02d:%05.2f', iter, iternum, hrs, mins, secs);
  else
    timername = utiltbx_timer_name{tid};
    printf('%s: iteration %d / %d, estimated remaining time: %02d:%02d:%05.2f', timername, iter, iternum, hrs, mins, secs);
  end
   
% invalid number of outputs
else
  error('Invalid number of output arguments');
end

