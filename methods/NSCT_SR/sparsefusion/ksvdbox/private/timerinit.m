function tid = timerinit(par1,par2)
%TIMERINIT Initialize a new timer.
%   TID = TIMERINIT() initializes a new timer for counting elapsed time,
%   and returns its id.
%
%   TID = TIMERINIT('TIMERNAME') sets the timer name to the specified
%   string for display purposes.
%
%   TID = TIMERINIT(ITERNUM) initializes a new ETA timer for a process with
%   ITERNUM iterations. An ETA timer can be used for both counting elapsed
%   time and estimating remaining time.
%
%   TID = TIMERINIT('TIMERNAME',ITERNUM) sets the ETA timer name to the
%   specified string for display purposes.
%
%   Example:
%
%     tid = timerinit(100); 
%     for i = 1:100
%       pause(0.07);
%       timereta(tid,i,1);
%     end
%     timereta(tid,i);
%
%   See also TIMERETA, TIMERCLEAR.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  June 2008


global utiltbx_timer_start_times    % start times
global utiltbx_time_lastdisp        % last display times
global utiltbx_timer_iternums       % iteration numbers
global utiltbx_timer_lastiter       % last queried iteration numbers
global utiltbx_timer_name           % timer names
global utiltbx_timer_callfun        % timer calling functions


% parse function arguments %

if (nargin==0)

  iternum = -1;
  timername = '';

elseif (nargin==1)

  if (ischar(par1))
    iternum = -1;
    timername = par1;

  elseif (isnumeric(par1) && numel(par1)==1 && par1>0)
    iternum = par1;
    timername = '';

  else
    error('Invalid number of iterations');
  end

elseif (nargin==2)

  if (ischar(par1) && isnumeric(par2))
    if (numel(par2)==1 && par2>0)
      timername = par1;
      iternum = par2;
    else
      error('Invalid number of iterations');
    end
  else
    error('Invalid function syntax');
  end

else
  error('Too many function parameters');
end


% register the timer %

if (isempty(utiltbx_timer_start_times))
  utiltbx_timer_start_times = clock;
  utiltbx_time_lastdisp = utiltbx_timer_start_times;
  utiltbx_timer_iternums = double(iternum);
  utiltbx_timer_lastiter = 0;
  utiltbx_timer_name = { timername };
  utiltbx_timer_callfun = {};
  tid = 1;
else
  utiltbx_timer_start_times(end+1,:) = clock;
  utiltbx_time_lastdisp(end+1,:) = utiltbx_timer_start_times(end,:);
  utiltbx_timer_iternums(end+1) = double(iternum);
  utiltbx_timer_lastiter(end+1) = 0;
  utiltbx_timer_name{end+1} = timername;
  tid = size(utiltbx_timer_start_times,1);
end


% detect timer calling function %

st = dbstack;
if (length(dbstack) >= 2)
  utiltbx_timer_callfun{end+1} = st(2).name;
else
  utiltbx_timer_callfun{end+1} = '';
end
