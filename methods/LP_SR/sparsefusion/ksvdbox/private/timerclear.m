function timerclear()
%TIMERCLEAR Clear all timers.
%   TIMERCLEAR clears all currenly registered timers, invalidating all
%   timer ids.
%
%   Note: since registered timers do not consume CPU power except for when
%   the TIMER<*> functions are called, this function is only useful in
%   situations where a large number of timers have been initialized, and
%   there is a need to reclaim memory.
%
%   See also TIMERINIT, TIMERETA.


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


% clear all timers %

utiltbx_timer_start_times = [];
utiltbx_time_lastdisp = [];
utiltbx_timer_iternums = [];
utiltbx_timer_lastiter = [];
utiltbx_timer_name = [];
utiltbx_timer_callfun = [];
