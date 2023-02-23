function str = printf(varargin)
%PRINTF Print formatted text to screen.
%  PRINTF(FMT,VAL1,VAL2,...) formats the data in VAL1,VAL2,... according to
%  the format string FMT, and prints the result to the screen.
%
%  The call to PRINTF(FMT,VAL1,VAL2,...) simply invokes the call
%  DISP(SPRINTF(FMT,VAL1,VAL2,...)). For a complete description of the
%  format string options see function SPRINTF.
%
%  STR = PRINTF(...) also returns the formatted string.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2008


if (nargout>0)
  str = sprintf(varargin{:});
  disp(str);
else
  disp(sprintf(varargin{:}));
end
