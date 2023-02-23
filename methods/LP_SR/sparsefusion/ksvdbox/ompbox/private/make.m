function make
%MAKE Build the OMPBox package.
%  MAKE compiles all OMPBox MEX functions, using Matlab's default MEX
%  compiler. If the MEX compiler has not been set-up before, please run
%
%    mex -setup
%
%  before using this MAKE file.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009


% detect platform 

compstr = computer;
is64bit = strcmp(compstr(end-1:end),'64');


% compilation parameters

compile_params = cell(0);
if (is64bit)
  compile_params{1} = '-largeArrayDims';
end


% Compile files %

ompsources = {'ompcore.c','omputils.c','myblas.c','ompprof.c'};

disp('Compiling ompmex...');
mex('ompmex.c', ompsources{:},compile_params{:});

disp('Compiling omp2mex...');
mex('omp2mex.c',ompsources{:},compile_params{:});

