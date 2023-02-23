function make
%MAKE Build the KSVDBox MEX support files.
%  MAKE compiles the KSVDBox supporting MEX functions, using Matlab's
%  default MEX compiler. If the MEX compiler has not been set-up before,
%  please run 
%
%    mex -setup
%
%  before using this MAKE file.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% detect platform 

compstr = computer;
is64bit = strcmp(compstr(end-1:end),'64');


% compilation parameters

compile_params = cell(0);
if (is64bit)
  compile_params{1} = '-largeArrayDims';
end


% Compile files %

sourcefiles = {{'addtocols.c'}, {'collincomb.c'}, {'rowlincomb.c'}, {'sprow.c','mexutils.c'}};

for i = 1:length(sourcefiles)
  printf('Compiling %s...', sourcefiles{i}{1});
  mex(sourcefiles{i}{:},compile_params{:});
end
