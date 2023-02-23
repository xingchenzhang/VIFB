function v = ompver(history)
%OMPVER OMP toolbox version information.
%  OMPVER displays the current OMP toolbox version information.
%
%  OMPVER('history') also displays history information about the previous
%  versions of the OMP toolbox and their change logs.
%
%  V = OMPVER returns the version number of the current OMP toolbox, and
%  does not display any information.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  May 2009


ver = 9;

if (nargout>0)
  if (nargin>0)
    error('Invalid number of parameters.');
  end
  v = ver;
  
else
  
  if (nargin==0 || (nargin==1 && strcmpi(history,'history')))
    
    disp(' ');
    disp('---------------------------------');
    printf('    OMP Toolbox version %d        ',ver);
    disp('---------------------------------');
    disp(' ');
    
  else
    error('Unknown parameters.');
  end
  
  if (nargin>0)
    
    
    %% changes only displayed since first public release %%
    
    
%     disp(' ');
%     disp('OMP Toolbox version 1, 28.10.07');
%     disp('--------------------------------');
%     disp(' ');
%     disp('Initial release.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%     
%     disp('OMP Toolbox version 2, 30.10.07');
%     disp('--------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. OMP and OMP2 can now operate without specifying G, computing it on the fly.');
%     disp('  2. Improved documentation.');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. OMP2 uses much less memory, by allocating space for at most N coefficients');
%     disp('     for each signal rather than M in v1.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%  
%     disp('OMP Toolbox version 3, 16.12.07');
%     disp('--------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. Support for sparse output, and an optional parameter for selecting');
%     disp('     between sparse and full oupput.');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. When selected atoms were dependent, OMP would cause a NaN result.');
%     disp('     Added a stopping criterion to handle this case.');
%     disp(' ');
%     disp('Internal changes:');
%     disp(' ');
%     disp('  1. Dynamic memory allocations now done using mxMalloc/mxCalloc.');
%     disp('  2. Integer types converted to mwSize/mwIndex to allow compilation');
%     disp('     under 64bit systems using the -largeArrayDims parameter.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%     
%     disp('OMP Toolbox version 4, 03.01.08');
%     disp('--------------------------------');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. Sparse Gamma''s are now represented correctly, with ascending row numbers.');
%     disp(' ');
%     disp('Internal changes:');
%     disp(' ');
%     disp('  1. Added MAKE function for automated compilation of all OMP functions.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%  
%     disp('OMP Toolbox version 5, 18.01.08');
%     disp('--------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. Significantly improved performance by using BLAS/LAPACK functions.');
%     disp('  2. More meaningful profiling reports.');
%     disp(' ');
%     disp('Internal changes:');
%     disp(' ');
%     disp('  1. Changes to code structure: extracted core omp code to separate');
%     disp('     callable functions for code re-use.');
%     disp('  2. Improved MAKE file identifies compilation platform automatically.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%      
%     disp('OMP Toolbox version 6, 31.01.08');
%     disp('--------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. Implemented efficient incremental update of the residual norm in OMP2,');
%     disp('     eliminating the need for explicit computation of the residual itself.');
%     disp('  2. Added ''maxatoms'' parameter to OMP2.');
%     disp('  3. Updated call syntax to support the changes in (1) and (2).');
%     disp('  4. Slightly improved documentation.');
%     disp('  5. Added the OMPVER function.');
%     disp('  6. Added the Contents.m file.');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. Better memory management in OMP2. Much less memory is consumed by using');
%     disp('     dynamic reallocation of temporary helper matrices which were originaly');
%     disp('     preallocated to a size proportional to G, even when G was not specified.');
%     disp('  2. OMP2: atom selection loop limited to N atoms (loop could continue beyond');
%     disp('     this limit if the error target was not achieved, causing a runtime error.');
%     disp('     This could occur when the error limit is very close to 0).');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%     
%     disp('OMP Toolbox version 7, 15.04.08');
%     disp('--------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. OMP2: Improved incremental computation of the residual norm requires less');
%     disp('     memory and is more efficient.');
%     disp('  2. Eliminated separate profiling OMP functions, now unified into the OMP');
%     disp('     and OMP2 functions.');
%     disp('  3. New ''messages'' option for displaying coding progress and remaining time.');
%     disp(' ');
%     disp('Internal changes:');
%     disp(' ');
%     disp('  1. New code structure: unified all omp and omp2 variants into ompmex.c and')
%     disp('     omp2mex.c. Parameter parsing is now done in the Matlab omp.m and omp2.m');
%     disp('     functions, allowing simpler and more advanced parsing.');
%     disp('  2. Much simpler MAKE function, and new create_blasfuncs.m function for');
%     disp('     managing the BLAS function names.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
    
    
    disp(' ');
    disp(' ');
    disp('OMP Toolbox version 7, 15.04.08');
    disp('--------------------------------');
    disp(' ');
    disp('First public release.');
    disp(' ');
    disp(' ');
    disp(' ');


    disp('OMP Toolbox version 8, 29.05.08');
    disp('--------------------------------');
    disp(' ');
    disp('Bug fixes:');
    disp(' ');
    disp('  1. Added OMP stopping criterion: when maximal inner product between the residual');
    disp('     and the dictionary is too small.');
    disp(' ');
    disp(' ');
    disp(' ');
     
    
    disp('OMP Toolbox version 9, 21.05.09');
    disp('--------------------------------');
    disp(' ');
    disp('Features:');
    disp(' ');
    disp('  1. Source code now available for compilation on any Matlab platform.');
    disp('  2. Dictionary normalization is now verified before performing OMP.');
    disp('  3. Added OMPDEMO and OMPSPEEDTEST.');
    disp('  4. Improved profiling information.');
    disp(' ');
    disp('Internal changes:');
    disp(' ');
    disp('  1. Completely restructured code.');
    disp('  2. When G is not specified, OMP-Cholesky is used rather than computing G');
    disp('     and using Batch-OMP (faster).');
    disp('  3. Estimated remaining time now computed with sub-second accuracy.');
    disp(' ');
    disp(' ');
    
  end
end
