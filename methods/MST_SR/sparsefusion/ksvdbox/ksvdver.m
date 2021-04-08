function v = ksvdver(history)
%KSVDVER K-SVD toolbox version information.
%  KSVDVER displays the current KSVD toolbox version information.
%
%  KSVDVER('history') also displays history information about the previous
%  versions of the KSVD toolbox and their change logs.
%
%  V = KSVDVER returns the version number of the current KSVD toolbox, and
%  does not display any information.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


ver = 12;

if (nargout>0)
  if (nargin>0)
    error('Invalid number of parameters.');
  end
  v = ver;
  
else
  
  if (nargin==0 || (nargin==1 && strcmpi(history,'history')))
    
    disp(' ');
    disp('---------------------------------');
    printf('    KSVD Toolbox version %d       ',ver);
    disp('---------------------------------');
    disp(' ');
    
  else
    error('Unknown parameters.');
  end
  
  if (nargin>0)
    
    
    %% changes only displayed since first public release %%
    
    
%     disp(' ');
%     disp('KSVD Toolbox version 1, 24.10.07');
%     disp('---------------------------------');
%     disp(' ');
%     disp('Initial release.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%     
%     disp('KSVD Toolbox version 2, 13.12.07');
%     disp('---------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. KSVD: support for sparsity and error based coding.');
%     disp('  2. New OMPDENOISE function for generic OMP-based denoising.');
%     disp('  3. ''memusage'' parameter introduced for controlling memory usage.');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. KSVD: when initializing the dictionary with a set of signals, ');
%     disp('     only selects non-zero ones.');
%     disp('  2. KSVD: only negates atoms with first element negative, instead of');
%     disp('     multiplying each atom by its sign (which nullified atoms whose');
%     disp('     first element was zero).');
%     disp('  3. KSVD: tracks training signals used to replace unused atoms, so');
%     disp('     the same one is not used twice.');
%     disp(' ');
%     disp('Internal changes:');
%     disp('  1. Eliminated the use of function cleardict in KSVD, following bugfix (3).');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%  
%     disp('KSVD Toolbox version 3, 16.01.08');
%     disp('---------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. Option for specifying a list of initial signals for the dictionary using');
%     disp('     the parameter initdict.')
%     disp('  2. Changed error computation from Frobenius norm to MSE.');
%     disp('  3. Added the VERBOSE parameter for controlling printed output.');
%     disp('  4. Added the muthresh parameter to KSVD.');
%     disp('  5. Added parameter maxval to OMPDENOISE and KSVDDENOISE.');
%     disp(' ');
%     disp('Internal changes:');
%     disp(' ');
%     disp('  1. Added the sparsecode method for cleaner code.');
%     disp('  2. Re-introduced the use of function cleardict.');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. Function cleardict now checks correlation between atoms in absolute value.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%     
%     disp('KSVD Toolbox version 4, 27.01.08');
%     disp('---------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. Added functions OMPDENOISE2 and OMPDENOISE3 for faster denoising of');
%     disp('     2-D and 3-D signals.');
%     disp('  2. Added gain and lambda parameters to OMPDENOISE and KSVDDENOISE.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%  
%     disp('KSVD Toolbox version 5, 03.02.08');
%     disp('---------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. KSVD: atom update is now in random order rather than sequential.');
%     disp('  2. KSVD: improved error computation: error is now computed after the');
%     disp('     dictionary update step, eliminating the need for an additional sparse');
%     disp('     coding step at the end of the algorithm for error computation.');
%     disp('     Computed error is now RMSE for sparsity-based minimization, and average');
%     disp('     atom count for error-based minimization.');
%     disp('  3. KSVD: cleardict now replaces every atom not used by at least 4 training');
%     disp('     signals. Atoms already replaced in the dictionary update step are');
%     disp('     excluded from this criterion.');
%     disp('  4. Added the ''exact'' parameter to KSVD/KSVDDENOISE, to allow using the exact');
%     disp('     SVD solution rather than the default quicker approximation method.');
%     disp('  5. Changed behavior of VERBOSE=3 mode in KSVD, and added VERBOSE=4 mode.');
%     disp('  6. VERBOSE in OMPDENOISE/2/3 now prints a message every VERBOSE seconds,');
%     disp('     not every VERBOSE blocks.');
%     disp('  7. OMPDENOISE: faster computation of the matrix cnt (a counter matrix');
%     disp('     containing the number of blocks covering each sample in the signal.');
%     disp('  8. KSVDDENOISE: new default values for gain and lambda.');
%     disp('  9. KSVDDENOISE: reverts to function OMPDENOISE when memusage is ''low''');
%     disp(' 10. Set a limit on the number of atoms used to code each block in OMPDENOISE');
%     disp('     and KSVDDENOISE. Added the maxatoms parameter to control this limit.');
%     disp(' 11. Improved documentation.');
%     disp(' 12. Added the KSVDVER function.');
%     disp(' ');
%     disp('Internal changes:');
%     disp(' ');
%     disp('  1. Support for the new call syntax of OMP Toolbox v6.');
%     disp('  2. OMPDENOISE3 now computes the matrix cnt using much less memory. Instead');
%     disp('     of allocating three matrices the size of the input signal, it only');
%     disp('     allocates one such matrix.');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. OMPDENOISE2/3 now compute the matrix cnt correctly.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%     
%     disp('KSVD Toolbox version 6, 14.02.08');
%     disp('---------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. New stepsize parameter in OMPDENOISE/2/3 and KSVDDENOISE.');
%     disp('  2. ''blocksize'' parameter in OMPDENOISE/2/3 and KSVDDENOISE can now be a scalar.');
%     disp('  3. Added the Contents.m file.');
%     disp('  4. Revised VERBOSE settings in KSVDDENOISE, and slightly improved display format');
%     disp('     in KSVD.');
%     disp(' ');
%     disp('Bug fixes:');
%     disp(' ');
%     disp('  1. OMPDENOISE/2/3: fixed a bug in the computation of the matrix cnt when');
%     disp('     the signal is not square.');
%     disp(' ');
%     disp(' ');
%     disp(' ');
%     
%     disp('KSVD Toolbox version 7, 18.04.08');
%     disp('---------------------------------');
%     disp(' ');
%     disp('Features:');
%     disp(' ');
%     disp('  1. OMPDENOISE now prints the estimated remaining time.');
%     disp('  2. New ''msgdelta'' parameter in KSVD and KSVDDENOISE for printing');
%     disp('     periodic status updates.');
%     disp(' ');
%     disp('Internal changes:');
%     disp(' ');
%     disp('  1. Compatible with OMP Toolbox v7.');
%     disp(' ');
%     disp(' ');
%     disp(' ');


    disp(' ');
    disp(' ');
    disp('KSVD Toolbox version 7, 18.04.08');
    disp('---------------------------------');
    disp(' ');
    disp('First public release.');
    disp(' ');
    disp(' ');
    disp(' ');
    
    
    disp('KSVD Toolbox version 8, 29.05.08');
    disp('--------------------------------');
    disp(' ');
    disp('Features:');
    disp(' ');
    disp('  1. Added ''sigma'' and ''noisemode'' parameters to OMPDENOISE and KSVDDENOISE.');
    disp('  2. Added the overcomplete DCT dictionary as an optional initial dictionary');
    disp('     in KSVDDENOISE, by setting the initdict parameter to the string ''odct''.');
    disp('  3. Improved default gain and lambda parameters in OMPDENOISE.');
    disp(' ');
    disp('Bug fixes:');
    disp(' ');
    disp('  1. KSVDDENOISE: corrected a bug which caused incorrect gain and lambda');
    disp('     parameters to be sent to OMPDENOISE in some cases.');
    disp(' ');
    disp('Internal changes:');
    disp(' ');
    disp('  1. When replacing an unused atom, the search for an example with maximal');
    disp('     error is limited to a fixed number of signals.');
    disp('  2. Improved implementation of KSVDVER.');
    disp('  3. Uses the TIMER package for time estimation.');
    disp(' ');
    disp(' ');
    disp(' ');
    
        
    disp('KSVD Toolbox version 9, 30.11.08');
    disp('--------------------------------');
    disp(' ');
    disp('Features:');
    disp(' ');
    disp('  1. OMPDENOISE and KSVDDENOISE can now return the average number of');
    disp('     non-zeros in the denoised block representations.');
    disp(' ');
    disp('Bug fixes:');
    disp(' ');
    disp('  1. Corrected an error when parameter MSGDELTA was not specified.');
    disp(' ');
    disp('Internal changes:');
    disp(' ');
    disp('  1. The COUNTCOVER function counts the number of times a pixel in a');
    disp('     multi-dimensional signal is covered by overlapping blocks.');
    disp('  2. Function CLEARDICT now computes the error norm in blocks to reduce');
    disp('     memory consumption.');
    disp(' ');
    disp(' ');
    disp(' ');
      
        
    disp('KSVD Toolbox version 10, 21.05.09');
    disp('---------------------------------');
    disp(' ');
    disp('Features:');
    disp(' ');
    disp('  1. Source code for supporting C functions now available for compilation');
    disp('     on any Matlab platform.');
    disp('  2. Added KSVDDEMO and KSVDDENOISEDEMO.');
    disp('  3. Improved performance of OMPDENOISE2, OMPDENOISE3 and KSVD.');
    disp('  4. New SHOWDICT function for displaying a dictionary of image patches.');
    disp('  5. KSVD: improved VERBOSE/MSGDELTA parameter syntax.');
    disp('  6. KSVD: messages are now displayed by default.');
    disp('  7. KSVD: only one of the parameters dictsize/initdict needs to be specified.');
    disp('  8. KSVD: eliminated the data saving option.');
    disp('  9. OMPDENOISE: messages are now displayed by default.');
    disp(' 10. KSVDDENOISE: simplified MSGDELTA parameter.');
    disp(' 11. KSVDDENOISE: default initial dictionary is now overcomplete DCT.');
    disp(' 12. KSVDDENOISE: parameter initdict can be the string ''data'' to indicate');
    disp('     initialization from training examples.');
    disp(' ');
    disp(' ');
    disp(' ');
    
    
    disp('KSVD Toolbox version 11, 03.08.09');
    disp('---------------------------------');
    disp(' ');
    disp('Features:');
    disp(' ');
    disp('  1. Added OMPDENOISE1 for fast OMP denoising of 1-D signals.');
    disp('  2. Added support for 1-D signals in KSVDDENOISE.');
    disp(' ');
    disp('Bug fixes:');
    disp('  1. Corrected a bug in OMPDENOISE/2/3 when memusage is low.');
    disp(' ');
    disp('Internal changes:');
    disp('  1. KSVD: only compute sum(data.*data) when memusage is high.');
    disp(' ');
    disp(' ');
    disp(' ');
    
    
    disp('KSVD Toolbox version 12, 24.08.09');
    disp('---------------------------------');
    disp(' ');
    disp('Features:');
    disp(' ');
    disp('  1. KSVD: significantly faster atom update computation via optimized');
    disp('     sparse-matrix code.');
    disp('  2. KSVD: reduced memory consumption by computing one column of G');
    disp('     at a time when clearing the dictionary, and using block processing');
    disp('     to compute representation error.');
    disp(' ');
    disp('Internal changes:');
    disp('  1. KSVD: when replacing an atom in cleardict, the search for an example');
    disp('     with maximal error is limited to a fixed number of signals.');
    disp('  2. KSVDDENOISEDEMO: memusage is set to ''high'' for faster performance.');
    disp(' ');
    disp(' ');
    
   end
end
