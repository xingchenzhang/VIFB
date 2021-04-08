%This is the Matlab interface to the OMP MEX implementation.
%The function is not for independent use, only through omp.m.


%OMPMEX Matlab interface to the OMP MEX implementation.
%  GAMMA = OMPMEX(D,X,DtX,G,L,SPARSE_G,MSGDELTA,PROFILE) invokes the OMP
%  MEX function according to the specified parameters. Not all the
%  parameters are required. Those among D, X, DtX and G which are not
%  specified should be passed as [].
%
%  L - the target sparsity.
%  SPARSE_G - returns a sparse GAMMA when nonzero, full GAMMA when zero.
%  MSGDELTA - the delay in secs between messages. Zero means no messages.
%  PROFILE - nonzero means that profiling information should be printed.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009
