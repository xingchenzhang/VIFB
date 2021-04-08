%This is the Matlab interface to the OMP2 MEX implementation.
%The function is not for independent use, only through omp2.m.


%OMP2MEX Matlab interface to the OMP2 MEX implementation.
%  GAMMA = OMP2MEX(D,X,DtX,XtX,G,EPSILON,SPARSE_G,MSGDELTA,MAXATOMS,PROFILE)
%  invokes the OMP2 MEX function according to the specified parameters. Not
%  all the parameters are required. Those among D, X, DtX, XtX and G which
%  are not specified should be passed as [].
%
%  EPSILON - the target error.
%  SPARSE_G - returns a sparse GAMMA when nonzero, full GAMMA when zero.
%  MSGDELTA - the delay in secs between messages. Zero means no messages.
%  MAXATOMS - the max number of atoms per signal, negative for no max.
%  PROFILE - nonzero means that profiling information should be printed.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  April 2009
