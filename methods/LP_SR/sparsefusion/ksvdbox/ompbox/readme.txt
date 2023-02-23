
OMPBox v9 README
May 21, 2009



OMPBox installation:
--------------------

1. Unpack the contents of the compressed file to a new directory, named e.g. "ompbox".
2. If you have not done so before, configure Matlab's MEX compiler by entering
    >> mex -setup
   prior to using MAKE. For optimal performance, it is recommended that you select a compiler
   that performs optimizations. For instance, in Windows, MS Visual Studio is preferred to Lcc.
3. Within Matlab, navigate to the OMPBox directory, and then to the "private" directory within it,
   and enter MAKE to run the compilation script.
4. Add the OMPBox package directory to the Matlab path (you can use the ADDPATH command for this).
   Do not add the private directory to the path.


OMPBox quick start:
-------------------

1. Enter "ompdemo" at the Matlab command prompt to run a short demo of the package.
2. Enter "ompspeedtest" to test the speed of the various OMP implementations.
3. For a complete list of functions in the package, enter
 >> help ompbox
This assumes the package was installed to a directory named "ompbox". If not, replace ompbox
in the above with the (unqualified) name of the OMPBox installation directory.

Also see faq.txt for some frequently asked questions about the package.
