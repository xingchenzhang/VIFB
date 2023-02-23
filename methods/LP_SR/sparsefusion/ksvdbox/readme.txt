
KSVDBox v12 README
August 24, 2009



KSVDBox installation:
---------------------

1. Make sure OMPBox v9 is installed prior to installing this package.
2. Unpack the contents of the compressed file to a new directory, named e.g. "ksvdbox".
3. If you have not done so before, configure Matlab's MEX compiler by entering
    >> mex -setup
   prior to using MAKE. For optimal performance, it is recommended that you select a compiler
   that performs optimizations. For instance, in Windows, MS Visual Studio is preferred to Lcc.
4. Within Matlab, navigate to the KSVDBox directory, and then to the "private" directory within it,
   and enter MAKE to run the compilation script.
5. Add the KSVDBox package directory to the Matlab path (you can use the ADDPATH command for this).
   Do not add the private directory to the path.


KSVDBox quick start:
--------------------

1. Enter "ksvddemo" and "ksvddenoisedemo" at the Matlab command prompt to run some
   demonstrations of the package.
2. For a complete list of functions in the package, enter
 >> help ksvdbox
This assumes the package was installed to a directory named "ksvdbox". If not, replace ksvdbox
in the above with the (unqualified) name of the KSVDBox installation directory.

Also see faq.txt for some frequently asked questions about the package.
