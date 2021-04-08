Dual-Tree Complex Wavelet Transform Pack - version 4.3

Nick Kingsbury, Cambridge University, June 2003.

This pack (version 4.2) includes the following functions for performing
the Dual Tree Complex wavelet Transform on 1-D and 2-D signals:

The main DT CWT functions:

  dtwavexfm => 1D DTCWT decomposition
  dtwaveifm => 1D DTCWT reconstruction

  dtwavexfm2 => 2D DTCWT decomposition
  dtwaveifm2 => 2D DTCWT reconstruction

Lower level functions:

  colfilter => Column filtering of a matrix with symmetric extension 
  coldfilt  => Column filtering with decimation by 2, using 2 filters.
  colifilt  => Column filtering with interpolation (upsampling) by 2, using 2 filters.
  reflect   => Reflect a vector about max and min limits (used for sym extension).

There are also equivalent functions for the usual Discrete Wavelet Transform
(included mainly for comparison purposes):

  wavexfm => 1D DWT decomposition
  waveifm => 1D DWT reconstruction

  wavexfm2 => 2D DWT decomposition
  waveifm2 => 2D DWT reconstruction

Lower level DWT functions:

  coldwtfilt  => Column filtering with decimation by 2.
  coliwtfilt  => Column filtering with interpolation (upsampling) by 2.

There are 2 test M-files which show how the main functions 
can be used, and also demonstrate the shift invariance properties
of the DT CWT, compared with the DWT, in both 1-D and 2-D:

  shift_test_1D => Demonstrate shift invariance in 1-D
  shift_test_2D => Demonstrate shift invariance in 1-D

These M-files use the following ancilliary functions:

  draw      => Draw an image in a correctly sized figure window.
  cimage5   => Draw a complex subimage using a colour palette for the complex numbers.
  drawcirc  => Generate a circular disc image.
  setfig    => Set up a predefined figure window.
  settitle  => Set the title of a figure window.

Various .MAT files contain the wavelet filter coefficients.

To test the 1-D DTCWT package use the following Matlab commands:

X = rand(512,1); 
figure; plot(X); drawnow
[Yl,Yh] = dtwavexfm(X,5,'near_sym_b','qshift_b');
Z = dtwaveifm(Yl,Yh,'near_sym_b','qshift_b');
figure; plot(Z)
dtcwt_error = max(abs(Z-X)) % Error should be < 1e-12

To test the 2-D DTCWT package use the following Matlab commands:

load lenna
figure; draw(X); drawnow
[Yl,Yh] = dtwavexfm2(X,4,'near_sym_b','qshift_b');
Z = dtwaveifm2(Yl,Yh,'near_sym_b','qshift_b');
figure; draw(Z)
dtcwt_error = max(abs(Z(:)-X(:))) % Error should be < 1e-12

For more interesting tests, run the M-files shift_test_1D and
shift_test_2D.

********************************

Further information on the DT CWT can be obtained from papers
downloadable from my website (given below). The best tutorial is in
the 1999 Royal Society Paper. In particular this explains the conversion
between 'real' quad-number subimages and pairs of complex subimages. 
The Q-shift filters are explained in the ICIP 2000 paper and in more detail
in the May 2001 paper for the Journal on Applied and Computational 
Harmonic Analysis.

This code is copyright and is supplied free of charge for research
purposes only. In return for supplying the code, all I ask is that, if
you use the algorithms, you give due reference to this work in any
papers that you write and that you let me know if you find any good
applications for the DT CWT. If the applications are good, I would be
very interested in collaboration. I accept no liability arising from use
of these algorithms.

Nick Kingsbury, 
Cambridge University, June 2003.

***********************************************************
Dr N G Kingsbury,
  Dept. of Engineering, University of Cambridge,
    Trumpington St., Cambridge CB2 1PZ, UK.
                      or
    Trinity College, Cambridge CB2 1TQ, UK.
Phone: (0 or +44) 1223 338514 / 332647;  Home: 1954 211152;
Fax: 1223 338564 / 332662;  E-mail: ngk@eng.cam.ac.uk
Web home page: http://www.eng.cam.ac.uk/~ngk/
***********************************************************

Version 4.3: Bug corrected at lines 45 and 48 of waveifm2.m 
(LoLo instead of Lo).  June 2003.



