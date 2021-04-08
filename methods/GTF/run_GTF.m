% This is the code for GTF algorithm:
% J. Ma, C. Chen, C. Li, and J. Huang, ¡°Infrared and visible image fusion via gradient transfer 
% and total variation minimization,¡± Information Fusion, vol. 31, pp. 100¨C109, 2016.
%
% The code of GTF is provided by the authors of GTF.
% The interface is created by the authors of VIFB.

function img = run_GTF(imgVI, imgIR, visualization)

    addpath(genpath(cd));

    % IR image
    I = double(imread(imgIR.img))/255;

    % VI image
    V = double(imread(imgVI.img))/255;

    nmpdef;
    pars_irn = irntvInputPars('l1tv');

    pars_irn.adapt_epsR   = 1;
    pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
    pars_irn.adapt_epsF   = 1;
    pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff
    pars_irn.pcgtol_ini = 1e-4;
    pars_irn.loops      = 5;
    pars_irn.U0         = I-V;
    pars_irn.variant       = NMP_TV_SUBSTITUTION;
    pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;
    pars_irn.pcgtol_ini    = 1e-2;
    pars_irn.adaptPCGtol   = 1;

    tic;
    U = irntv(I-V, {}, 4, pars_irn);
    toc;

    X=U+V;
    X=im2uint8(X);
    img = X;

end
