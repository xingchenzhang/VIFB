function [h0, h1] = dfilters(fname, type)
% DFILTERS	Generate directional 2D filters
%
%	[h0, h1] = dfilters(fname, type)
%
% Input:
%	fname:	Filter name.  Available 'fname' are:
%		'haar':		the "Haar" filters
%		'vk':		McClellan transformed of the filter from
%				    the VK book
%		'ko':		orthogonal filter in the Kovacevic's paper
%		'kos':		smooth 'ko' filter
%		'lax':		17 x 17 by Lu, Antoniou and Xu
%		'sk':		9 x 9 by Shah and Kalker
%		'cd':		7 and 9 McClellan transformed by 
%				                       Cohen and Daubechies
%		'pkva':		ladder filters by Phong et al.
%		'oqf_362':	regular 3 x 6 filter
%       'dvmlp'  :  regular linear phase biorthogonal filter with 3 dvm
%		'sinc':		ideal filter (*NO perfect recontruction*)
%       'dmaxflat': diamond maxflat filters obtained from a three stage ladder
%
%	     type:	'd' or 'r' for decomposition or reconstruction filters
%
% Output:
%	h0, h1:	diamond filter pair (lowpass and highpass)

% To test those filters (for the PR condition for the FIR case), verify that:
% conv2(h0, modulate2(h1, 'b')) + conv2(modulate2(h0, 'b'), h1) = 2
% (replace + with - for even size filters)
%
% To test for orthogonal filter
% conv2(h, reverse2(h)) + modulate2(conv2(h, reverse2(h)), 'b') = 2

% The diamond-shaped filter pair
switch fname
    case {'haar'}
	if lower(type(1)) == 'd'
	    h0 = [1, 1] / sqrt(2);
	    h1 = [-1, 1] / sqrt(2);
	else
	    h0 = [1, 1] / sqrt(2);
	    h1 = [1, -1] / sqrt(2);
	end
	
    case 'vk'	% in Vetterli and Kovacevic book
	if lower(type(1)) == 'd'
	    h0 = [1, 2, 1] / 4;
	    h1 = [-1, -2, 6, -2, -1] / 4;
	else
	    h0 = [-1, 2, 6, 2, -1] / 4;
	    h1 = [-1, 2, -1] / 4;
	end
	
	% McClellan transfrom to obtain 2D diamond filters
	t = [0, 1, 0; 1, 0, 1; 0, 1, 0] / 4;	% diamond kernel
	h0 = ftrans2(h0, t);
	h1 = ftrans2(h1, t);
	
    case 'ko'	% orthogonal filters in Kovacevic's thesis
	a0 = 2; a1 = 0.5; a2 = 1;
	
	h0 = [0,      -a1, -a0*a1, 0;
	      -a2, -a0*a2, -a0,    1;
	      0, a0*a1*a2, -a1*a2, 0];
	
	% h1 = qmf2(h0);
	h1 = [0, -a1*a2, -a0*a1*a2, 0;
	      1,     a0, -a0*a2, a2;
	      0, -a0*a1,  a1, 0];
	
	% Normalize filter sum and norm;
	norm = sqrt(2) / sum(h0(:));
	
	h0 = h0 * norm;	
	h1 = h1 * norm;
	      
	if type == 'r'
	    % Reverse filters for reconstruction
	    h0 = h0(end:-1:1, end:-1:1);
	    h1 = h1(end:-1:1, end:-1:1);
	end
    
    case 'kos'	% Smooth orthogonal filters in Kovacevic's thesis
	a0 = -sqrt(3); a1 = -sqrt(3); a2 = 2+sqrt(3);
	
	h0 = [0,      -a1, -a0*a1, 0;
	      -a2, -a0*a2, -a0,    1;
	      0, a0*a1*a2, -a1*a2, 0];
	
	% h1 = qmf2(h0);
	h1 = [0, -a1*a2, -a0*a1*a2, 0;
	      1,     a0, -a0*a2, a2;
	      0, -a0*a1,  a1, 0];
	
	% Normalize filter sum and norm;
	norm = sqrt(2) / sum(h0(:));
	
	h0 = h0 * norm;	
	h1 = h1 * norm;
	      
	if type == 'r'
	    % Reverse filters for reconstruction
	    h0 = h0(end:-1:1, end:-1:1);
	    h1 = h1(end:-1:1, end:-1:1);
	end
	
    case 'lax'	% by Lu, Antoniou and Xu
	h = [-1.2972901e-5,  1.2316237e-4, -7.5212207e-5,  6.3686104e-5, ...
	      9.4800610e-5, -7.5862919e-5,  2.9586164e-4, -1.8430337e-4; ...
	    
	      1.2355540e-4, -1.2780882e-4, -1.9663685e-5, -4.5956538e-5, ...
	     -6.5195193e-4, -2.4722942e-4, -2.1538331e-5, -7.0882131e-4; ...
	    
	     -7.5319075e-5, -1.9350810e-5, -7.1947086e-4,  1.2295412e-3, ...
	      5.7411214e-4,  4.4705422e-4,  1.9623554e-3,  3.3596717e-4; ...
	    
	      6.3400249e-5, -2.4947178e-4,  4.4905711e-4, -4.1053629e-3, ...
	     -2.8588307e-3,  4.3782726e-3, -3.1690509e-3, -3.4371484e-3; ...
	    
	      9.6404973e-5, -4.6116254e-5,  1.2371871e-3, -1.1675575e-2, ...
	      1.6173911e-2, -4.1197559e-3,  4.4911165e-3,  1.1635130e-2; ...
	    
	     -7.6955555e-5, -6.5618379e-4,  5.7752252e-4,  1.6211426e-2, ...
	      2.1310378e-2, -2.8712621e-3, -4.8422645e-2, -5.9246338e-3; ...
	    
	      2.9802986e-4, -2.1365364e-5,  1.9701350e-3,  4.5047673e-3, ...
	     -4.8489158e-2, -3.1809526e-3, -2.9406153e-2,  1.8993868e-1; ...
	    
	     -1.8556637e-4, -7.1279432e-4,  3.3839195e-4,  1.1662001e-2, ...
	     -5.9398223e-3, -3.4467920e-3,  1.9006499e-1,  5.7235228e-1];
	
	h0 = sqrt(2) * [h, h(:, end-1:-1:1); ...
			h(end-1:-1:1, :), h(end-1:-1:1, end-1:-1:1)];	
	
	h1 = modulate2(h0, 'b');
	
    case 'sk'	% by Shah and Kalker
	h = [ 0.621729,    0.161889,  -0.0126949, -0.00542504, 0.00124838; ...
	      0.161889,   -0.0353769, -0.0162751, -0.00499353, 0; ...
	     -0.0126949,  -0.0162751,  0.00749029, 0, 0; ...
	     -0.00542504,  0.00499353, 0, 0, 0; ...
	      0.00124838, 0, 0, 0, 0];
	
	h0 = sqrt(2) * [h(end:-1:2, end:-1:2), h(end:-1:2, :); ...
			h(:, end:-1:2), h];
	
	h1 = modulate2(h0, 'b');
	
    case 'dvmlp'
    q= sqrt(2); b=.02; b1 = b*b;
    h  = [b/q 0 -2*q*b 0 3*q*b 0 -2*q*b 0 b/q;
          0 -1/(16*q) 0 9/(16*q) 1/q  9/(16*q) 0 -1/(16*q) 0;
          b/q 0 -2*q*b 0 3*q*b 0 -2*q*b 0 b/q];     
    g0 = [-b1/q   0   4*b1*q   0    -14*q*b1    0   28*q*b1    0       -35*q*b1  0    28*q*b1    0     -14*q*b1   0    4*b1*q  0     -b1/q;
           0     b/(8*q)   0  -13*b/(8*q) b/q  33*b/(8*q) -2*q*b -21*b/(8*q) 3*q*b -21*b/(8*q)   -2*q*b 33*b/(8*q)  b/q -13*b/(8*q) 0  b/(8*q) 0;
          -q*b1  0  -1/(256*q) + 8*q*b1  0   9/(128*q) - 28*q*b1   -1/(q*16)   -63/(256*q) + 56*q*b1   9/(16*q)   87/(64*q)-70*q*b1      9/(16*q)   -63/(256*q) + 56*q*b1    -1/(q*16)   9/(128*q) - 28*q*b1   0   -1/(256*q) + 8*q*b1   0   -q*b1;
          0     b/(8*q)   0  -13*b/(8*q) b/q  33*b/(8*q) -2*q*b -21*b/(8*q) 3*q*b -21*b/(8*q)   -2*q*b 33*b/(8*q)  b/q -13*b/(8*q) 0  b/(8*q) 0;
          -b1/q   0   4*b1*q   0    -14*q*b1    0   28*q*b1    0       -35*q*b1  0    28*q*b1    0     -14*q*b1   0    4*b1*q  0     -b1/q];
    h1 = modulate2(g0,'b' );
    h0 = h;
    if lower(type(1)) == 'r'
        h1 = modulate2(h ,'b');
        h0 = g0;    
    end
         
    
    
    case {'cd', '7-9'}	% by Cohen and Daubechies
	% 1D prototype filters: the '7-9' pair
	h0 = [0.026748757411, -0.016864118443, -0.078223266529, ...
	      0.266864118443, 0.602949018236, 0.266864118443, ...
	      -0.078223266529, -0.016864118443, 0.026748757411];
	g0 = [-0.045635881557, -0.028771763114, 0.295635881557, ...
	      0.557543526229, 0.295635881557, -0.028771763114, ...
	      -0.045635881557];
	
	if lower(type(1)) == 'd'
	    h1 = modulate2(g0, 'c');
	else
	    h1 = modulate2(h0, 'c');
	    h0 = g0;
	end

	% Use McClellan to obtain 2D filters
	t = [0, 1, 0; 1, 0, 1; 0, 1, 0] / 4;	% diamond kernel
	h0 = sqrt(2) * ftrans2(h0, t);		
	h1 = sqrt(2) * ftrans2(h1, t);
	
    case {'pkva', 'ldtest'}	% Filters from the ladder structure	
	% Allpass filter for the ladder structure network
	beta = ldfilter(fname);
	
	% Analysis filters
	[h0, h1] = ld2quin(beta);
	
	% Normalize norm
	h0 = sqrt(2) * h0;
	h1 = sqrt(2) * h1;
	
	% Synthesis filters
	if lower(type(1)) == 'r'
	    f0 = modulate2(h1, 'b');
	    f1 = modulate2(h0, 'b');
	    
	    h0 = f0;
	    h1 = f1;
	end	
	
 case {'pkva-half4'}	% Filters from the ladder structure	
	% Allpass filter for the ladder structure network
	beta = ldfilterhalf( 4);
	
	% Analysis filters
	[h0, h1] = ld2quin(beta);
	
	% Normalize norm
	h0 = sqrt(2) * h0;
	h1 = sqrt(2) * h1;
	
	% Synthesis filters
	if lower(type(1)) == 'r'
	    f0 = modulate2(h1, 'b');
	    f1 = modulate2(h0, 'b');
	    
	    h0 = f0;
	    h1 = f1;
	end	
	
    case {'pkva-half6'}	% Filters from the ladder structure	
	% Allpass filter for the ladder structure network
	beta = ldfilterhalf( 6);
	
	% Analysis filters
	[h0, h1] = ld2quin(beta);
	
	% Normalize norm
	h0 = sqrt(2) * h0;
	h1 = sqrt(2) * h1;
	
	% Synthesis filters
	if lower(type(1)) == 'r'
	    f0 = modulate2(h1, 'b');
	    f1 = modulate2(h0, 'b');
	    
	    h0 = f0;
	    h1 = f1;
	end	
    
    case {'pkva-half8'}	% Filters from the ladder structure	
	% Allpass filter for the ladder structure network
	beta = ldfilterhalf( 8);
	
	% Analysis filters
	[h0, h1] = ld2quin(beta);
	
	% Normalize norm
	h0 = sqrt(2) * h0;
	h1 = sqrt(2) * h1;
	
	% Synthesis filters
	if lower(type(1)) == 'r'
	    f0 = modulate2(h1, 'b');
	    f1 = modulate2(h0, 'b');
	    
	    h0 = f0;
	    h1 = f1;
	end	
    
        
    case 'sinc'	% The "sinc" case, NO Perfect Reconstruction
	% Ideal low and high pass filters
	flength = 30;
	
	h0 = fir1(flength, 0.5);
	h1 = modulate2(h0, 'c');
	
	% Use McClellan to obtain 2D filters
	t = [0, 1, 0; 1, 0, 1; 0, 1, 0] / 4;	% diamond kernel
	h0 = sqrt(2) * ftrans2(h0, t);	
	h1 = sqrt(2) * ftrans2(h1, t);	

    case {'oqf_362'}	% Some "home-made" filters!
        h0 = sqrt(2) / 64 * ...
             [	sqrt(15),	-3,	0; ...
                 0,		5,	-sqrt(15); ...
                 -2*sqrt(15),	30,	0; ...
                 0,		30,	2*sqrt(15); ...
                 sqrt(15),	5, 	0; ...
                 0, 		-3,	-sqrt(15)]';

         h1 = -reverse2(modulate2(h0, 'b'));
	
	if type == 'r'
	    % Reverse filters for reconstruction
	    h0 = h0(end:-1:1, end:-1:1);
	    h1 = h1(end:-1:1, end:-1:1);
	end
    
    
    	
    case {'test'}	% Only for the shape, not for PR
	h0 = [0, 1, 0; 1, 4, 1; 0, 1, 0];
	h1 = [0, -1, 0; -1, 4, -1; 0, -1, 0];
	
    case {'testDVM'}	% Only for directional vanishing moment
	h0 = [1, 1; 1, 1] / sqrt(2);
	h1 = [-1, 1; 1, -1] / sqrt(2);	
	
 	
    case 'qmf'	% by Lu, Antoniou and Xu
	% ideal response
    % window
    m=2;
    n=2;
    w=[];
    w1d=kaiser(4*m+1,2.6);
    for n1=-m:m
        for n2 = -n:n
            w(n1+m+1,n2+n+1) = w1d(2*m+n1+n2+1)*w1d(2*m+n1-n2+1) ;
        end
    end
    h=[];
    for n1=-m:m
        for n2 = -n:n
            h(n1+m+1,n2+n+1) =  .5*sinc((n1+n2)/2)*.5*sinc((n1-n2)/2);
        end
    end
    c=sum(sum(h));
    h = sqrt(2)*h/c;
    h0=h.*w;
    h1 = modulate2(h0,'b');
    
    %h0 = modulate2(h,'r');
    %h1 = modulate2(h,'b');
    
    
     case 'qmf2'	% by Lu, Antoniou and Xu
	% ideal response
    % window
    
    h=[-.001104 .002494 -0.001744 0.004895 -0.000048 -.000311;
    0.008918 -0.002844 -0.025197 -0.017135 0.003905 -0.000081;
    -0.007587 -0.065904 00.100431 -0.055878 0.007023 0.001504;
    0.001725 0.184162 0.632115 0.099414 -0.027006 -0.001110;
    -0.017935 -0.000491 0.191397 -0.001787 -0.010587 0.002060;
    .001353 0.005635 -0.001231 -0.009052 -0.002668 0.000596];
    h0 = h./sum(sum(h));
    h1 = modulate2(h0,'b');
     
    %h0 = modulate2(h,'r');
    %h1 = modulate2(h,'b');   
    
   case 'dmaxflat4'
       M1 = 1/sqrt(2); M2=M1;
       k1=1-sqrt(2);k3=k1;
       k2 = M1;          
       h  = [.25*k2*k3 .5*k2 1+.5*k2*k3]*M1; h = [h fliplr(h(1:end-1))];
       g  = [-.125*k1*k2*k3 0.25*k1*k2 (-0.5*k1-0.5*k3-0.375*k1*k2*k3) 1+ .5*k1*k2]*M2;
       g = [g fliplr(g(1:end-1))];
       B  = dmaxflat(4,0);
       h0 = mctrans(h,B);
       g0 = mctrans(g,B);
       h0 = sqrt(2)*(h0./sum(h0(:)));
       g0 = sqrt(2)*(g0./sum(g0(:)));
       
       h1 = modulate2(g0,'b' );
    if lower(type(1)) == 'r'
       h1 = modulate2(h0 ,'b');
       h0 = g0;    
    end
   case 'dmaxflat5'
       M1 = 1/sqrt(2); M2=M1;
       k1=1-sqrt(2);k3=k1;
       k2 = M1;          
       h  = [.25*k2*k3 .5*k2 1+.5*k2*k3]*M1; h = [h fliplr(h(1:end-1))];
       g  = [-.125*k1*k2*k3 0.25*k1*k2 (-0.5*k1-0.5*k3-0.375*k1*k2*k3) 1+ .5*k1*k2]*M2;
       g = [g fliplr(g(1:end-1))];
       B  = dmaxflat(5,0);
       h0 = mctrans(h,B);
       g0 = mctrans(g,B);
       h0 = sqrt(2)*(h0./sum(h0(:)));
       g0 = sqrt(2)*(g0./sum(g0(:)));
       
       h1 = modulate2(g0,'b' );
    if lower(type(1)) == 'r'
       h1 = modulate2(h0 ,'b');
       h0 = g0;    
    end
    
   case 'dmaxflat6'
       M1 = 1/sqrt(2); M2=M1;
       k1=1-sqrt(2);k3=k1;
       k2 = M1;          
       h  = [.25*k2*k3 .5*k2 1+.5*k2*k3]*M1; h = [h fliplr(h(1:end-1))];
       g  = [-.125*k1*k2*k3 0.25*k1*k2 (-0.5*k1-0.5*k3-0.375*k1*k2*k3) 1+ .5*k1*k2]*M2;
       g = [g fliplr(g(1:end-1))];
       B  = dmaxflat(6,0);
       h0 = mctrans(h,B);
       g0 = mctrans(g,B);
       h0 = sqrt(2)*(h0./sum(h0(:)));
       g0 = sqrt(2)*(g0./sum(g0(:)));
       
       h1 = modulate2(g0,'b' );
    if lower(type(1)) == 'r'
       h1 = modulate2(h0 ,'b');
       h0 = g0;    
    end
    
   case 'dmaxflat7'
       M1 = 1/sqrt(2); M2=M1;
       k1=1-sqrt(2);k3=k1;
       k2 = M1;          
       h  = [.25*k2*k3 .5*k2 1+.5*k2*k3]*M1; h = [h fliplr(h(1:end-1))];
       g  = [-.125*k1*k2*k3 0.25*k1*k2 (-0.5*k1-0.5*k3-0.375*k1*k2*k3) 1+ .5*k1*k2]*M2;
       g = [g fliplr(g(1:end-1))];
       B  = dmaxflat(7,0);
       h0 = mctrans(h,B);
       g0 = mctrans(g,B);
       h0 = sqrt(2)*(h0./sum(h0(:)));
       g0 = sqrt(2)*(g0./sum(g0(:)));
       
       h1 = modulate2(g0,'b' );
    if lower(type(1)) == 'r'
       h1 = modulate2(h0 ,'b');
       h0 = g0;    
    end
      
    
    otherwise
	% Assume the "degenerated" case: 1D wavelet filters
	[h0, h1] = wfilters(fname, type); 
end