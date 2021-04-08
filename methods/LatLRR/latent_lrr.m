%% Latent Low-Rank Representation
function [Z,L,E] = latent_lrr(X,lambda)
% Latent Low-Rank Representation for Subspace Segmentation and Feature Extraction
% Guangcan Liu, Shuicheng Yan. ICCV 2011.
% Problem:
%   min_Z,L,E ||Z||_* + ||L||_* +¡¡lambda||E||_1,
%           s.t. X = XZ + LX +E.
% Solning problem by Inexact ALM

A = X;
tol = 1e-6;
rho = 1.1;
max_mu = 1e6;
mu = 1e-6;
maxIter = 1e6;
[d n] = size(X);
m = size(A,2);
atx = X'*X;
inv_a = inv(A'*A+eye(m));
inv_b = inv(A*A'+eye(d));
%% Initializing optimization variables
J = zeros(m,n);
Z = zeros(m,n);
L = zeros(d,d);
S = zeros(d,d);

% E = sparse(d,n);
E = zeros(d,n);

Y1 = zeros(d,n);
Y2 = zeros(m,n);
Y3 = zeros(d,d);

%% Start main loop
iter = 0;
disp('initial');

while iter<maxIter
    iter = iter + 1;
%     disp(['iter====>>>>' num2str(iter)]);
    %updating J by the Singular Value Thresholding(SVT) operator
    temp_J = Z + Y2/mu;
    [U_J,sigma_J,V_J] = svd(temp_J,'econ');
    sigma_J = diag(sigma_J);
    svp_J = length(find(sigma_J>1/mu));
    if svp_J>=1
        sigma_J = sigma_J(1:svp_J)-1/mu;
    else
        svp_J = 1;
        sigma_J = 0;
    end
    J = U_J(:,1:svp_J)*diag(sigma_J)*V_J(:,1:svp_J)';
    
    %updating S by the Singular Value Thresholding(SVT) operator
    temp_S = L + Y3/mu;
    [U_S,sigma_S,V_S] = svd(temp_S,'econ');
    sigma_S = diag(sigma_S);
    svp_S = length(find(sigma_S>1/mu));
    if svp_S>=1
        sigma_S = sigma_S(1:svp_S)-1/mu;
    else
        svp_S = 1;
        sigma_S = 0;
    end
    S = U_S(:,1:svp_S)*diag(sigma_S)*V_S(:,1:svp_S)';
    
    %udpate Z
    Z = inv_a*(atx-X'*L*X-X'*E+J+(X'*Y1-Y2)/mu);
    
    %udpate L
    L = ((X-X*Z-E)*X'+S+(Y1*X'-Y3)/mu)*inv_b;
    
    %update E
    xmaz = X-X*Z-L*X;
    temp = xmaz+Y1/mu;
    E = max(0,temp - lambda/mu)+min(0,temp + lambda/mu);
    
    leq1 = xmaz-E;
    leq2 = Z-J;
    leq3 = L-S;
    max_l1 = max(max(abs(leq1)));
    max_l2 = max(max(abs(leq2)));
    max_l3 = max(max(abs(leq3)));
    
    stopC1 = max(max_l1, max_l2);
    stopC = max(stopC1, max_l3);
    if stopC<tol
        disp('LRR done.');
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        Y3 = Y3 + mu*leq3;
        mu = min(max_mu,mu*rho);
    end
    
end
end