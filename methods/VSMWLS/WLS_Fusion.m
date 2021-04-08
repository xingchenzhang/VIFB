function fused = WLS_Fusion(I1,I2)

%% multi-scale decomposition using RGF and Gaussian filter
nLevel = 4; % levels
G1 = cell(1, nLevel + 1);
L1 = cell(1, nLevel + 1);
G1{1} = I1;
sigma_s = 2;
sigma_r = [0.05 0.05 0.05 0.05];
iteration = [4 4 4 4];
for i = 2 : nLevel
    G1{i} = RollingGuidanceFilter_Guided(G1{i-1},sigma_s,sigma_r(i-1),iteration(i-1));
    L1{i-1} = G1{i-1} - G1{i};
    sigma_s = 2 * sigma_s;
end
sigma_s = 2;
G1{nLevel + 1} = gaussFilter(G1{nLevel},sigma_s);
L1{nLevel} = G1{nLevel} - G1{nLevel + 1};
L1{nLevel+1} = G1{nLevel+1};


G2 = cell(1, nLevel + 1);
L2 = cell(1, nLevel + 1);
G2{1} = I2;
sigma_s = 2;
sigma_r = [0.05 0.05 0.05 0.05];
for i = 2 : nLevel
    G2{i} = RollingGuidanceFilter_Guided(G2{i-1},sigma_s,sigma_r(i-1),iteration(i-1));
    L2{i-1} = G2{i-1} - G2{i};
    sigma_s = 2 * sigma_s;
end
sigma_s = 2;
G2{nLevel + 1} = gaussFilter(G2{nLevel},sigma_s);
L2{nLevel} = G2{nLevel} - G2{nLevel + 1};
L2{nLevel+1} = G2{nLevel+1};

%% base layer fusion
weight1 = Visual_Weight_Map(I1);
weight2 = Visual_Weight_Map(I2);
BF = (0.5+0.5*(weight1-weight2)).*L1{nLevel+1} + (0.5+0.5*(weight2-weight1)).*L2{nLevel+1};

%% detail layer fusion
sigma0 = 2;
w = floor(3*sigma0);
C_0 = double(abs(L1{1}) < abs(L2{1}));
DF = C_0.*L2{1} + (1-C_0).*L1{1};

for i =  nLevel : -1 : 2
    
    w = floor(3*sigma0);
    h = fspecial('gaussian', [2*w+1, 2*w+1], sigma0);   
    C_0 = double(abs(L1{i}) < abs(L2{i}));
    C_0 = imfilter(C_0, h, 'symmetric');
    M = C_0.*L2{i} + (1-C_0).*L1{i};
    
    lambda = 0.01;  
    dd = Solve_Optimal(M,L1{i},L2{i},lambda);

    DF = DF + dd;
end


%% fused image
fused = BF + DF;


end


