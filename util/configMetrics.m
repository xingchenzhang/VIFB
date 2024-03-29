% This is used to set the evaluation metrics that to be computed
%
% Author:Xingchen Zhang, Ping Ye, Gang Xiao
% Contact: xingchen.zhang@imperial.ac.uk

function metrics=configMetrics

 % the same type of metrics are put together
  metricVIFB={struct('name','Cross_entropy'),...
     struct('name','Entropy'),...
     struct('name','Mutinf'),...
     struct('name','Psnr'),... % Information theory-based
     struct('name','Avg_gradient'),...
     struct('name','Edge_intensity'),...
     struct('name','Qabf'),...
     struct('name','Variance'),...
     struct('name','Spatial_frequency'),...% Image feature-based
     struct('name','Rmse'),...
     struct('name','Ssim'),... & Structural similarity-based
     struct('name','Qcb'),...
     struct('name','Qcv'),... % human perception-inspired
     };

 metrics = [metricVIFB];
