% This code can be used to compute values of selected metrics for selected
% algorithms on selected image pair
%
% VIFB is the first benchmark in the field of visible-infrared image fsuion, and also the first benchmark in 
% the field of image fusion.
%
% Note: Please change the path in line 30 to your own path before running
%
% If you use this code, please site the following paper:
%
% X. Zhang, P. Ye, G. Xiao. VIFB:A Visible and Infrared Benchmark. In
% Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern
% Recognition Workshops, 2020.
%
% Thanks a lot!
%
% For more information, please see https://github.com/xingchenzhang/VIFB
%
% Contact: xingchen.zhang@imperial.ac.uk

close all
clear
clc
warning off all;

addpath('.\metrics');
addpath('.\util');
addpath('.\methods');

path = 'Your own path\VIFB\output';
fusedPath = [path '\fused_images\'];
outputPath = [path '\evaluation_metrics\'];
outputPathSingle = [path '\evaluation_metrics_single\'];

if ~exist(outputPath,'dir')
    mkdir(outputPath);
end

if ~exist(outputPathSingle,'dir')
    mkdir(outputPathSingle);
end

imgsVI = configImgsVI;
imgsIR = configImgsIR;
methods = configMethods;
metrics = configMetrics;

numImgsVI = length(imgsVI);
numImgsIR = length(imgsIR);
numMethods = length(methods);
numMetrics = length(metrics);

% output information
fid = fopen(strcat(path, '\information.txt'),'w');
fprintf(fid,'%15s\r\n','The image paris are:');

for i=1:numImgsVI
    fprintf(fid ,'%15s\r\n',imgsVI{i}.name);
end

fprintf(fid,'%15s\r\n','');
fprintf(fid,'%15s\r\n','The methods are:');
for i=1:numMethods
    fprintf(fid,'%15s\r\n', methods{i}.name);
end

fprintf(fid,'%15s\r\n','');
fprintf(fid,'%15s\r\n','The metrics are:');
for i=1:numMetrics
    fprintf(fid,'%15s\r\n', metrics{i}.name);
end
fclose(fid);

visualization = 0;
resultsMetrics = zeros(numImgsVI, numMethods, numMetrics);
            
for idxMethod = 1:numMethods
    m = methods{idxMethod};

    for idxImgs = 1:length(imgsVI)
        sVI = imgsVI{idxImgs};
        sIR = imgsIR{idxImgs};

        sVI.img = strcat(sVI.path,sVI.name, '.',sVI.ext);
        sIR.img = strcat(sIR.path,sIR.name, '.',sIR.ext);

        imgVI = imread(sVI.img);
        imgIR = imread(sIR.img);

        [imgH_VI,imgW_VI,chVI] = size(imgVI);
        [imgH_IR,imgW_IR,chIR] = size(imgIR);
        
        for idxMetrics = 1:numMetrics
            
            sMetrics = metrics{idxMetrics};
        
            fusedName = [fusedPath sVI.name '_' m.name '.jpg'];
            if exist([fusedPath sVI.name '_' m.name '.jpg'])
                sFused = imread(fusedName);              
                % check whether the result exists
                if exist(strcat(outputPathSingle,sVI.name, '_', m.name,'_',sMetrics.name ,'.txt'))                    
                    A = importdata(strcat(outputPathSingle,sVI.name, '_', m.name,'_',sMetrics.name ,'.txt'));              
                    resultsMetrics(idxImgs, idxMethod, idxMetrics) = A;
                    continue;
                end
                
                disp([num2str(idxMethod) '_' m.name ', ' num2str(idxImgs) '_' sVI.name ', ' num2str(idxMetrics) '_' sMetrics.name])       

                funcName = ['res = metrics' sMetrics.name '(imgVI, imgIR, sFused);'];
                disp(funcName)

                try
                    cd(['./metrics/']);
                    addpath(genpath('./'))

                    eval(funcName);

                catch err
                    disp('error');
                    rmpath(genpath('./'))
                    cd('../../')
                    continue;
                end
                
                resultsMetrics(idxImgs, idxMethod, idxMetrics) = res;
                
                outputFileSingle = strcat(outputPathSingle,sVI.name, '_', m.name,'_',sMetrics.name ,'.txt');

                dlmwrite(outputFileSingle,res)
                cd('../');
                
            else                
               str=['The fused image ' fusedName ' does not exists, please check'];
               disp(str)           
            end            
        end               
    end  
end

outputFile = strcat(outputPath, 'evaluationMetrics.mat');
save(outputFile,'resultsMetrics');

% compute the average value of each metric on all image pairs
resultsMetricsAverageImg = nanmean(resultsMetrics,1); 
outputFileAverage = strcat(outputPath, 'evaluationMetricsAverageImg.mat');
save(outputFileAverage,'resultsMetricsAverageImg');
