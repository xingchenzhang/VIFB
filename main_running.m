% This is the main program of VIFB and can be used to produce fused images. 
%
% VIFB is the first benchmark in the field of visible-infrared image fsuion and also the first benchmark in 
% the field of image fusion.
%
% Note: Please change the path in line 28 to your own path before running
%
% If you use this code, please cite the following paper:
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

addpath('.\util');
addpath(['.\methods'])

path = 'Your own path\VIFB\output';
outputPath = [path '\fused_images\'];

imgsVI = configImgsVI;
imgsIR = configImgsIR;

methods = configMethods;

numImgsVI = length(imgsVI);
numImgsIR = length(imgsIR);
numMethods = length(methods);

if ~exist(outputPath,'dir')
    mkdir(outputPath);
end

visualization = 0;
            
for idxMethod = 1:numMethods
    m = methods{idxMethod};
    t1 = clock;

    j = 0;
    for idxImgs = 1:length(imgsVI)
        sVI = imgsVI{idxImgs};
        sIR = imgsIR{idxImgs};

        sVI.img = strcat(sVI.path,sVI.name, '.',sVI.ext);
        sIR.img = strcat(sIR.path,sIR.name, '.',sIR.ext);

        imgVI = imread(sVI.img);
        imgIR = imread(sIR.img);

        [imgH_VI,imgW_VI,chVI] = size(imgVI);
        [imgH_IR,imgW_IR,chIR] = size(imgIR);
        
        % check whether the result exists
        if exist([outputPath sVI.name '_' m.name '.jpg'])
            continue;
        end
        
        disp([num2str(idxMethod) '_' m.name ', ' num2str(idxImgs) '_' sVI.name])       

        funcName = ['img = run_' m.name '(sVI, sIR, visualization);'];

        try

            cd(['./methods/' m.name]);
            addpath(genpath('./'))
            
            eval(funcName);
            j=j+1;
            
        catch err
            disp('error');
            rmpath(genpath('./'))
            cd('../../')
            continue;
        end
        
        imwrite(img, [outputPath '/' sVI.name '_' m.name '.jpg']); 
        cd('../../');
    end
    
    t2 = clock;
    runtimeAverage = etime(t2,t1)./j;
        
    str=['The total runtime of ' m.name ' is: ' num2str(etime(t2,t1)) 's'];
    disp(str)
    
    str=['The average runtime of ' m.name ' per image is: ' num2str(runtimeAverage) 's'];
    disp(str)
end
