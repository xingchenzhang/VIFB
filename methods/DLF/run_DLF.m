% This is the code for DLF image fusion algorithm:
% H. Li, X.-J. Wu, and J. Kittler, ¡°Infrared and visible image
% fusion using a deep learning framework,¡± 24th International
% Conference on Pattern Recognition, 2018.
% More details about DLF can be found from the README.md file in DLF folder
%
% The interface is created by the authors of VIFB. Necessary changes are
% also made to be integrated in VIFB.

function img = run_DLF(imgVI, imgIR, visualization)

    % add the matconvnet path
    addpath(genpath('Your own path\matconvnet-1.0-beta25\matlab'));
    vl_setupnn ;

    %load vgg19
    net = load('imagenet-vgg-verydeep-19.mat');
    net = vl_simplenn_tidy(net);

    n = 21;
        
    % IR image
    path1 = imgIR.img;

    % VI image
    path2 = imgVI.img;  
        
    image1 = imread(path1);
    image2 = imread(path2);
                
    % IR
    image1 = im2double(image1);
    
    % VI
    image2 = im2double(image2);
   
    tic;    
    if size(image2,3)==1         
        img = DLF(image1, image2, net);  
    elseif size(image1,3) == 1    
        img = zeros(size(image2)); 
        for i=1:3
            img(:,:,i) = DLF(image1, image2(:,:,i), net);  
        end    
    else
        img = zeros(size(image2)); 
        for i=1:3
            img(:,:,i) = DLF(image1(:,:,i), image2(:,:,i), net);  
        end    
    end
    toc;  
    img=im2uint8(img);     
end

    
function res = DLF(image1, image2, net)   
        
    % Highpass filter test image
    npd = 16;
    fltlmbd = 5;    
        
    [I_lrr1, I_saliency1] = lowpass(image1, fltlmbd, npd);
    [I_lrr2, I_saliency2] = lowpass(image2, fltlmbd, npd);

    %% fuison lrr parts
    F_lrr = (I_lrr1+I_lrr2)/2;
    %figure;imshow(F_lrr);

    %% fuison saliency parts use VGG19
%    disp('VGG19-saliency');
    saliency_a = make_3c(I_saliency1);
    saliency_b = make_3c(I_saliency2);
%      saliency_a = I_saliency1;
%      saliency_b = I_saliency2;
    
    saliency_a = single(saliency_a) ; % note: 255 range
    saliency_b = single(saliency_b) ; % note: 255 range

    res_a = vl_simplenn(net, saliency_a);
    res_b = vl_simplenn(net, saliency_b);

    %% relu1_1
    out_relu1_1_a = res_a(2).x;
    out_relu1_1_b = res_b(2).x;
    unit_relu1_1 = 1;

    l1_featrues_relu1_a = extract_l1_feature(out_relu1_1_a);
    l1_featrues_relu1_b = extract_l1_feature(out_relu1_1_b);

    [F_saliency_relu1, l1_featrues_relu1_ave_a, l1_featrues_relu1_ave_b] = ...
                fusion_strategy(l1_featrues_relu1_a, l1_featrues_relu1_b, I_saliency1, I_saliency2, unit_relu1_1);

    %% relu2_1
    out_relu2_1_a = res_a(7).x;
    out_relu2_1_b = res_b(7).x;
    unit_relu2_1 = 2;

    l1_featrues_relu2_a = extract_l1_feature(out_relu2_1_a);
    l1_featrues_relu2_b = extract_l1_feature(out_relu2_1_b);

    [F_saliency_relu2, l1_featrues_relu2_ave_a, l1_featrues_relu2_ave_b] = ...
                fusion_strategy(l1_featrues_relu2_a, l1_featrues_relu2_b, I_saliency1, I_saliency2, unit_relu2_1);

    %% relu3_1
    out_relu3_1_a = res_a(12).x;
    out_relu3_1_b = res_b(12).x;
    unit_relu3_1 = 4;

    l1_featrues_relu3_a = extract_l1_feature(out_relu3_1_a);
    l1_featrues_relu3_b = extract_l1_feature(out_relu3_1_b);

    [F_saliency_relu3, l1_featrues_relu3_ave_a, l1_featrues_relu3_ave_b] = ...
                fusion_strategy(l1_featrues_relu3_a, l1_featrues_relu3_b, I_saliency1, I_saliency2, unit_relu3_1);

    %% relu4_1
    out_relu4_1_a = res_a(21).x;
    out_relu4_1_b = res_b(21).x;
    unit_relu4_1 = 8;

    l1_featrues_relu4_a = extract_l1_feature(out_relu4_1_a);
    l1_featrues_relu4_b = extract_l1_feature(out_relu4_1_b);

    [F_saliency_relu4, l1_featrues_relu4_ave_a, l1_featrues_relu4_ave_b] = ...
                fusion_strategy(l1_featrues_relu4_a, l1_featrues_relu4_b, I_saliency1, I_saliency2, unit_relu4_1);

    %% fusion strategy
    F_saliency = max(F_saliency_relu1, F_saliency_relu2);
    F_saliency = max(F_saliency, F_saliency_relu3);
    F_saliency = max(F_saliency, F_saliency_relu4);

    fusion_im = F_lrr + F_saliency;

    res = fusion_im; 

end
