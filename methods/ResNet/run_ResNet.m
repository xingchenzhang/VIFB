% This is the code of ResNet image fusion algorithm:
% H. Li, X.-J. Wu, and T. S. Durrani, ¡°Infrared and visible image fusion with resnet and zero-phase component analysis,¡±
% Infrared Physics & Technology, vol. 102, p. 103039, 2019.
%
% The interface is created by the authors of VIFB.

function img = run_ResNet(imgVI, imgIR, visualization)

    % add the matconvnet path
    addpath(genpath('your own path\matconvnet-1.0-beta25\matlab'));
    vl_setupnn ;

    %% ResNet50 + ZCA & norm(l1, l2, nuclear)
    % load the pre-trained model - ResNet-50
    model_path = '.\models\';
    % http://www.vlfeat.org/matconvnet/pretrained/
    net_ = load([model_path, 'imagenet-resnet-50-dag.mat']);
    net = dagnn.DagNN.loadobj(net_);
    net.mode = 'test';
    %% remove layers - ResNet50
    % Conv5 - res5cx
    for i = 173:175
        net.removeLayer(net.layers(173).name);
    end
    net_res5cx = net;
    % Conv4 - res4fx
    net = dagnn.DagNN.loadobj(net_);
    net.mode = 'test';
    for i = 141:175
        net.removeLayer(net.layers(141).name);
    end
    net_res4fx = net;

    %% Start
    tic;
    index = i;
    disp(num2str(index));
        
    % IR image
    path1 = imgIR.img;

    % VI image
    path2 = imgVI.img;  

    image1 = imread(path1);
    image2 = imread(path2);
        
    I1 = double(image1);
    I2 = double(image2);

    %% Extract features, run the net - ResNet50
    disp('ResNet');
    if size(image1, 3)<3
        I1 = make_3c(image1);
    end
    if size(image2, 3)<3
        I2 = make_3c(image2);
    end
    I1 = single(I1) ; % note: 255 range
    I2 = single(I2) ; % note: 255 range

    tic;
    % I1
    disp('run the ResNet - I1');
    net_res4fx.eval({'data', I1}) ;
    output4_1 = net_res4fx.vars(net_res4fx.getVarIndex('res4fx')).value ;
    net_res5cx.eval({'data', I1}) ;
    output5_1 = net_res5cx.vars(net_res5cx.getVarIndex('res5cx')).value ;
    % I2
    disp('run the ResNet - I2');
    net_res4fx.eval({'data', I2}) ;
    output4_2 = net_res4fx.vars(net_res4fx.getVarIndex('res4fx')).value ;
    net_res5cx.eval({'data', I2}) ;
    output5_2 = net_res5cx.vars(net_res5cx.getVarIndex('res5cx')).value ;

    %% extract features - ZCA & l1-norm operation
    disp('extract features(whitening operation) - I1');
    feature4_1 = whitening_norm(output4_1);
    feature5_1 = whitening_norm(output5_1);
    disp('extract features(whitening operation) - I2');
    feature4_2 = whitening_norm(output4_2);
    feature5_2 = whitening_norm(output5_2);

    %% fusion strategy - resize to original size and soft-max
    disp('fusion strategy(weighting)');
    % output4 - 1024
    %[F_relu4, weight4_a, weight4_b] = fusion_strategy(feature4_1, feature4_2, image1, image2);
    % output5 - 2048

    if size(I2, 3) == 1
        [F_relu5, weight5_a, weight5_b] = fusion_strategy(feature5_1, feature5_2, I1, I2);
    elseif size(I1,3) == 1
        F_relu5 = zeros(size(I2));
        for i=1:3
            [F_relu5(:,:,i), weight5_a, weight5_b] = fusion_strategy(feature5_1, feature5_2, I1, I2(:,:,i));
        end       
    else
        F_relu5 = zeros(size(I2));
        for i=1:3
            [F_relu5(:,:,i), weight5_a, weight5_b] = fusion_strategy(feature5_1, feature5_2, I1(:,:,i), I2(:,:,i));            
        end     
    end
    toc;
    img = uint8(F_relu5);
end


