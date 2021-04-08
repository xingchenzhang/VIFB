% This is the code for CNN image fusion algorithm:
% Y. Liu, X. Chen, J. Cheng, H. Peng, and Z. Wang, ¡°Infrared and visible image fusion with convolutional 
% neural networks,¡± International Journal of Wavelets, Multiresolution and Information Processing, vol. 16, 
% no. 03, p. 1850018, 2018.
%
% The code of CNN is provided by authors of CNN.
% The interface is created by authors of VIFB. Besides, necessary changes are made to adapt to VIFB

function img = run_CNN(imgVI, imgIR, visualization)
    
    %% load model
    addpath(genpath('pyramid'));
    model =  'model/cnnmodel.mat';
    load(model);
    
    [conv1_patchsize2,conv1_filters] = size(weights_b1_1);  %9*64
    conv1_patchsize = sqrt(conv1_patchsize2);
    [conv2_channels,conv2_patchsize2,conv2_filters] = size(weights_b1_2);  %64*9*128
    conv2_patchsize = sqrt(conv2_patchsize2);
    [conv3_channels,conv3_patchsize2,conv3_filters] = size(weights_b1_3);  %128*9*256
    conv3_patchsize = sqrt(conv3_patchsize2);
    [conv4_channels,conv4_patchsize2,conv4_filters] = size(weights_output);  %512*64*2
    conv4_patchsize = sqrt(conv4_patchsize2);
        
    % IR image
    img1 = imread(imgIR.img);

    % VI image
    img2 = imread(imgVI.img);

    if visualization == 1
        figure;imshow(img1);
        figure;imshow(img2);
    end
    
    % three channels
    if size(img1,3)>1
        imgIR_1 = img1(:,:,1);
        imgIR_2 = img1(:,:,2);
        imgIR_3 = img1(:,:,3);
    else
        imgIR_1 = img1;
        imgIR_2 = img1;
        imgIR_3 = img1;
    end

    % three channels
    if size(img2,3)>1
        imgVI_1 = img2(:,:,1);
        imgVI_2 = img2(:,:,2);
        imgVI_3 = img2(:,:,3);
    else
        imgVI_1 = img2;
        imgVI_2 = img2;
        imgVI_3 = img2;
    end

    tic;

    img1 = double(img1)/255;
    img2 = double(img2)/255;

    imgIR_1 = double(imgIR_1)/255;
    imgIR_2 = double(imgIR_2)/255;
    imgIR_3 = double(imgIR_3)/255;

    imgVI_1 = double(imgVI_1)/255;
    imgVI_2 = double(imgVI_2)/255;
    imgVI_3 = double(imgVI_3)/255;

    if size(img1,3)>1
        img1_gray=rgb2gray(img1);
    else
        img1_gray=img1;
    end

    if size(img2,3)>1
        img2_gray=rgb2gray(img2);
    else
        img2_gray=img2;
    end

    [hei, wid] = size(img1_gray);
    
    %% conv1
    weights_conv1 = reshape(weights_b1_1, conv1_patchsize, conv1_patchsize, conv1_filters);
    conv1_data1 = zeros(hei, wid, conv1_filters,'single');
    conv1_data2 = zeros(hei, wid, conv1_filters,'single');
    for i = 1 : conv1_filters
        conv1_data1(:,:,i) = conv2(img1_gray, rot90(weights_conv1(:,:,i),2), 'same');
        conv1_data1(:,:,i) = max(conv1_data1(:,:,i) + biases_b1_1(i), 0);
        conv1_data2(:,:,i) = conv2(img2_gray, rot90(weights_conv1(:,:,i),2), 'same');
        conv1_data2(:,:,i) = max(conv1_data2(:,:,i) + biases_b1_1(i), 0);
    end
    %% conv2
    conv2_data1 = zeros(hei, wid, conv2_filters,'single');
    conv2_data2 = zeros(hei, wid, conv2_filters,'single');
    for i = 1 : conv2_filters
        for j = 1 : conv2_channels
            conv2_subfilter = rot90(reshape(weights_b1_2(j,:,i), conv2_patchsize, conv2_patchsize),2);
            conv2_data1(:,:,i) = conv2_data1(:,:,i) + conv2(conv1_data1(:,:,j), conv2_subfilter, 'same');
            conv2_data2(:,:,i) = conv2_data2(:,:,i) + conv2(conv1_data2(:,:,j), conv2_subfilter, 'same');
        end
        conv2_data1(:,:,i) = max(conv2_data1(:,:,i) + biases_b1_2(i), 0);
        conv2_data2(:,:,i) = max(conv2_data2(:,:,i) + biases_b1_2(i), 0);
    end

    %% max-pooling2
    conv2_data1_pooling=zeros(ceil(hei/2), ceil(wid/2), conv2_filters,'single');
    conv2_data2_pooling=zeros(ceil(hei/2), ceil(wid/2), conv2_filters,'single');
    for i = 1 : conv2_filters    
        conv2_data1_pooling(:,:,i) = maxpooling_s2(conv2_data1(:,:,i));
        conv2_data2_pooling(:,:,i) = maxpooling_s2(conv2_data2(:,:,i));
    end
    %% conv3
    conv3_data1 = zeros(ceil(hei/2), ceil(wid/2), conv3_filters,'single');
    conv3_data2 = zeros(ceil(hei/2), ceil(wid/2), conv3_filters,'single');
    for i = 1 : conv3_filters
        for j = 1 : conv3_channels
            conv3_subfilter = rot90(reshape(weights_b1_3(j,:,i), conv3_patchsize, conv3_patchsize),2);
            conv3_data1(:,:,i) = conv3_data1(:,:,i) + conv2(conv2_data1_pooling(:,:,j), conv3_subfilter, 'same');
            conv3_data2(:,:,i) = conv3_data2(:,:,i) + conv2(conv2_data2_pooling(:,:,j), conv3_subfilter, 'same');
        end
        conv3_data1(:,:,i) = max(conv3_data1(:,:,i) + biases_b1_3(i), 0);
        conv3_data2(:,:,i) = max(conv3_data2(:,:,i) + biases_b1_3(i), 0);
    end

    %% feature layer
    conv3_data=cat(3,conv3_data1,conv3_data2);
    conv4_data=zeros(ceil(hei/2)-conv4_patchsize+1,ceil(wid/2)-conv4_patchsize+1,conv4_filters,'single');
    for i = 1 : conv4_filters
        for j = 1 : conv4_channels
            conv4_subfilter = rot90((reshape(weights_output(j,:,i), conv4_patchsize, conv4_patchsize)),2);
            conv4_data(:,:,i) = conv4_data(:,:,i) + conv2(conv3_data(:,:,j), conv4_subfilter, 'valid');
        end
    end 
    %% softmax ouput layer
    conv4_data=double(conv4_data);
    output_data=zeros(ceil(hei/2)-conv4_patchsize+1,ceil(wid/2)-conv4_patchsize+1,conv4_filters);
    output_data(:,:,1)=exp(conv4_data(:,:,1))./(exp(conv4_data(:,:,1))+exp(conv4_data(:,:,2)));
    output_data(:,:,2)=1-output_data(:,:,1);
    outMap=output_data(:,:,2);

    %% focus map generation
    sumMap=zeros(hei,wid);
    cntMap=zeros(hei,wid);
    patch_size=16;
    temp_size_y=patch_size;
    temp_size_x=patch_size;
    stride=2;
    y_bound=hei-patch_size+1;
    x_bound=wid-patch_size+1;

    [h,w]=size(outMap);
    for j=1:h
        jj=(j-1)*stride+1;
        if jj<=y_bound
            temp_size_y=patch_size;
        else
            temp_size_y=hei-jj+1;
        end
        for i=1:w
            ii=(i-1)*stride+1;
            if ii<=x_bound
                temp_size_x=patch_size;
            else
                temp_size_x=wid-ii+1;
            end
            sumMap(jj:jj+temp_size_y-1,ii:ii+temp_size_x-1)=sumMap(jj:jj+temp_size_y-1,ii:ii+temp_size_x-1)+outMap(j,i);
            cntMap(jj:jj+temp_size_y-1,ii:ii+temp_size_x-1)=cntMap(jj:jj+temp_size_y-1,ii:ii+temp_size_x-1)+1;
        end
    end

    focusMap=sumMap./cntMap;

    %% LP
    %if size(img1,3)>1
    %     weightMap=repmat(focusMap,[1 1 3]);
    %else
        weightMap=focusMap;
    %end

    pyr = gaussian_pyramid(zeros(hei,wid));
    nlev = length(pyr);

    pyrW=gaussian_pyramid(weightMap,nlev);
    % pyrI1=laplacian_pyramid(img1,nlev);
    % pyrI2=laplacian_pyramid(img2,nlev);

    if (size(img1,3)>1) || (size(img2,3)>1)
        disp('at least one figure has 3 channels') ;
        % the first channel
        pyrI1=laplacian_pyramid(imgIR_1,nlev);
        pyrI2=laplacian_pyramid(imgVI_1,nlev);

        for l = 1:nlev
           pyr{l}=band_fuse(pyrI1{l},pyrI2{l},pyrW{l},0.6);
        end

        % reconstruct
        imgf_1 = reconstruct_laplacian_pyramid(pyr);

        % the second channel
        pyrI3=laplacian_pyramid(imgIR_2,nlev);
        pyrI4=laplacian_pyramid(imgVI_2,nlev);

        for l = 1:nlev
           pyr{l}=band_fuse(pyrI3{l},pyrI4{l},pyrW{l},0.6);
        end

        % reconstruct
        imgf_2 = reconstruct_laplacian_pyramid(pyr);

        % the third channel
        pyrI5=laplacian_pyramid(imgIR_3,nlev);
        pyrI6=laplacian_pyramid(imgVI_3,nlev);

        for l = 1:nlev
           pyr{l}=band_fuse(pyrI5{l},pyrI6{l},pyrW{l},0.6);
        end

        % reconstruct
        imgf_3 = reconstruct_laplacian_pyramid(pyr);

        imgf(:,:,1) = imgf_1;
        imgf(:,:,2) = imgf_2;
        imgf(:,:,3) = imgf_3;

    else
        disp('both figures have only 1 channel') ;
        pyrI1=laplacian_pyramid(img1,nlev);
        pyrI2=laplacian_pyramid(img2,nlev);

        for l = 1:nlev
           pyr{l}=band_fuse(pyrI1{l},pyrI2{l},pyrW{l},0.6);
        end

        % reconstruct
        imgf = reconstruct_laplacian_pyramid(pyr);
    end

    toc;
    if visualization == 1
        figure;imshow(uint8(imgf*255));
    end
    img = uint8(imgf*255);
end
