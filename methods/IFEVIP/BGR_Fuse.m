function result = BGR_Fuse(imgVis, imgIR, QuadNormDim, QuadMinDim, GaussScale, MaxRatio, StdRatio)
% This is the core function of the paper "Infrared and Visual Image Fusion 
% through Infrared Feature Extraction and Visual Information.
% Implemented by Zhang Yu (uzeful@163.com).

    % Reconstruct background image using quadtree decomposition
    I = imresize(imgIR, [QuadNormDim QuadNormDim]);
    S = qtdecomp(I,.2,QuadMinDim);

    % reconstruct image pixels
    minImg = imerode(I, strel('square',QuadMinDim));
    bgImg = QuadReconstructRefined(S, minImg, QuadMinDim);
    bgImg = imresize(bgImg, size(imgIR));

    % Gaussian filtering
    h = fspecial('gaussian', GaussScale, GaussScale / 2);
    bgImg = imfilter(bgImg, h);

    % Extract infrared bright feature
    img_Inf = imgIR - uint8(bgImg);

    % entropy based feature fusion 
    addFeature = double(img_Inf) + (double(imgIR) - double(imgVis)) * (entropy(imgIR) / entropy(imgVis));

    % graylevel driven feature fusion
    addedVals = double(addFeature) .* (addFeature > 0) + double(imgVis);
    maxVals = sort(addedVals(:), 'descend');
    maxMean = mean(maxVals(1 : round(MaxRatio * length(maxVals))));
    ratio = min(255 / maxMean, StdRatio);

    result = imgVis + uint8(addFeature * ratio); 
end