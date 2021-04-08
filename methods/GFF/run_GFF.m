% This is the code for GFF image fusion algorithm:
% S. Li, X. Kang, and J. Hu, ¡°Image fusion with guided filtering,¡± IEEE Transactions on Image processing, vol. 22, no. 7,
% pp. 2864¨C2875, 2013
%
% The code of GFF is provided by the authors of GFF.
% The interface is created by the authors of VIFB.

function img=run_GFF(imgVI, imgIR,  visualization)
   
    % IR image
    I1 = imread(imgIR.img);

    % VI image
    I2 = imread(imgVI.img);
    
    tic;    
    path = 'Your own path\methods\GFF\sourceimages\colourset';
    
    if ~exist(path,'dir')
        mkdir(path);
    end
    
    imwrite(I2, [path '/' imgVI.name '1.jpg']); 
    imwrite(I1, [path '/' imgVI.name '2.jpg']); 
       
    I = load_images( './sourceimages/colourset',1); 
    
    rmdir(path, 's') %rmdir
       
    F = GFF(I);
    toc;
    if visualization == 1
        figure,imshow(F);
    end
    img = F;
end