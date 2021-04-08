function bgImg = QuadReconstructRefined(S, img, minDim)
% This function is used to reconstruct the background image using bezier
% interpolation interpolation on the quadtree structure S.
% Implemented by Zhang Yu (uzeful@163.com).

    % double image 
    img = double(img);

    % Matrix M
    M = [
        -1 3 -3 1; 
       	3 -6 3 0;
      	-3 3 0 0;
    	1 0 0 0
        ];

    % preprocess the obtained quadtree structure S
    newS = S + (S > 0);

    % max level of the quadtree structure
    newS = full(newS);
    maxDim = max(newS(:));
    dim = maxDim;

    % pad image for the following computation
    newS = padarray(newS, [1 1], 'replicate', 'post');
    newImg = padarray(img, [1 1], 'replicate', 'post');

    % temporal reconstructed image
    tempReconstImg = zeros(size(newImg));

    % reconstruct background image from each scale blocks of the quadtree structure
    while (dim >= minDim + 1)
        % block number of the current scale blocks
        len = length(find(newS == dim));
        if len ~= 0
            % Extrat the corresponding blocks at the current level
            [blks, Sind] = qtgetblk(newImg, newS, dim);

            % pixel locations
            subDim = (dim - 1) / 4;
            row = [1, subDim * 2, subDim * 3, subDim * 4 + 1];
            
            xx = [row; row; row; row];
            yy = xx';
            inds = sub2ind([dim, dim], yy, xx);

            % generate u and v
            u = ([1 : dim]' - 1) / (dim - 1);
            U = [u.^3, u.^2, u, u.^0];  % V = U';

            % precompute the left matrix(LM) U*M and right matrix(RM) M'V'
            LM = U * M;
            RM = M * U'; % M'=M and V = U
            
            % reconstructed blocks
            reBlkSeq = zeros(dim, dim, len);

            for ii = 1 : len
                blockVal = blks( : , : , ii);
                blockVal = blockVal(inds);
                reblkVal = LM * blockVal * RM;
                reBlkSeq(:,:,ii) = reblkVal;
            end

            % Set the reconstructed blocks into tempReconstImg
            tempReconstImg = qtsetblk(tempReconstImg,newS,dim,reBlkSeq);
        end
        dim = (dim - 1) / 2 + 1;
    end
    % Final reconstrcucted background image
    bgImg = tempReconstImg(1 : end - 1, 1 : end - 1);
end
