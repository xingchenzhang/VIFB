function fuse_im=cbf_ieeeconf2011f(inp,detail,cov_wsize)

%%% method3_ieeeconf2011: Fuses the 2 images using the method discussed in "An Efficient Adaptive Fusion Scheme
%%% for Multifocus Images in Wavelet Domain using Statistical Properties of Neighborhood", 
%%% in Proc. of 14th Int. Conf. on Information Fusion, Jul. 2011, pp. 1-7.
%%%
%%% Details are obtained by subtracting original image by cross bilateral filter output. (My Method)
%%% These details are used to find weights (Edge Strength) for fusing the images contained in "inp".

%%% Author : B. K. SHREYAMSHA KUMAR
%%% Created on 15-02-2012.
%%% Updated on 16-02-2012.

half_wsize=(cov_wsize-1)/2; %%% Half window size excluding centre pixel.

temp1=per_extn_im_fn(detail{1},cov_wsize); 
temp2=per_extn_im_fn(detail{2},cov_wsize); 
[MM,NN]=size(temp1);
for ii=half_wsize+1:MM-half_wsize
   for jj=half_wsize+1:NN-half_wsize
      %%% 1st Detail Image.
      tt1=temp1(ii-half_wsize:ii+half_wsize,jj-half_wsize:jj+half_wsize);
      hr_cov_mat1=covarf(tt1,cov_wsize);
      ver_cov_mat1=covarf(tt1',cov_wsize);
      hor_es1=sum(eig(hr_cov_mat1));
      ver_es1=sum(eig(ver_cov_mat1));
      wt1(ii-half_wsize,jj-half_wsize)=hor_es1+ver_es1;
      %%% 2nd Detail Image.
      tt2=temp2(ii-half_wsize:ii+half_wsize,jj-half_wsize:jj+half_wsize);         
      hr_cov_mat2=covarf(tt2,cov_wsize);
      ver_cov_mat2=covarf(tt2',cov_wsize);         
      hor_es2=sum(eig(hr_cov_mat2));
      ver_es2=sum(eig(ver_cov_mat2));
      wt2(ii-half_wsize,jj-half_wsize)=hor_es2+ver_es2;
   end
end
[a b]=find(wt1==0);
for kk=1:length(a);
   wt1(a(kk),b(kk))=eps;
end
[a b]=find(wt2==0);
for kk=1:length(a);
   wt2(a(kk),b(kk))=eps;
end
fuse_im=(double(inp{1}).*wt1+double(inp{2}).*wt2)./(wt1+wt2);