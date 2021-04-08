function F_BHE_g = BHEg_gray( F )

% clear all;
% close all;
% clc;
x0=F;
% figure(1);
% imshow(x0);
%x0=imread('om.bmp');
%x0=imread('11.bmp');
% x0=imread('3.bmp');
% y0=rgb2ycbcr(x0);
 luma=x0;%y0(:,:,1);
[m,n]=size(x0);
% cb=y0(:,:,2);
% cr=y0(:,:,3);
k=mean(mean(x0));
XM=k; % Mean of luminense
r=round(k);
l=length(luma(luma<=r));
        listindex=find(luma<=r);
 k1=reshape(luma(listindex),l,1);
 sum=0;
 %for i=0:r
     xpdf=hist(k1,[0:r]);%pdf from 0:r
     xpdf=xpdf/l;%normalized pdf to get nk/n,l=sum of xpdf,total no of pixels form 0 to r.
%      plot(xpdf);
%      xlabel('gray levels up to mean');
%      ylabel('pdf up to mean');
%      title('histogram for half an image up to mean');
     sk=xpdf*triu(ones(r+1));
%      figure(2);
%      plot(sk);
%      xlabel('gray levels upto mean');
%      ylabel('cdf upto mean');
%      title('cdf for half  of an image up to mean');
     for k2=0:r
         if(xpdf(k2+1)>0)
             list1=find(k1==k2);%find value in an vector i.e converted from matrix
             list(list1)=sk(k2+1)*(r+1);%map
             ert(k2+1)=sk(k2+1)*(r+1);
         end
     end
%      for i=0:l-1
          p=zeros(m,n);
%          if (p(listindex))
%              p(:)=list;
%          end
%      end
             
     p(listindex)= list;
%     figure(3);
%     imshow(p);
%     xlabel('gray levels up to mean');
%      ylabel('luma component equilized image up to mean');
%      title('processed luma image up to mean');
lupper=length(luma(luma>r));
%for i=0:r
    %if(luma(luma<=r))
        listindexupper=find(luma>r);
   % end
%end
 k1upper=reshape(luma(listindexupper),lupper,1);
 sum=0;
 %for i=0:r
     xpdfupper=hist(k1upper,[r+1:255]);%pdf from r+1:255
     xpdfupper=xpdfupper/lupper;%normalized pdf to get nk/n,l=sum of xpdf,total no of pixels form r+1 to 255.
%      figure(4);
%      plot(xpdfupper);
%      xlabel('gray levels after mean');
%      ylabel('pdf after mean');
%      title('histogram for upper half an image after mean');
     skupper=xpdfupper*triu(ones(255-r));
%      figure(5);
%      plot(skupper);
%      xlabel('gray levels after mean');
%      ylabel('cdf after mean');
%      title('cdf for upper half  of an image after mean');
     k2u=1;
     for k2upper=(r+1):255
         %if(xpdfupper(k2upper)>r)
             list1upper=find(k1upper==k2upper);%find value in an vector i.e converted from matrix
             %for k2u=1:58
             listnew(list1upper)=(r+1+skupper(k2u)*(255-r));%map
             ert(k2upper)=(r+1+skupper(k2u)*(255-r));
         %end
         k2u=k2u+1;
     end
     
%      for i=0:l-1
          p1=zeros(m,n);
%          if (p(listindex))
%              p(:)=list;
%          end
%      end
             
     p1(listindexupper)= listnew;
%     figure(6);
%     imshow(uint8(p1));
%     xlabel('gray levels after mean');
%      ylabel('luma component equilized image after mean');
%      title('processed luma image after mean');
    ommmmm=p1+p;
    F_BHE_g=ommmmm;
   % figure(7);
%     imshow(ommmmm);
%      colormap('gray');
%       xlabel('gray level');
%      ylabel('combined lower and upper half luma component equilized image');
%      title('combined luma image');
%     figure(8);
%     imshow(uint8(ommmmm));
% % %     cat1=cat(3,ommmmm,cb,cr);
% % %     figure(8);
% % %     imshow(cat1);
% % %     xlabel('gray level(ycbcr)');
% % %      ylabel('combined lower and upper half luma,cromablue,croma red component equilized image');
% % %      title('luma croma b and r color processed image');
% % %     catconversion=ycbcr2rgb(cat1);
% % %     figure(9);
% % %     imshow(catconversion);
%     xlabel('gray level(rgb)');
%      ylabel('combined lower and upper half  RGB component equilized image');
%      title('converted from ycbcr2rgb color(RGB) processed image');
% YM=mean(mean(ommmmm));
%      AMBE=abs(XM-YM);
%      disp('Absolute mean brightness error=');
%      disp(AMBE);
%      MSE = MeanSquareError(luma, ommmmm);
%     % figure(16);
% disp('Mean Square Error = ');
% disp(MSE);
% PSNR = PeakSignaltoNoiseRatio(luma, ommmmm);
% %figure(17);
% disp('Peak Signal to Noise Ratio = ');
% disp(PSNR);
% E=entropy(uint8(ommmmm));
% disp('Entropy=');
% disp(E);
%      figure(16);
%      plot(imhist(ommmmm));
end
