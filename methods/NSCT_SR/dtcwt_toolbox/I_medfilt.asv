function I=I_medfilt(I1,R,N)

[r,c]=size(I1);
I=zeros(r,c);
se=strel('line',R,N);
[rr,cc]=size(getheight(se));
nhood = getnhood(se);
for i=(floor(rr/2)+1):(r-floor(rr/2))
    for j=(floor(cc/2)+1):(c-floor(cc/2))
       B=I1(i-floor(rr/2):i+floor(rr/2),j-floor(cc/2):j+floor(cc/2));
      index=find(nhood==1);
      S=B(index);
       %S=sort(S(:));
       I(i,j)=median(S); 
    end
end




%         I(i,:)=I(floor(rr/2)+1,:);
% 
% for i=floor(rr/2)+1:r
%     I(i,:)=I(r-floor(rr/2),:);
% end




