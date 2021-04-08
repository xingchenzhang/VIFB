function y=nsct_fuse(I1,I2,nlevels)
%    NSCT
%    Input:
%    I1 - input image A
%    I2 - input image B
%    nlevels - number of directions in each decomposition level
%    Output:
%    y  - fused image   
%
%    The code is edited by Yu Liu, 01-09-2014.


%nlevels = [2,3,3,4] ;       
pfilter = 'pyrexc' ;              
dfilter = 'vk' ; 
 
I1=double(I1);
I2=double(I2);
coeffs_1 = nsctdec( I1, nlevels, dfilter, pfilter );
coeffs_2 = nsctdec( I2, nlevels, dfilter, pfilter );


[m,n]=size(I1);
coeffs=coeffs_2;
for i=2:numel(nlevels)+1

    if nlevels(i-1)==0
        E1=abs(coeffs_1{i});
        E2=abs(coeffs_2{i});
        %             map=E1>E2;
         um=3;
    A1 = ordfilt2(abs(es2(E1,floor(um/2))), um*um, ones(um));
  	A2 = ordfilt2(abs(es2(E2,floor(um/2))), um*um, ones(um));
    % second step
  	map= (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
        coeffs{i}(map)=coeffs_1{i}(map);
    else
        for j=1:(2^nlevels(i-1))
            E1=abs(coeffs_1{i}{j});
            E2=abs(coeffs_2{i}{j});
%             map=E1>E2;
    um=3;
    A1 = ordfilt2(abs(es2(E1,floor(um/2))), um*um, ones(um));
  	A2 = ordfilt2(abs(es2(E2,floor(um/2))), um*um, ones(um));
    % second step
  	map= (conv2(double(A1 > A2), ones(um), 'valid')) > floor(um*um/2);
            coeffs{i}{j}(map)=coeffs_1{i}{j}(map);
        end
    end
end

% energy based method

coeffs{1}=(coeffs_1{1}+coeffs_2{1})/2;

y= nsctrec( coeffs, dfilter, pfilter ) ;


