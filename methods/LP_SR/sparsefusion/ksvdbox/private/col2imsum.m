function a = col2imsum(b,sz,imgsz,delta)
%COL2IMSUM Rearrange matrix columns as overlapping blocks with summation
%  A = COL2IMSUM(B,[M N],[MM NN]) rearranges the columns of B into sliding
%  M-by-N blocks, producing the matrix A of size MM-by-NN. B is usually the
%  result of calling IM2COL(...,'sliding'), and must have (MM-M+1)*(NN-N+1)
%  columns. Overlapping blocks in A are summed.
%
%  A = COL2IMSUM(B,[M N],[MM NN],[D1 D2]) places the blocks in A with a
%  distance of [D1 D2] between them. If not specified, the default distance
%  is [1 1]. A value of [M N] implies no overlap.
%
%  See also COL2IM.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009

if (nargin<4), delta = [1 1]; end

a = zeros(imgsz);

k = 1;
for j = 1:sz(2)
  for i = 1:sz(1)
    ids_i = i+(0:imgsz(1)/delta(1)-sz(1)/delta(1))*delta(1);
    ids_j = j+(0:imgsz(2)/delta(2)-sz(2)/delta(2))*delta(2);
    a(ids_i,ids_j) = a(ids_i,ids_j) + reshape(b(k,:),[length(ids_i) length(ids_j)]);
    k = k+1;
  end
end
