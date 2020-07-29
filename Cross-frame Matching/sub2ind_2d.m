% input
% dim1 - first dimension of the array (number of rows)
% v1 - row index (vector :,1)
% v2 - column index (vector :,1)
% size(v1,1)==size(v2,1) = 1

function [lindex]=sub2ind_2d(dim1, v1, v2)

lindex=repmat(dim1,size(v1,1),1).*(v2-1)+v1;

end