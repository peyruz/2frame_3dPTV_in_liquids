% Calculate norms of N vectors containted in the input matrix.
% size(matrix)=[3,N]
% size(NM)=[1,N]

function [NM]=normMat(matrix)
        [NM]=sqrt(matrix(1,:).^2+matrix(2,:).^2+matrix(3,:).^2);
end
