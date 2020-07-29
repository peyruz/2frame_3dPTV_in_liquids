% Calculate norms of N vectors containted in the input matrix.
% size(matrix)=[2,N]
% size(NM)=[1,N]

function [NM]=normMat2d(matrix)
        [NM]=sqrt(matrix(1,:).^2+matrix(2,:).^2);     
end
