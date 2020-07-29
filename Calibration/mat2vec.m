% Convert matrix data into a vertical vector
% Input
% mat - matrix <M,N>
% vec - vertical vector <M*N,1>

function [vec]=mat2vec(mat);

vec=mat(:);

end