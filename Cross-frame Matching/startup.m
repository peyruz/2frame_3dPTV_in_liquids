clear
clc
N1=500;
unmatched=50;
N2=N1+unmatched;
pos1=rand(N1,3);
pos1(:,3)=0;
pos2=pos1+repmat([0.02, 0.025, 0],size(pos1,1),1);
pos2=[pos2; [rand(unmatched,2) zeros(unmatched,1)]];
ix=randperm(N2,N2);
pos2=pos2(ix,1:3);

Rn=0.04;
Rs=0.05;
Rq=0.1*Rs;
A=0.3;
B=3;
[NeighborMatrix]=BuildNeighborMatrix(pos1, Rn);
[MatchMatrix, NoMatchP]=BuildInitialMatchMatrix(pos1, pos2, Rs);