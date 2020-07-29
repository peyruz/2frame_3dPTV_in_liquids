function [s]=GenerateRecPattern(xmin,xmax,xn,ymin,ymax,yn,depth)
% 
% xmin=110
% xmax=300
% xn=20
% 
% ymin=-150
% ymax=-50
% yn=10

sx=linspace(xmin,xmax,xn);
sy=linspace(ymax,ymin,yn);

[sxx syy]=meshgrid(sx,sy);

sx=sxx(:);
sy=syy(:);

sz=depth*ones(1,size(sx,1));

sx=sx';
sy=sy';

s=[sx;sy;sz];
