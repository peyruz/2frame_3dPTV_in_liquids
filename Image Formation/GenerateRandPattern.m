% Script for generating pattern of random particle positions

function s=GenerateRandPattern(xmin,xmax,ymin,ymax,zmin,zmax,N);

% N=1000  % Total number

% xmin=110
% xmax=300
%     
% ymin=-150
% ymax=-50
% 
% zmin=-2
% zmax=-1

s=[xmin+rand(1,N)*(xmax-xmin)
   ymin+rand(1,N)*(ymax-ymin)
   zmin+rand(1,N)*(zmax-zmin)];