% part of 2d Empirical calibration of 3dPTV package
% Construct the matrix for Pseudo-Inverse regression
% 
% Inverse mapping:
% Image RF -> World RF (z=0 plane)
% 
% Input
% ix, iy - 1d vectors - image coordinates of points, px
%
% Output
% PolyCalSystem - matrix based on 7th order polynomials
%
% --- Changelog
% version 2, Jun 5, 2018
% Improved speed through vectorization (especially for larger arrays)
% 
% by Peyruz Gasimov Sep, 2017 - Jun 5, 2018

function [invPolyCalSystem]=invPolyCal2d7ver2(ix,iy)

% Preallocate memory
invPolyCalSystem=zeros(length(ix),36);

% Convert the pixels into mm for a better conditioning
ix=(ix-3300)*5.5e-3;   
iy=(iy-2200)*5.5e-3;

invPolyCalSystem(:,1)=1;
invPolyCalSystem(:,2)=ix.^4;
invPolyCalSystem(:,3)=iy.^4;
invPolyCalSystem(:,4)=ix.^3.*iy;
invPolyCalSystem(:,5)=iy.^3.*ix;
invPolyCalSystem(:,6)=ix.^2.*iy.^2;
invPolyCalSystem(:,7)=ix.^3;
invPolyCalSystem(:,8)=iy.^3;
invPolyCalSystem(:,9)=ix.^2.*iy;
invPolyCalSystem(:,10)=iy.^2.*ix;
invPolyCalSystem(:,11)=ix.^2;
invPolyCalSystem(:,12)=iy.^2;
invPolyCalSystem(:,13)=ix.*iy;
invPolyCalSystem(:,14)=ix;
invPolyCalSystem(:,15)=iy;
invPolyCalSystem(:,16)=ix.^5;
invPolyCalSystem(:,17)=iy.^5;
invPolyCalSystem(:,18)=ix.^4.*iy;
invPolyCalSystem(:,19)=iy.^4.*ix;
invPolyCalSystem(:,20)=ix.^3.*iy.^2;
invPolyCalSystem(:,21)=iy.^3.*ix.^2;
invPolyCalSystem(:,22)=ix.^6;
invPolyCalSystem(:,23)=iy.^6;
invPolyCalSystem(:,24)=ix.^5.*iy;
invPolyCalSystem(:,25)=iy.^5.*ix;
invPolyCalSystem(:,26)=ix.^4.*iy.^2;
invPolyCalSystem(:,27)=iy.^4.*ix.^2;
invPolyCalSystem(:,28)=iy.^3.*ix.^3;
invPolyCalSystem(:,29)=ix.^7;
invPolyCalSystem(:,30)=iy.^7;
invPolyCalSystem(:,31)=ix.^6.*iy;
invPolyCalSystem(:,32)=iy.^6.*ix;
invPolyCalSystem(:,33)=ix.^5.*iy.^2;
invPolyCalSystem(:,34)=iy.^5.*ix.^2;
invPolyCalSystem(:,35)=ix.^4.*iy.^3;
invPolyCalSystem(:,36)=iy.^4.*ix.^3;


