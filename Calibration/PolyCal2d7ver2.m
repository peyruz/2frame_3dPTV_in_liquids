% part of 2d Empirical calibration of 3dPTV package
% Construct the matrix for Pseudo-Inverse regression
% 
% Planar (z=0) World RF -> Image RF
% 
% Input
% xw, yw - 1d vectors - world coordinates of points, mm
%
% Output
% PolyCalSystem - matrix based on 7th order polynomials
%
% --- Changelog
% version 2, Jun 5, 2018
% Improved speed through vectorization (especially for larger arrays)
% 
% by Peyruz Gasimov, Jun 5, 2018

function [PolyCalSystem]=PolyCal2d7ver2(xw,yw)

% Preallocate memory
PolyCalSystem=zeros(size(xw,1),36);

% Scale and shift the coordinates for better conditioning
xw=(xw-100)*1e-2;
yw=(yw+100)*1e-2;

PolyCalSystem(:,1)=1;
PolyCalSystem(:,2)=xw.^4;
PolyCalSystem(:,3)=yw.^4;
PolyCalSystem(:,4)=xw.^3.*yw;
PolyCalSystem(:,5)=yw.^3.*xw;
PolyCalSystem(:,6)=xw.^2.*yw.^2;
PolyCalSystem(:,7)=xw.^3;
PolyCalSystem(:,8)=yw.^3;
PolyCalSystem(:,9)=xw.^2.*yw;
PolyCalSystem(:,10)=yw.^2.*xw;
PolyCalSystem(:,11)=xw.^2;
PolyCalSystem(:,12)=yw.^2;
PolyCalSystem(:,13)=xw.*yw;
PolyCalSystem(:,14)=xw;
PolyCalSystem(:,15)=yw;
PolyCalSystem(:,16)=xw.^5;
PolyCalSystem(:,17)=yw.^5;
PolyCalSystem(:,18)=xw.^4.*yw;
PolyCalSystem(:,19)=yw.^4.*xw;
PolyCalSystem(:,20)=xw.^3.*yw.^2;
PolyCalSystem(:,21)=yw.^3.*xw.^2;
PolyCalSystem(:,22)=xw.^6;
PolyCalSystem(:,23)=yw.^6;
PolyCalSystem(:,24)=xw.^5.*yw;
PolyCalSystem(:,25)=yw.^5.*xw;
PolyCalSystem(:,26)=xw.^4.*yw.^2;
PolyCalSystem(:,27)=yw.^4.*xw.^2;
PolyCalSystem(:,28)=yw.^3.*xw.^3;
PolyCalSystem(:,29)=xw.^7;
PolyCalSystem(:,30)=yw.^7;
PolyCalSystem(:,31)=xw.^6.*yw;
PolyCalSystem(:,32)=yw.^6.*xw;
PolyCalSystem(:,33)=xw.^5.*yw.^2;
PolyCalSystem(:,34)=yw.^5.*xw.^2;
PolyCalSystem(:,35)=xw.^4.*yw.^3;
PolyCalSystem(:,36)=yw.^4.*xw.^3;
