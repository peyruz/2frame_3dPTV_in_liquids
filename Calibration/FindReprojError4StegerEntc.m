% Calculate reprojection error of a calibration

% Input
% xw -      <3,N> physical coordinates of the calibration dots, mm
% CamCal -  <1,20> vector containing the calibration parameters of the camera
% xi -      <2,N> image coordinates of the calibration dots, px

% Output
% ReprojError - Reprojection Error for each input dot

% Example usage:
% Output only norm of the Reprojection Error
%   [ReprojErrorNorm]=FindReprojError(xw,CamCal,xi)
% Output additionally the actual Reprojection Error values.
%   [ReprojErrorNorm,ReprojError]=FindReprojError(xw,CamCal,xi)

function [ReprojErrorNorm,varargout]=FindReprojError4StegerEntc(xw,CamCal,xi)

if size(xw,1)>size(xw,2)
    xw=xw';
end

if size(xi,1)>size(xi,2)
    xi=xi';
end

if size(xw,2)~=size(xi,2)
    error('The number of input physical points and image points needs to be the same.')
end
%% Assign parameters
k1=CamCal(1);     % Radial Distortion terms
k2=CamCal(2);

tau1=CamCal(3);   % Tangential distortion terms
tau2=CamCal(4);

A1=CamCal(5);           % Affinity distortion 
A2=CamCal(6);
B1=CamCal(7);        % Shear distortion
B2=CamCal(8);

q=quaternion(CamCal(9),CamCal(10),CamCal(11),CamCal(12));
t=[CamCal(13); CamCal(14); CamCal(15)];

tiltAngle1=CamCal(16);   % rho, in virtual plane, rad
tiltAngle2=CamCal(17);   % tau, lens tilt angle, rad

fv=CamCal(18);
od=CamCal(19);         
                    

pxs=[CamCal(20); CamCal(21)];

distCen=[CamCal(22); CamCal(23)];

%% From the world RF side
xc=RotateVector(q,xw)+t;

pxuW(1,:)=xc(1,:)*fv./xc(3,:);       % Undistorted and unscaled image space
pxuW(2,:)=xc(2,:)*fv./xc(3,:);

%% From the Image RF side
% Tilt Rotation Matrix
A = [cos(tiltAngle1)^2*(1-cos(tiltAngle2))+cos(tiltAngle2),    cos(tiltAngle1)*sin(tiltAngle1)*(1-cos(tiltAngle2)),     sin(tiltAngle1)*sin(tiltAngle2);...
     cos(tiltAngle1)*sin(tiltAngle1)*(1-cos(tiltAngle2)),      sin(tiltAngle1)^2*(1-cos(tiltAngle2))+cos(tiltAngle2),  -cos(tiltAngle1)*sin(tiltAngle2);...
    -sin(tiltAngle1)*sin(tiltAngle2),                   cos(tiltAngle1)*sin(tiltAngle2),                  cos(tiltAngle2)];

% Perspective camera projection matrix
H = [A(1,1)*A(3,3)-A(1,3)*A(3,1), A(2,1)*A(3,3)-A(2,3)*A(3,1),    0;...
     A(1,2)*A(3,3)-A(1,3)*A(3,2), A(2,2)*A(3,3)-A(2,3)*A(3,2),    0;...
     A(1,3)/od,                   A(2,3)/od,                      A(3,3)];

pxdt=(xi-distCen).*pxs;     % sensor plane, local, distorted, unscaled

pxd=H\[pxdt; ones(1, size(pxdt,2))];   % sensor
pxd=pxd./pxd(3,:);

% Find undistorted coordinates
r=(pxd(1,:).^2+pxd(2,:).^2).^(0.5);
pxuIm(1,:)=pxd(1,:).*(1+k1*r.^2+k2*r.^4)+tau1*(r.^2+2*pxd(1,:).^2)+2*tau2*pxd(1,:).*pxd(2,:)+A1*pxd(1,:)+B1*pxd(2,:);
pxuIm(2,:)=pxd(2,:).*(1+k1*r.^2+k2*r.^4)+tau2*(r.^2+2*pxd(2,:).^2)+2*tau1*pxd(1,:).*pxd(2,:)+A2*pxd(2,:)+B2*pxd(1,:);
%           ^radial distortion          ^tangential distortion                             ^Affinity   ^Shear

% Scale to pixels
pxuImpx=pxs.^(-1).*pxuIm;
pxuWpx=pxs.^(-1).*pxuW;

%% Calculate the error
ReprojErrorNorm=sqrt((pxuImpx(1,:)-pxuWpx(1,:)).^2+(pxuImpx(2,:)-pxuWpx(2,:)).^2);
RepErrRMS=sqrt(mean((pxuImpx(1,:)-pxuWpx(1,:)).^2+(pxuImpx(2,:)-pxuWpx(2,:)).^2));
fprintf('Mean Reprojection Error: %2.4f px \n',mean(ReprojErrorNorm));
fprintf('Reprojection Error RMS: %2.4f px \n',RepErrRMS);

% Variable output
if nargout==2
   varargout{1}=[pxuWpx(1,:); pxuWpx(2,:); pxuImpx(1,:)-pxuWpx(1,:); pxuImpx(2,:)-pxuWpx(2,:)];
end

% Visualization
tri=delaunay(pxuImpx(1,:),pxuImpx(2,:));
RepErrorFig3(pxuImpx(1,:)+distCen(1),pxuImpx(2,:)+distCen(2),pxuImpx(1,:)-pxuWpx(1,:),pxuImpx(2,:)-pxuWpx(2,:), distCen(1),distCen(2),800,...
                   tri, ReprojErrorNorm, ReprojErrorNorm);