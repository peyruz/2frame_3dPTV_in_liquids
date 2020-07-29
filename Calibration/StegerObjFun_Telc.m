function [ObjFunVal]=StegerObjFun_Telc(x,xw,xi)
% This function calculates the value of the objection function to fit
% Stegen Image-telecentric Camera model

%% Assign parameters
k1=x(1);     % Radial Distortion terms
k2=x(2);
k3=x(5);
tau1=x(3);   % Tangential distortion terms
tau2=x(4);

% A1=x(5);           % Affinity distortion 
A1=0;
A2=x(6);
B1=x(7);        % Shear distortion
B2=x(8);

q=quaternion(x(9),x(10),x(11),x(12));
t=[x(13); x(14); x(15)];

tiltAngle1=x(16);   % rho, in virtual plane, rad
tiltAngle2=x(17);   % tau, lens tilt angle, rad

fv=x(18);
% od=x(19);         % od=inf for telecentric lenses and os therefore
                    % irrelevant here

pxs=[x(20); x(21)];

distCen=[x(22); x(23)];

%% From the world RF side
xc=RotateVector(q,xw)+t;

pxuW(1,:)=xc(1,:)*fv./xc(3,:);       % Undistorted and unscaled image space
pxuW(2,:)=xc(2,:)*fv./xc(3,:);

%% From the Image RF side
% Tilt Rotation Matrix
A = [cos(tiltAngle1)^2*(1-cos(tiltAngle2))+cos(tiltAngle2),    cos(tiltAngle1)*sin(tiltAngle1)*(1-cos(tiltAngle2)),     sin(tiltAngle1)*sin(tiltAngle2);...
     cos(tiltAngle1)*sin(tiltAngle1)*(1-cos(tiltAngle2)),      sin(tiltAngle1)^2*(1-cos(tiltAngle2))+cos(tiltAngle2),  -cos(tiltAngle1)*sin(tiltAngle2);...
    -sin(tiltAngle1)*sin(tiltAngle2),                   cos(tiltAngle1)*sin(tiltAngle2),                  cos(tiltAngle2)];

% Telecentric camera projection matrix
Hinv=   [A(1,1) A(1,2) 0;
         A(2,1) A(2,2) 0;
         0      0      1];

pxdt=(xi-distCen).*pxs;     % sensor plane, local, distorted, unscaled

pxd=Hinv*[pxdt; ones(1, size(pxdt,2))];   % sensor
pxd=pxd./pxd(3,:);

% Find undistorted coordinates
r=(pxd(1,:).^2+pxd(2,:).^2).^(0.5);
pxuIm(1,:)=pxd(1,:).*(1+k1*r.^2+k2*r.^4+k3*r.^6)+tau1*(r.^2+2*pxd(1,:).^2)+2*tau2*pxd(1,:).*pxd(2,:)+A1*pxd(1,:)+B1*pxd(2,:);
pxuIm(2,:)=pxd(2,:).*(1+k1*r.^2+k2*r.^4+k3*r.^6)+tau2*(r.^2+2*pxd(2,:).^2)+2*tau1*pxd(1,:).*pxd(2,:)+A2*pxd(2,:)+B2*pxd(1,:);
%           ^radial distortion                     ^tangential distortion                            ^Affinity   ^Shear

%% Objective function calculation (for e.g. fminsearch)
% ObjFunVal=sum([(pxuIm(1,:)-pxuW(1,:)).^2, (pxuIm(2,:)-pxuW(2,:)).^2, (norm(q)-1)^2]);
   
% %% Objective function calculation (for e.g. fmincon)
ObjFunVal=[pxuIm(1,:)-pxuW(1,:), pxuIm(2,:)-pxuW(2,:), 1000*(norm(q)-1)];
   
   
   
   
   
   