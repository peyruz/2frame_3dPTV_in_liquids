% This setup is based on the coordinate system of the actual 3dPTV cameras
% clear
% clc
% Define the setup - Tsai camera model

% Standard Basis for object space
e1=[1;0;0];
e2=[0;1;0];
e3=[0;0;1];

c1=[-70;-75;230];    % Positions of the camera. They are adjusted (by trial and error mostly) to provide the maximum overlap 
c2=[250;-75;230];

Hs=24.2; % mm, height of sensor
Ws=36.3; % mm, width of sensor

phi1=degtorad(30);   % angle between the sensor and PoF
phi2=phi1;
tau1=deg2rad(4);     % tilt of lens
tau2=deg2rad(4);  

% tau=0;

Rs1=[-cos(phi1) 0 sin(phi1); 0 1 0; -sin(phi1) 0 -cos(phi1)];     % Transformation from sensor space into the object space for the first camera 
                                                            % (hint: for inverse transformation one can use the inverse (=transpose, in this case) of this matrix. Possible due to the orthogonality of the matrix; inv=tranp)
Rs2=[-cos(phi2) 0 -sin(phi2); 0 1 0; sin(phi2) 0 -cos(phi2)];     % .. for the second camera. The normal for the sensor is directed away from the world origin (to keep the reference frame right-handed). 

Rv1=vecAngle2RotMatrix([0 1 0],pi-phi1+tau1);
Rv2=vecAngle2RotMatrix([0 1 0],pi+phi2-tau2);


% Standard Basis for sensor space
es11=Rs1(1:3,1);
es12=Rs1(1:3,2);
es13=Rs1(1:3,3);

es21=Rs2(1:3,1);
es22=Rs2(1:3,2);
es23=Rs2(1:3,3);

% Standard Basis for virtual image space
ev11=Rv1(1:3,1);
ev12=Rv1(1:3,2);
ev13=Rv1(1:3,3);

ev21=Rv2(1:3,1);
ev22=Rv2(1:3,2);
ev23=Rv2(1:3,3);


f1=-45;     % mm. The distance along sensor z-axis from the pinhole to the sensor
f2=-45;

od1=f1/cos(tau1); % mm. The distance along the optical axis to sensor from the projection center
od2=f2/cos(tau2);

fv1=-20;    % mm. The virtual principal distance (ie distance to the virtual image plane form the projection center)
fv2=-20;

% fv=od;    % pinhole model

% Having the camera constants defined above from the position of the sensor
% we may calculate the position of the origin of the sensor coordinates

As1=c1+f1*es13;    % A1 and A2 are principal points of the sensor (world RF)
As2=c2+f2*es23;    

% Difference between principal point (PP) and image center (IC) ie t=IC-PP (to estimate in pixels: t/pixelSize)
% t=[-2 0]; % mm
% t=[-0.5 0];
t=[0 0];

oi1=As1+t(1)*es11+t(2)*es12-Ws/2*es11-Hs/2*es12; % Image origin point, physical space
oi2=As2+t(1)*es21+t(2)*es22-Ws/2*es21-Hs/2*es22;


% ps=5.5e-3;         % Pixel size, mm

n=1.5;   % Refractive index

close all
PlotCameras5;