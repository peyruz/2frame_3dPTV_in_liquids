% This setup is based on the coordinate system of the actual 3dPTV cameras

% Define the setup - Tsai camera model

% Standard Basis for object space
e1=[1;0;0];
e2=[0;1;0];
e3=[0;0;1];

c1=[-70;-75;230];    % Positions of the camera. They are adjusted (by trial and error mostly) to provide the maximum overlap 
c2=[250;-75;230];

Hs=24.2; % mm, height of sensor
Ws=36.3; % mm, width of sensor

phi=degtorad(30);   % angle between the sensor and PoF

R1=[-cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 -cos(phi)];     % Transformation from sensor space into the object space for the first camera 
                                                            % (hint: for inverse transformation one can use the inverse (=transpose, in this case) of this matrix. Possible due to the orthogonality of the matrix; inv=tranp)
R2=[-cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 -cos(phi)];     % .. for the second camera. The normal for the sensor is directed away from the world origin (to keep the reference frame right-handed). 

% Standard Basis for sensor space

es11=R1(1:3,1);
es12=R1(1:3,2);
es13=R1(1:3,3);

es21=R2(1:3,1);
es22=R2(1:3,2);
es23=R2(1:3,3);


f=-45;     % mm. The distance from the pinhole to the sensor

% x0=[3296;2200]; % Principal point (pixels)

% Having the camera constants defined above from the position of the sensor
% we may calculate the position of the origin of the sensor coordinates

A1=c1+f*es13;    % A1 and A2 are principal points of the sensor (physical coord's)
A2=c2+f*es23;    

% Difference between principal point (PP) and image center (IC) ie t=IC-PP (to estimate in pixels: t/pixelSize)
% t=[-2 0]; % mm
t=[0 0];

oi1=A1+t(1)*es11+t(2)*es12-Ws/2*es11-Hs/2*es12; % Image origin point, physical space
oi2=A2+t(1)*es21+t(2)*es22-Ws/2*es21-Hs/2*es22;

ps=5.5e-3;         % Pixel size, mm

n=1.5;   % Refractive index

close all
PlotCameras;