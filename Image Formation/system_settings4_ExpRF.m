% Generate a the system in an Experimental RF


c1ExpRF=[-70;-75;230+e_dRF(k)];    % Positions of the camera with the misalignment of experimental and calibration reference frames (e_dRF) included
c2ExpRF=[250;-75;230+e_dRF(k)];

A1ExpRF=c1ExpRF+f*es13;   % A1 and A2 are principal points of the sensor (physical coord's)
A2ExpRF=c2ExpRF+f*es23;    

oi1ExpRF=A1ExpRF+t(1)*es12+t(2)*es11-Ws/2*es11-Hs/2*es12; % Image origin point, physical space
oi2ExpRF=A2ExpRF-t(1)*es21+t(2)*es22-Ws/2*es21-Hs/2*es22;
