% Calculate unquantified distortion in y for Camera 1
% Based on the real calibration
% x is normalized by mean 4000 and std 1393
% y is normalized by mean 2042 and std 1022

function [ydist]=distortInYCam1(x,y)

x=(x-4000)/1393;
y=(y-2042)/1022;

p00 =    -0.02516;
p10 =     0.06445;
p01 =     0.04938;
p20 =     0.06253;
p11 =     0.09046;
p02 =    0.001029;
p30 =    -0.02914;
p21 =     -0.1936;
p12 =    -0.02798;
p03 =     0.02946;
p40 =    -0.01701;
p31 =   -0.006681;
p22 =    -0.02811;
p13 =    -0.02492;
p04 =    0.007845;
p50 =    0.005477;
p41 =     0.03481;
p32 =   -0.006797;
p23 =     0.03211;
p14 =  -0.0001477;
p05 =   -0.009599;

ydist = p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 + p21.*x.^2.*y ...
        + p12.*x.*y.^2 + p03.*y.^3 + p40.*x.^4 + p31.*x.^3.*y + p22.*x.^2.*y.^2 ...
        + p13.*x.*y.^3 + p04.*y.^4 + p50.*x.^5 + p41.*x.^4.*y + p32.*x.^3.*y.^2 ...
        + p23.*x.^2.*y.^3 + p14.*x.*y.^4 + p05.*y.^5;

