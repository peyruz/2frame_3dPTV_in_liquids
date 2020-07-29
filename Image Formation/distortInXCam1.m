% Calculate unquantified distortion in x for Camera 1
% Based on the real calibration
% x is normalized by mean 4000 and std 1393
% y is normalized by mean 2042 and std 1022

function [xdist]=distortInXCam1(x,y)

x=(x-4000)/1393;
y=(y-2042)/1022;

p00 =     -0.0392;
p10 =      0.2954;
p01 =     0.03333;
p20 =    -0.02545;
p11 =     0.03167;
p02 =     0.05814;
p30 =     -0.3642;
p21 =    -0.02949;
p12 =    -0.06227;
p03 =  -9.712e-05;
p40 =     0.01027;
p31 =    -0.01558;
p22 =    0.005552;
p13 =   -0.008111;
p04 =   -0.002297;
p50 =     0.07929;
p41 =   -0.004498;
p32 =     0.05784;
p23 =   -0.005178;
p14 =    0.007905;
p05 =    0.003011;

xdist = p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 + p21.*x.^2.*y ...
        + p12.*x.*y.^2 + p03.*y.^3 + p40.*x.^4 + p31.*x.^3.*y + p22.*x.^2.*y.^2 ...
        + p13.*x.*y.^3 + p04.*y.^4 + p50.*x.^5 + p41.*x.^4.*y + p32.*x.^3.*y.^2 ...
        + p23.*x.^2.*y.^3 + p14.*x.*y.^4 + p05.*y.^5;


