% This script deliniates the mismatch between the Calibration and the
% Experimental Reference Frames (dRF)

dRFini=0 % initial guess
    
g = @(dRF)dRFobjFun(dRF,sz01,sz02,c1p,c2p,n)
options = optimset('Display','iter','MaxFunEvals',60000,'MaxIter',60000,'TolFun',1e-10,'TolX',1e-10);
[dRF(k), fval] = fminsearch(g,dRFini,options)