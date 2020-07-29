function [c1Ref,c2Ref,dRF]=selfCalibOpt1(c1Est, c2Est,sz01,sz02,xi1r, xi2r,PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y,n,e3)
%% Adjustment of the dRF value through self-calibration with data.
% We minimize the total reprojection error by adjusting the dRF.
% 
% Input:
% c1Est, c2Est  - <1,3> current estimates of the camera positions, available from the calibtation
% sz01, sz02    - <2,N> interface images of the particles
% xi1r, xi2r    - <2,N> image coordinates of the particles
% PolyCoef1x, 
% PolyCoef1y, 
% PolyCoef2x, 
% PolyCoef2y    - <36,1> empirical calibration polynomial coefficients
% n             - <1,1> refractive index of the interface
% e3            - <1,3> basis vector
% 
% Output
% c1Ref, c2Ref  - <1,3> - refined positions of the cameras
% dRF           - <1,1> - refined dRF
% 
% by Peyruz Gasimov
% Last update: Nov 9, 2018

%% Initial guess for the dRF
xini=0;

%% Nelder-Mead optimizatoin
% g = @(x)selfCalibOptObjFun1(x,c1Est,c2Est,sz01,sz02,xi1r, xi2r,PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y,e3,n);
% options = optimset('Display','iter','MaxFunEvals',60000,'MaxIter',60000,'TolFun',1e-8,'TolX',1e-8);
% x = fminsearch(g,xini,options);

%% Trust-region reflective optimimization
optionsLsq = optimoptions(@lsqnonlin,'Display','iter',...
    'TypicalX', 0.001,...
    'MaxFunctionEvaluations',50000,'MaxIterations', 500,'StepTolerance',1e-5,'FunctionTolerance',1e-8,'UseParallel',true);

f = @(x)selfCalibOptObjFun1(x,c1Est,c2Est,sz01,sz02,xi1r, xi2r,PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y,e3,n);

% Define the search bounds
lb=-2;
ub=2;

[x] = lsqnonlin(f,xini,lb,ub,optionsLsq);

c1Ref=c1Est+[0; 0; x(1)];
c2Ref=c2Est+[0; 0; x(1)];
dRF=x(1);

%% Objective function for 1-variable optimization
function [ObjFunVal]=selfCalibOptObjFun1(x,c1Est,c2Est,sz01,sz02,xi1r, xi2r,PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y,e3,n)
dRF=x;
c1=c1Est+[0 0 dRF]';
c2=c2Est+[0 0 dRF]';

[~,ObjFunVal]=optimTriangulation4(PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y, sz01, sz02, xi1r, xi2r,c1,c2,n,e3,dRF);

end

end
