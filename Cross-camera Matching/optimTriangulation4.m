% Optimal triangulation based on minimization of reprojection errors
% The algorithm is based on the idea stated in Kanatani on optimal
% triangulation via minimization of reprojection errors. 
%
% Input:
% PolyCoef1x,PolyCoef1y,
% PolyCoef2x,PolyCoef2y -   Empirical Calibration polynomial coefficients
% sz01, sz02            -   Interface images of the detected particles
% xi1, xi2              -   Image coordinates of the detected particles
% c1, c2                -   Camera positions
% n                     -   Refractive index
% e3                    -   3rd basis vector
% dRF                   -   z-difference between the calibration and
%                           experimental reference frames
% 
% Output:
% x3All                 -   3d position of the particles obtained via the
%                           optimal triangulation
% SumSqRepErr           -   Sum of 
% Take dRF into account

% Changelog:
% ver 4, Jan, 2018
% - downhill simplex optimization algorithm is replaced by the
%   trust-region-reflective (dramatic speedup)

function [x3dAll,cumResnorm]=optimTriangulation4(PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y, sz01, sz02,xi1,xi2,c1,c2,n,e3,dRF)

% Initial Guess (calculated in the Experimental RF)
[x3dini]=TriangulateFromDewarpedImageVec(sz01,sz02,c1,c2,e3,n,dRF);

% Slice the variables for pararallization
x3diniX=x3dini(1,:);
x3diniY=x3dini(2,:);
x3diniZ=x3dini(3,:);

xi1X=xi1(1,:);
xi1Y=xi1(2,:);
xi2X=xi2(1,:);
xi2Y=xi2(2,:);

cumResnorm=0;

parfor ii=1:size(sz01,2)
    xi1cur=[xi1X(ii);xi1Y(ii)];
    xi2cur=[xi2X(ii);xi2Y(ii)];
    
    x3dinicur=[x3diniX(ii);x3diniY(ii);x3diniZ(ii)];
    
    optionsLsq = optimoptions(@lsqnonlin,'Display','none',...
        'TypicalX', x3dinicur,...
        'MaxFunctionEvaluations',50000,'MaxIterations', 500,'StepTolerance',1e-10,'FunctionTolerance',1e-8,'UseParallel',false);
    
    f = @(x3d)optimTrianObjFun2(x3d, c1, c2, xi1cur, xi2cur, n, PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y,dRF);
    
    searchMargin=1;
    lb=[x3dinicur(1)-searchMargin,x3dinicur(2)-searchMargin,-4.5];
    ub=[x3dinicur(1)+searchMargin,x3dinicur(2)+searchMargin,-3];
    
    [x3d,resnorm] = lsqnonlin(f,x3dinicur,lb,ub,optionsLsq);
    
    x3dX(ii)=x3d(1);
    x3dY(ii)=x3d(2);
    x3dZ(ii)=x3d(3);
    
    cumResnorm=cumResnorm+resnorm;
end

x3dAll=[x3dX;x3dY;x3dZ];

end

%% Objective function
function [ObjFunVal]=optimTrianObjFun2(x3d, c1, c2, xi1, xi2, n, PolyCoef1x,PolyCoef1y,PolyCoef2x,PolyCoef2y, dRF)

% Convert the camera coordinates to the Experimental RF
c1(3)=c1(3)+dRF;
c2(3)=c2(3)+dRF;

%% Solve inverse refraction problem
% Camera 1
ini1=x3d-c1;
ini1=ini1/norm(ini1);   % initial guess for u_i1
ini1=ini1(2);

solx1 = fzero(@(i12)InvRefract(i12,c1,x3d(1),x3d(2),x3d(3),n), ini1);
u_i12=real(double(solx1));
u_i11=u_i12*(x3d(1)-c1(1))/(x3d(2)-c1(2));
u_i13=-sqrt(1-u_i11^2-u_i12^2);

u_i1=[u_i11; u_i12; u_i13];

sz01=c1+u_i1*abs((c1(3)-dRF)/u_i1(3)); 
sz01=sz01(1:2,:);

% Camera 2
ini2=x3d-c2;
ini2=ini2/norm(ini2);   % initial guess for u_i2
ini2=ini2(2);

solx2 = fzero(@(i22)InvRefract(i22,c2,x3d(1),x3d(2),x3d(3),n), ini2);
u_i22=real(double(solx2));
u_i21=u_i22*(x3d(1)-c2(1))/(x3d(2)-c2(2));
u_i23=-sqrt(1-u_i21^2-u_i22^2);

u_i2=[u_i21; u_i22; u_i23];

sz02=c2+u_i2*abs((c2(3)-dRF)/u_i2(3));
sz02=sz02(1:2,:);

%% Calculate Objective Function
PolyCalMatrix1=PolyCal2d7ver2(sz01(1),sz01(2));
PolyCalMatrix2=PolyCal2d7ver2(sz02(1),sz02(2));

ObjFunVal=  [PolyCalMatrix1*PolyCoef1x-xi1(1) PolyCalMatrix1*PolyCoef1y-xi1(2),...
             PolyCalMatrix2*PolyCoef2x-xi2(1) PolyCalMatrix2*PolyCoef2y-xi2(2)];

end
