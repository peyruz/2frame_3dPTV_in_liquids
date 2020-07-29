% Script to execute Monte Carlo simulations of 3d-PTV measurements.
% This version updates the calibration procedure
clc


NoT=3300
start_k=3001; % Number of Monte-Carlo trials
for k=start_k:NoT
clc
clearvars -except Error k NoT e_ls e_ds e_um c1ERR c2ERR e_dRF ErrorCORR dRF start_k
    
% We start off by creating the environment (Calibration Reference Frame)
system_settings4;

% Create a rectangular array of dots

xmin=-25;   % Set up the dimensions of the calibration pattern
xmax=170;
ymin=-125;
ymax=-15;
dotNx=70;   % Number of dots in x/y
dotNy=50;

e_ls(k)=normrnd(0, deg2rad(0.05));  % Level sensor error
% e_ls(k)=0

f1 = @(x1)RotMinmz(x1,e_ls(k));
xini=[1 0 0 0];
options = optimset('Display','iter','MaxFunEvals',60000,'MaxIter',60000,'TolFun',1e-10,'TolX',1e-10);
[x1, fval] = fminsearch(f1,xini,options);
q_ls=quaternion(x1);

dz=[1,-9];   % Distance between calibration planes

e_ds(k)=normrnd(0,0.001);  % Error of displacement sensor
% e_ds(k)=0

[s2c]=GenerateRecPattern(xmin,xmax,dotNx,ymin,ymax,dotNy,0);
[s1c]=s2c+repmat(RotateVector(q_ls,[0,0,dz(1)+e_ds(k)]'),1,size(s2c,2));    
[s3c]=s2c+repmat(RotateVector(q_ls,[0,0,dz(2)+e_ds(k)]'),1,size(s2c,2));    % This imitates translation up and down of the calibration target. 
                                                                        % Ideally the translation is normal to the calibration target however in reality 
                                                                        % we need to take into account the error of the level sensor.

sc=[s1c,s2c,s3c];

% Capture calibration images
RInt=0;     % Rint - boolean to denote presense of refracting interface at z=0 plane(=1 if present)
[xi1c, xi2c]=Image_Formation_ver4(c1, c2, A1,A2, R1, R2, es13,es23, sc,oi1,oi2,RInt,n);
[xi1s2c, xi2s2c]=Image_Formation_ver4(c1, c2, A1,A2, R1, R2, es13,es23, s2c,oi1,oi2,RInt,n);

    
%% Physics-based calibration
% Camera 1
% Initial guesses
pxsyini=5.5e-3;
pxsxini=5.5e-3;

% It is a good idea to not simply assume the principal point to be in the
% middle of the sensor. Rather, we know that because of the tilt the
% principal point will be shifted. In ths particular experiment, we can
% assume that the point will be shifted towards the origin in x.
x0ini1=[3100;2200]; 

% Initial estimates
[fini1,skini1,tini1,qini1]=TsaiLin(xi1c,sc,pxsxini, pxsyini,x0ini1);
qvini1=[real(qini1) (vector(qini1))'];

xini1=[0,0,0,0,qvini1(1),qvini1(2),qvini1(3),qvini1(4),tini1(1), tini1(2), tini1(3),fini1,pxsxini,pxsyini,x0ini1(1), x0ini1(2), 0 0 0 0];

% Non-linear optimization
f1 = @(x1)TsaiObjFunVec_v1_1(x1,sc,xi1c);
lb1=[-0.001 -0.00001 -0.01 -0.01 -1 -1 -1 -1 -200 -200 0 -55 0.0054 0.0054 2000 1000 -0.1 -0.1 -0.1 -0.1];
ub1=[0.001 0.00001 0.01 0.01 1 1 1 1 200 200 400 -40 0.0056 0.0056 4000 3000 0.1 0.1 0.1 0.1];

% Nelder-Mead method
% optionsNM = optimset('Display','iter','MaxFunEvals',400000,'MaxIter',200000,'TolFun',1e-10,'TolX',1e-10);
% [x1, fval1] = fminsearchbnd(f1,xini1,lb1,ub1,optionsNM);

% Trust-region-reflective algorithm
optionsLsq = optimoptions('lsqnonlin','Display','iter',...
    'TypicalX', [1e-6 1e-8, 1e-8 1e-8 0.5 0.5 0.5 0.5 10 10 100 50 0.005 0.005  3000 2000 0.01 0.01 1e-4 1e-4],...
    'MaxFunctionEvaluations',50000,'MaxIterations', 1500,'StepTolerance',1e-8);
x1 = lsqnonlin(f1,xini1,lb1,ub1,optionsLsq);

q1=quaternion([x1(5) x1(6) x1(7) x1(8)]);
t1=[x1(9);x1(10);x1(11)];
c1p=-RotateVector(inverse(q1),t1);

% % No distortion case (otherwise - uncomment)
% x1=[0 0 0 0 x1 0 0 0 0];

% Calculate reprojection error
[ReprojErrorNorm1]=FindReprojError(s2c,x1,xi1s2c);

%%
% Camera 2 
% For the second camera, we expect the principal point to have, on
% the contrary, higher x value due to the tilt.
x0ini2=[3350;2200];

% Initial estimates
[fini2,skini2,tini2,qini2]=TsaiLin(xi2c,sc,pxsxini, pxsyini,x0ini2);
qvini2=[real(qini2) (vector(qini2))'];

xini2=[0,0,0,0,qvini2(1),qvini2(2),qvini2(3),qvini2(4),tini2(1), tini2(2), tini2(3),fini2,pxsxini,pxsyini,x0ini2(1), x0ini2(2),0,0,0,0];

% Nelder-Mead optimization
f2 = @(x2)TsaiObjFunVec_v1_1(x2,sc,xi2c);
lb2=[-1 -1 -1 -1 -inf -inf -inf -inf -200 -200 0 -55 0.0054 0.0054 3000 1000 -1 -1 -1 -1];
ub2=[0.001 0.001 0.001 0.001 inf inf inf inf 200 200 400 -40 0.0056 0.0056 5000 3000 1 1 1 1];

[x2, fval2] = fminsearchbnd(f2,xini2,lb2,ub2,options);

q2=quaternion([x2(5) x2(6) x2(7) x2(8)]);
t2=[x2(9);x2(10);x2(11)];
c2p=-RotateVector(inverse(q2),t2);

% % No distortion case (otherwise - uncomment)
% x2=[0 0 0 0 x2 0 0 0 0];

% Calculate reprojection error
[ReprojErrorNorm2]=FindReprojError(s2c,x2,xi2s2c);

%% Inverse Empirical Calibration (Image RF -> World RF)
% Pseudo-inverse linear regression
% Camera 1
[invPolyCalSystem1c]=invPolyCal2d7(xi1s2c(1,:),xi1s2c(2,:));
invPolyCoef1x=(invPolyCalSystem1c\eye(size(invPolyCalSystem1c,1)))*s2c(1,:)';
invPolyCoef1y=(invPolyCalSystem1c\eye(size(invPolyCalSystem1c,1)))*s2c(2,:)';
FindEmpInvReprojError(s2c',invPolyCalSystem1c,invPolyCoef1x,invPolyCoef1y);

% Camera 2
[invPolyCalSystem2c]=invPolyCal2d7(xi2s2c(1,:),xi2s2c(2,:));
invPolyCoef2x=(invPolyCalSystem2c\eye(size(invPolyCalSystem2c,1)))*s2c(1,:)';
invPolyCoef2y=(invPolyCalSystem2c\eye(size(invPolyCalSystem2c,1)))*s2c(2,:)';
FindEmpInvReprojError(s2c',invPolyCalSystem2c,invPolyCoef2x,invPolyCoef2y);

%% Direct Empirical Calibration (World RF -> Image RF)
% Pseudo-inverse linear regression
% Camera 1
[PolyCalSystem1c]=PolyCal2d7(s2c(1,:),s2c(2,:));
PolyCoef1x=(PolyCalSystem1c\eye(size(PolyCalSystem1c,1)))*xi1s2c(1,:)';
PolyCoef1y=(PolyCalSystem1c\eye(size(PolyCalSystem1c,1)))*xi1s2c(2,:)';
FindEmpReprojError(xi1s2c',PolyCalSystem1c,PolyCoef1x,PolyCoef1y);

% Camera 2
[PolyCalSystem2c]=PolyCal2d7(s2c(1,:),s2c(2,:));
PolyCoef2x=(PolyCalSystem2c\eye(size(PolyCalSystem2c,1)))*xi2s2c(1,:)';
PolyCoef2y=(PolyCalSystem2c\eye(size(PolyCalSystem2c,1)))*xi2s2c(2,:)';
FindEmpReprojError(xi2s2c',PolyCalSystem2c,PolyCoef2x,PolyCoef2y);
    
%% Experiment

e_ds2(k)=normrnd(0,0.001);
e_um(k)=normrnd(0,0.001);
e_dRF(k)=e_ds2(k)+e_um(k);  % Micrometer Error + Displacement Sensor Error = difference of z=0 plane position in Calibration and Experimental RFs.
system_settings4_ExpRF;   % Regenerate the environment (Experimental RF)

% PROG=k/NoT*100;
% PROG=(k-start_k)/(NoT-start_k)*100;
% fprintf('Progress: %2.1f %%',PROG);

zmin=-1;
zmax=-8;
s=GenerateRandPattern(xmin,xmax,ymin,ymax,zmin,zmax,3000);

RInt=1;
[xi1r, xi2r]=Image_Formation_ver4(c1ExpRF, c2ExpRF, A1ExpRF,A2ExpRF, R1, R2, es13,es23, s,oi1ExpRF,oi2ExpRF,RInt,n);

% Find the z=0 image of the particles
[invPolyCalSystem1r]=invPolyCal2d7(xi1r(1,:),xi1r(2,:));
[invPolyCalSystem2r]=invPolyCal2d7(xi2r(1,:),xi2r(2,:));

sz01(1,:)=(invPolyCalSystem1r*invPolyCoef1x)';
sz01(2,:)=(invPolyCalSystem1r*invPolyCoef1y)';

sz02(1,:)=(invPolyCalSystem2r*invPolyCoef2x)';
sz02(2,:)=(invPolyCalSystem2r*invPolyCoef2y)';

% Triangulation with no dRF correction
[OPtvec]=TriangulateFromDewarpedImageVec(sz01,sz02,c1p,c2p,e3,n);

% Triangulation with dRF correction
% Correct the RF mismatch error (dRF) by minimizing the Triangulation Error

Find_dRF;   % Script to find dRF which minimizes triangulation error


siz=size(sz01,2);
sz01=[sz01;repmat(dRF(k),1,siz)];
sz02=[sz02;repmat(dRF(k),1,siz)];

c1p=c1p+dRF(k);
c2p=c2p+dRF(k);

u_i1=(sz01-repmat(c1p,1,siz))./repmat(normMat(sz01-repmat(c1p,1,siz)),3,1);
u_i2=(sz02-repmat(c2p,1,siz))./repmat(normMat(sz02-repmat(c2p,1,siz)),3,1);

l1=repmat(dRF(k),1,siz)./dot(-u_i1,repmat(e3,1,siz));
l2=repmat(dRF(k),1,siz)./dot(-u_i2,repmat(e3,1,siz));

sz01p=sz01+repmat(l1,3,1).*u_i1;
sz02p=sz02+repmat(l2,3,1).*u_i2;

[OPtvecCORR,l3]=TriangulateFromDewarpedImageVec(sz01p,sz02p,c1p,c2p,e3,n);

% Calculation of error
Error{k}=(s-OPtvec);
ErrorCORR{k}=(s-OPtvecCORR);

end
