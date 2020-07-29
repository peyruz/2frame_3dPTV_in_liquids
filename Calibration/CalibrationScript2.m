% CalibrationScript for 2 cameras, 3 calibration image positions for each
% by Peyruz Gasimov, Nov 2017
% 
% Headnotes:
%   {1} :   it is important to limit the principal distance from it reaching zero,
%           since zero constitutes a sort of a trivial solution (since it's zero
%           principal distance we get all points projected to a single point because it
%           essentially means that the center of projection itself is on the virtual
%           image plane) resulting in very low value of the objective function. We
%           therefore set the upper bound  at -15 for example.

% -------- Changelog
% version 2, May 24, 2018
% - code adjustments and bug fixes
% - add a choice of lens type model and whether perspective correction
% should be performed

% Calibration Options
optionsCal.LensType='Telecentric';
optionsCal.CorrectPerspective=true;

%% Load Images
% Camera 1
% Position 1
[fileName,PathName,fIndex] = uigetfile({'*.tif*';'*.tiff'},'Load the Image batch. Camera 1; Position 1.','Multiselect','on');
if fIndex==0
    error('No images selected. Restart the script.');
end

if ischar(fileName)
    imNum11=1;
else
    imNum11=size(fileName,2);
end

% In case only one image is selected, generate the expected cell array
if imNum11==1
    foo=fileName;
    fileName=cell(1,1);
    fileName{1}=foo;
end

CalimPath11=cell(1,imNum11);

for ii=1:imNum11
CalimPath11{ii}=strcat(PathName,fileName{ii});
end

% Position 2
[fileName,PathName,fIndex] = uigetfile({'*.tif*';'*.tiff'},'Load the Image batch. Camera 1; Position 2.','Multiselect','on',PathName);
if fIndex==0
    error('No images selected. Restart the script.');
end

if ischar(fileName)
    imNum12=1;
else
    imNum12=size(fileName,2);
end

% In case only one image is selected, generate the expected cell array
if imNum12==1
    foo=fileName;
    fileName=cell(1,1);
    fileName{1}=foo;
end

CalimPath12=cell(1,imNum12);

for ii=1:imNum12
CalimPath12{ii}=strcat(PathName,fileName{ii});
end

% Position 3
[fileName,PathName,fIndex] = uigetfile({'*.tif*';'*.tiff'},'Load the Image batch. Camera 1; Position 3.','Multiselect','on',PathName);
if fIndex==0
    error('No images selected. Restart the script.');
end

if ischar(fileName)
    imNum13=1;
else
    imNum13=size(fileName,2);
end

% In case only one image is selected, generate the expected cell array
if imNum13==1
    foo=fileName;
    fileName=cell(1,1);
    fileName{1}=foo;
end

CalimPath13=cell(1,imNum13);

for ii=1:imNum13
CalimPath13{ii}=strcat(PathName,fileName{ii});
end

% Camera 2
% Position 1
[fileName,PathName,fIndex] = uigetfile({'*.tif*';'*.tiff'},'Load the Image batch. Camera 2; Position 1.','Multiselect','on',PathName);
if fIndex==0
    error('No images selected. Restart the script.');
end

if ischar(fileName)
    imNum21=1;
else
    imNum21=size(fileName,2);
end

% In case only one image is selected, generate the expected cell array
if imNum21==1
    foo=fileName;
    fileName=cell(1,1);
    fileName{1}=foo;
end

CalimPath21=cell(1,imNum21);

for ii=1:imNum21
CalimPath21{ii}=strcat(PathName,fileName{ii});
end

% Position 2
[fileName,PathName,fIndex] = uigetfile({'*.tif*';'*.tiff'},'Load the Image batch. Camera 2; Position 2.','Multiselect','on',PathName);
if fIndex==0
    error('No images selected. Restart the script.');
end

if ischar(fileName)
    imNum22=1;
else
    imNum22=size(fileName,2);
end

% In case only one image is selected, generate the expected cell array
if imNum22==1
    foo=fileName;
    fileName=cell(1,1);
    fileName{1}=foo;
end

CalimPath22=cell(1,imNum22);

for ii=1:imNum22
CalimPath22{ii}=strcat(PathName,fileName{ii});
end

% Position 3
[fileName,PathName,fIndex] = uigetfile({'*.tif*';'*.tiff'},'Load the Image batch. Camera 2; Position 3.','Multiselect','on',PathName);
if fIndex==0
    error('No images selected. Restart the script.');
end

if ischar(fileName)
    imNum23=1;
else
    imNum23=size(fileName,2);
end

% In case only one image is selected, generate the expected cell array
if imNum23==1
    foo=fileName;
    fileName=cell(1,1);
    fileName{1}=foo;
end

CalimPath23=cell(1,imNum23);

for ii=1:imNum23
CalimPath23{ii}=strcat(PathName,fileName{ii});
end

%% Calibration Dot Detection
% Configure options (see fields' explanation in the detectParticles.m help)
optionsIP=struct( 'batchName','placeHolder', ...
                'startNum',1,...
                'detectCenters',true,...
                'ThreshMethod','adaptive', ...
                'RSG_xPixelOverflowed',false, ...
                'minBlobDia',[], ...
                'maxBlobDia', [], ...
                'SubtractBackg',false, ...
                'Mask',true, ...
                'CenterFindingAlg','Centroid', ...
                'Filter','Gaussian',...
                'NumOfFilterIterations',2,...
                'Average',false,...
                'StrictnessFactor',[],...
                'detectContours',true, ...
                'write2file',true, ...
                'saveDir',[]);
            
% Detect edges only of perspective correction is requested
if optionsCal.CorrectPerspective
    optionsIP.detectContours=true;
else
    optionsIP.detectContours=false;
end

% Particle detection
% Camera 1
% Position 1
optionsIP.batchName='Calib_Cam01_Pos01';
[~, contours11]=detectParticles5_5(CalimPath11, optionsIP);

% Position 2
optionsIP.batchName='Calib_Cam01_Pos02';    
[~, contours12]=detectParticles5_5(CalimPath12, optionsIP);

% Position 3
optionsIP.batchName='Calib_Cam01_Pos03';
[~, contours13]=detectParticles5_5(CalimPath13, optionsIP);

% Camera 2
% Position 1
optionsIP.batchName='Calib_Cam02_Pos01';
[~, contours21]=detectParticles5_5(CalimPath21, optionsIP);

% Position 2
optionsIP.batchName='Calib_Cam02_Pos02';
[~, contours22]=detectParticles5_5(CalimPath22, optionsIP);

% Position 3
optionsIP.batchName='Calib_Cam02_Pos03';
[~, contours23]=detectParticles5_5(CalimPath23, optionsIP);

%% Load the Dot files
% Camera 1
% Position 1
[fileName,PathName,fIndex] = uigetfile({'*.3dp'},'Load Dot Files. Camera 1; Position 1.','Multiselect','on');
if fIndex==0
    error('No files selected. Try again.');
end

ParFilePath11=cell(1,imNum11);

for ii=1:imNum11
    ParFilePath11{ii}=strcat(PathName,fileName{ii});
end

% Position 2
[fileName,PathName,fIndex] = uigetfile({'*.3dp'},'Load Dot Files. Camera 1; Position 2.','Multiselect','on',PathName);
if fIndex==0
    error('No files selected. Try again.');
end

ParFilePath12=cell(1,imNum12);

for ii=1:imNum12
    ParFilePath12{ii}=strcat(PathName,fileName{ii});
end

% Position 3
[fileName,PathName,fIndex] = uigetfile({'*.3dp'},'Load Dot Files. Camera 1; Position 3.','Multiselect','on',PathName);
if fIndex==0
    error('No files selected. Try again.');
end

ParFilePath13=cell(1,imNum13);

for ii=1:imNum13
    ParFilePath13{ii}=strcat(PathName,fileName{ii});
end

% Camera 2
% Position 1
[fileName,PathName,fIndex] = uigetfile({'*.3dp'},'Load Dot Files. Camera 2; Position 1.','Multiselect','on',PathName);
if fIndex==0
    error('No files selected. Try again.');
end

ParFilePath21=cell(1,imNum21);

for ii=1:imNum21
    ParFilePath21{ii}=strcat(PathName,fileName{ii});
end

% Position 2
[fileName,PathName,fIndex] = uigetfile({'*.3dp'},'Load Dot Files. Camera 2; Position 2.','Multiselect','on',PathName);
if fIndex==0
    error('No files selected. Try again.');
end

ParFilePath22=cell(1,imNum22);

for ii=1:imNum22
    ParFilePath22{ii}=strcat(PathName,fileName{ii});
end

% Position 3
[fileName,PathName,fIndex] = uigetfile({'*.3dp'},'Load Dot Files. Camera 2; Position 3.','Multiselect','on',PathName);
if fIndex==0
    error('No files selected. Try again.');
end

ParFilePath23=cell(1,imNum23);

for ii=1:imNum23
    ParFilePath23{ii}=strcat(PathName,fileName{ii});
end

%% Process the dot data
colx=0;
coly=1;

% Camera 1
[CalCamMean11, CalCamMaskAll11,ixCol11, RowNum11]=MultReadCalim(CalimPath11{1},ParFilePath11,colx,coly);
[CalCamMean12, CalCamMaskAll12,ixCol12, RowNum12]=MultReadCalim(CalimPath12{1},ParFilePath12,colx,coly);
[CalCamMean13, CalCamMaskAll13,ixCol13, RowNum13]=MultReadCalim(CalimPath13{1},ParFilePath13,colx,coly);

% Camera 2
[CalCamMean21, CalCamMaskAll21,ixCol21, RowNum21]=MultReadCalim(CalimPath21{1},ParFilePath21,colx,coly);
[CalCamMean22, CalCamMaskAll22,ixCol22, RowNum22]=MultReadCalim(CalimPath22{1},ParFilePath22,colx,coly);
[CalCamMean23, CalCamMaskAll23,ixCol23, RowNum23]=MultReadCalim(CalimPath23{1},ParFilePath23,colx,coly);

%% Arrange the contours' data with accordance to the dot center data
if optionsCal.CorrectPerspective
    % Camera 1
    % Position 1
    [contSort11]=arrangeContours(contours11, ixCol11,CalCamMaskAll11);
    % Position 2
    [contSort12]=arrangeContours(contours12, ixCol12,CalCamMaskAll12);
    % Position 3
    [contSort13]=arrangeContours(contours13, ixCol13,CalCamMaskAll13);
    
    % Camera 2
    % Position 1
    [contSort21]=arrangeContours(contours21, ixCol21,CalCamMaskAll21);
    % Position 2
    [contSort22]=arrangeContours(contours22, ixCol22,CalCamMaskAll22);
    % Position 3
    [contSort23]=arrangeContours(contours23, ixCol23,CalCamMaskAll23);
end

%% Create correspondent array of physical dot coordinates
elev=[0 -2.8896 2.0362];
dg=2;

% Camera 1 
[CalPhys11]=GenPhysDotArray(dg,elev,CalCamMaskAll11,false,false,RowNum11,RowNum12,RowNum13);
[CalPhys12]=GenPhysDotArray(dg,elev,false,CalCamMaskAll12,false,RowNum11,RowNum12,RowNum13);
[CalPhys13]=GenPhysDotArray(dg,elev,false,false,CalCamMaskAll13,RowNum11,RowNum12,RowNum13);

% Camera 2
[CalPhys21]=GenPhysDotArray(dg,elev,CalCamMaskAll21,false,false,RowNum21,RowNum22,RowNum23);
[CalPhys22]=GenPhysDotArray(dg,elev,false,CalCamMaskAll22,false,RowNum21,RowNum22,RowNum23);
[CalPhys23]=GenPhysDotArray(dg,elev,false,false,CalCamMaskAll23,RowNum21,RowNum22,RowNum23);

CalPhys1=[CalPhys11; CalPhys12; CalPhys13];
CalPhys2=[CalPhys21; CalPhys22; CalPhys23];

%% Perspective and distortion bias correction; Unbiased Empirical Calibration
if optionsCal.CorrectPerspective
    iterNum=1;
    
    % Camera 1
    [CorrCalCamMean11,baseMask11, invPolyCoefx11, invPolyCoefy11, PolyCoefx11, PolyCoefy11]=UnbiasedEmpCal(CalCamMean11, CalPhys11, contSort11, iterNum);
    [CorrCalCamMean12,baseMask12, ~, ~, ~, ~]=UnbiasedEmpCal(CalCamMean12, CalPhys12, contSort12, iterNum);
    [CorrCalCamMean13,baseMask13, ~, ~, ~, ~]=UnbiasedEmpCal(CalCamMean13, CalPhys13, contSort13, iterNum);
    
    % Camera 2
    [CorrCalCamMean21,baseMask21, invPolyCoefx21, invPolyCoefy21, PolyCoefx21, PolyCoefy21]=UnbiasedEmpCal(CalCamMean21, CalPhys21, contSort21, iterNum);
    [CorrCalCamMean22,baseMask22, ~, ~, ~, ~]=UnbiasedEmpCal(CalCamMean22, CalPhys22, contSort22, iterNum);
    [CorrCalCamMean23,baseMask23, ~, ~, ~, ~]=UnbiasedEmpCal(CalCamMean23, CalPhys23, contSort23, iterNum);
    
    CorrCalCamMean1=[CorrCalCamMean11; CorrCalCamMean12; CorrCalCamMean13];
    CorrCalCamMean2=[CorrCalCamMean21; CorrCalCamMean22; CorrCalCamMean23];
    
    CalCamMean1=[CalCamMean11; CalCamMean12; CalCamMean13];
    CalCamMean2=[CalCamMean21; CalCamMean22; CalCamMean23];
    
    baseMask1=[baseMask11'; baseMask12'; baseMask13'];
    baseMask2=[baseMask21'; baseMask22'; baseMask23'];
    
else
    % Biased empirical calibration
    % Camera 1
    [PolyCalSystemCent1]=invPolyCal2d7ver2(CalCamMean11(:,1),CalCamMean11(:,2));
    invPolyCoefx11=(PolyCalSystemCent1\eye(size(PolyCalSystemCent1,1)))*CalPhys11(:,1);
    invPolyCoefy11=(PolyCalSystemCent1\eye(size(PolyCalSystemCent1,1)))*CalPhys11(:,2);

    PolyCalSystemW01=PolyCal2d7ver2(CalPhys11(:,1),CalPhys11(:,2));
    PolyCoefx11=(PolyCalSystemW01\eye(size(PolyCalSystemW01,1)))*CalCamMean11(:,1);
    PolyCoefy11=(PolyCalSystemW01\eye(size(PolyCalSystemW01,1)))*CalCamMean11(:,2);
    
    % Camera 2
    [PolyCalSystemCent2]=invPolyCal2d7ver2(CalCamMean21(:,1),CalCamMean21(:,2));
    invPolyCoefx21=(PolyCalSystemCent2\eye(size(PolyCalSystemCent2,1)))*CalPhys21(:,1);
    invPolyCoefy21=(PolyCalSystemCent2\eye(size(PolyCalSystemCent2,1)))*CalPhys21(:,2);

    PolyCalSystemW02=PolyCal2d7ver2(CalPhys21(:,1),CalPhys21(:,2));
    PolyCoefx21=(PolyCalSystemW02\eye(size(PolyCalSystemW02,1)))*CalCamMean21(:,1);
    PolyCoefy21=(PolyCalSystemW02\eye(size(PolyCalSystemW02,1)))*CalCamMean21(:,2);
    
end

%% Filter out the data outliers of the perspective correction
if optionsCal.CorrectPerspective
    % The outliers are found based on the difference between the raw and
    % corrected dot coordinates. Dots with too large correction are suspicious
    % (both the raw and the corrected ones) and will be excluded.
    
    % Compute Euclidian distance between the corected and raw coordinates and
    % use the cut-iff thresold to define the correction Mask:
    
    % Prompt input of correction cut-off value
    % Camera 1
    fh3=figure;
    histogram(normMat2d(CalCamMean1(baseMask1,:)'-CorrCalCamMean1(baseMask1,:)'));
    title('Difference between the raw and the corrected dot coordinates. Camera 1.');
    xlabel('Difference, px');
    ylabel('Count');
    corrThresh1=input('Please input the correction threshold for Camera 1 (px): ');
    if ishandle(fh3)
        close(fh3)
    end
    
    % Camera 2
    fh3=figure;
    histogram(normMat2d(CalCamMean2(baseMask2,:)'-CorrCalCamMean2(baseMask2,:)'));
    title('Difference between the raw and the corrected dot coordinates. Camera 2.');
    xlabel('Difference, px');
    ylabel('Count');
    corrThresh2=input('Please input the correction threshold for Camera 2 (px): ');
    if ishandle(fh3)
        close(fh3)
    end
    
    % Constrtuct the validation Mask based on difference magnitude
    magMask11=normMat2d(CalCamMean11(baseMask11,:)'-CorrCalCamMean11(baseMask11,:)')<corrThresh1;
    magMask12=normMat2d(CalCamMean12(baseMask12,:)'-CorrCalCamMean12(baseMask12,:)')<corrThresh1;
    magMask13=normMat2d(CalCamMean13(baseMask13,:)'-CorrCalCamMean13(baseMask13,:)')<corrThresh1;
    
    magMask21=normMat2d(CalCamMean21(baseMask21,:)'-CorrCalCamMean21(baseMask21,:)')<corrThresh2;
    magMask22=normMat2d(CalCamMean22(baseMask22,:)'-CorrCalCamMean22(baseMask22,:)')<corrThresh2;
    magMask23=normMat2d(CalCamMean23(baseMask23,:)'-CorrCalCamMean23(baseMask23,:)')<corrThresh2;
    
    % Incorporate magMask into baseMask
    baseMask11(baseMask11)=magMask11;
    baseMask12(baseMask12)=magMask12;
    baseMask13(baseMask13)=magMask13;
    
    baseMask21(baseMask21)=magMask21;
    baseMask22(baseMask22)=magMask22;
    baseMask23(baseMask23)=magMask23;
    
end

%% Validation of the perspective correction
if optionsCal.CorrectPerspective
    % Choose neighborhood size in pixels (see validateVelField help)
    Rn=150; % px
    maxDev=2; % [0->inf] decreasing strictness
    
    fprintf('--- Unbiased Calibration --- \n');
    % Camera 1
    fprintf('--- Camera 1    Position 1 \n');
    [validMask11]=validateVelField(CalCamMean11(baseMask11,:),CorrCalCamMean11(baseMask11,:)-CalCamMean11(baseMask11,:),Rn,maxDev,'Visualize','on');
    fprintf('--- Camera 1    Position 2 \n');
    [validMask12]=validateVelField(CalCamMean12(baseMask12,:),CorrCalCamMean12(baseMask12,:)-CalCamMean12(baseMask12,:),Rn,maxDev,'Visualize','on');
    fprintf('--- Camera 1    Position 3 \n');
    [validMask13]=validateVelField(CalCamMean13(baseMask13,:),CorrCalCamMean13(baseMask13,:)-CalCamMean13(baseMask13,:),Rn,maxDev,'Visualize','on');
    
    % Camera 2
    fprintf('--- Camera 2    Position 1 \n');
    [validMask21]=validateVelField(CalCamMean21(baseMask21,:),CorrCalCamMean21(baseMask21,:)-CalCamMean21(baseMask21,:),Rn,maxDev,'Visualize','on');
    fprintf('--- Camera 2    Position 2 \n');
    [validMask22]=validateVelField(CalCamMean22(baseMask22,:),CorrCalCamMean22(baseMask22,:)-CalCamMean22(baseMask22,:),Rn,maxDev,'Visualize','on');
    fprintf('--- Camera 2    Position 3 \n');
    [validMask23]=validateVelField(CalCamMean23(baseMask23,:),CorrCalCamMean23(baseMask23,:)-CalCamMean23(baseMask23,:),Rn,maxDev,'Visualize','on');
    
    % Incorporate the validMask into baseMask
    baseMask11(baseMask11)=validMask11;
    baseMask12(baseMask12)=validMask12;
    baseMask13(baseMask13)=validMask13;
    
    baseMask21(baseMask21)=validMask21;
    baseMask22(baseMask22)=validMask22;
    baseMask23(baseMask23)=validMask23;
    
    baseMask1=[baseMask11'; baseMask12'; baseMask13'];
    baseMask2=[baseMask21'; baseMask22'; baseMask23'];
    
    CorrCalCamMean1=[CorrCalCamMean11; CorrCalCamMean12; CorrCalCamMean13];
    CorrCalCamMean2=[CorrCalCamMean21; CorrCalCamMean22; CorrCalCamMean23];
    
else
    CorrCalCamMean1=CalCamMean1; % actually not corrected, but this simplifies the code
    CorrCalCamMean2=CalCamMean2;
    
    baseMask1=true(size(CorrCalCamMean1,1),1);
    baseMask2=true(size(CorrCalCamMean2,1),1);
end


%% Estimate principal point (equivalent to center of distortion in Steger camera models)

% Initial guess of the center of distortion
x0ini1=[3300; 2200];
x0ini2=[3300; 2200];

% Estimate using radial distortion
x0ini1=EstimDistCen(invPolyCoefx11, invPolyCoefy11,x0ini1(1),x0ini1(2));
x0ini2=EstimDistCen(invPolyCoefx21, invPolyCoefy21,x0ini2(1),x0ini2(2));

%% Physics-based Calibration
ms = MultiStart('UseParallel',true);

optionsLsq = optimoptions('lsqnonlin','Display','iter',...
                            'TypicalX', [1e-6 1e-8, 1e-8 1e-8 1e-10 0.01 1e-4 1e-4 0.5 0.5 0.5 0.5 10 10 100 4.5 0.1 50 0.1 0.005 0.005  3000 2000],...
                            'MaxFunctionEvaluations',50000,'MaxIterations', 200,'StepTolerance',1e-10,'FunctionTolerance',1e-8,'UseParallel',false);
% Camera 1
% Initial guesses
pxsxini=0.0055;
pxsyini=0.0055;

% Initial estimates
[fini1,skini1,tini1,qini1]=TsaiLin(CorrCalCamMean1(baseMask1,:)',CalPhys1(baseMask1,:)',pxsxini, pxsyini,x0ini1');

qvini1=[real(qini1) (vector(qini1))'];

% Non-linear optimization (see Headnote {1})
if strcmp(optionsCal.LensType,'Telecentric')
    f1 = @(x1)StegerObjFun_Telc(x1,CalPhys1(baseMask1,:)',CorrCalCamMean1(baseMask1,:)');
    xini1=  [   0       0       0       0       0 0 0 0 qvini1(1)     qvini1(2)       qvini1(3)       qvini1(4)       tini1(1)        tini1(2)        tini1(3)        0.5*pi      degtorad(3)     fini1   0 pxsxini   pxsyini     x0ini1(1)       x0ini1(2)];
    lb1=    [  -0.001  -0.001  -0.001  -0.001   0 0 0 0 xini1(9)-0.2  xini1(10)-0.2   xini1(11)-0.2   xini1(12)-0.2   tini1(1)-100    tini1(2)-100    tini1(3)-100    0           0              -100     0 0.0055    0.0055      x0ini1(1)-1000  x0ini1(2)-500];
    ub1=    [   0.001   0.001   0.001   0.001   0 0 0 0 xini1(9)+0.2  xini1(10)+0.2   xini1(11)+0.2   xini1(12)+0.2   tini1(1)+100    tini1(2)+100    tini1(3)+100    pi          deg2rad(6)     -35      0 0.0055    0.0055      x0ini1(1)+1000  x0ini1(2)+500];
elseif strcmp(optionsCal.LensType,'Entocentric')
    f1 = @(x1)StegerObjFun_Entc(x1,CalPhys1(baseMask1,:)',CorrCalCamMean1(baseMask1,:)');
    xini1=  [   0      0      0       0       0 0 0 0 qvini1(1)     qvini1(2)       qvini1(3)       qvini1(4)       tini1(1)        tini1(2)        tini1(3)        0.5*pi      degtorad(2)     fini1  -50      pxsxini   pxsyini     x0ini1(1)       x0ini1(2)];
    lb1=    [  -0.001 -0.001 -1      -1       0 0 0 0 xini1(9)-0.2  xini1(10)-0.2   xini1(11)-0.2   xini1(12)-0.2   tini1(1)-100    tini1(2)-100    tini1(3)-100    0.5*pi-0.2  deg2rad(0)     -100    -inf     0.0055    0.0055      x0ini1(1)-1000  x0ini1(2)-500];
    ub1=    [   0.001  0.001  0.001   0.001   0 0 0 0 xini1(9)+0.2  xini1(10)+0.2   xini1(11)+0.2   xini1(12)+0.2   tini1(1)+100    tini1(2)+100    tini1(3)+100    0.5*pi+0.2  deg2rad(4)     -35      inf     0.0055    0.0055      x0ini1(1)+1000  x0ini1(2)+500];
else
    error('Unknown Lens Type. Accepted types: ''Telecentric'' and ''Entocentric''. ');
end

problem1 = createOptimProblem('lsqnonlin','objective',...
                                f1,'x0',xini1,'lb',lb1,'ub',ub1,'options',optionsLsq);

NumOfStarts=20;    
tic
[x1,f1,exitflag1,output1,solutions1] = run(ms,problem1,NumOfStarts);
toc

q1=quaternion([x1(9) x1(10) x1(11) x1(12)]);
t1=[x1(13);x1(14);x1(15)];
c1=-RotateVector(inverse(q1),t1);

% Calculate reprojection error
if strcmp(optionsCal.LensType,'Telecentric')
    [ReprojError1]=FindReprojError4StegerTelc(CalPhys1(baseMask1,:)',...
                                                x1,CorrCalCamMean1(baseMask1,:)');
elseif strcmp(optionsCal.LensType,'Entocentric')
    [ReprojError1]=FindReprojError4StegerEntc(CalPhys11(baseMask1,:)',...
                                                x1,CorrCalCamMean11(baseMask1,:)');
end

% Camera 2
% Initial estimates
[fini2,skini2,tini2,qini2]=TsaiLin(CorrCalCamMean2(baseMask2,:)',CalPhys2(baseMask2,:)',pxsxini, pxsyini,x0ini2');

qvini2=[real(qini2) (vector(qini2))'];

% Non-linear optimization (see Headnote {1})
if strcmp(optionsCal.LensType,'Telecentric')
    f2 = @(x2)StegerObjFun_Telc(x2,CalPhys2(baseMask2,:)',CorrCalCamMean2(baseMask2,:)');
    xini2=  [ 0         0       0       0       0 0 0 0 qvini2(1)       qvini2(2)       qvini2(3)       qvini2(4)       tini2(1)        tini2(2)        tini2(3)        3/2*pi    degtorad(3)     fini2    0   pxsxini pxsyini x0ini2(1)       x0ini2(2)];
    lb2=    [-0.001    -0.001  -0.001  -0.001   0 0 0 0 xini2(9)-0.2    xini2(10)-0.2   xini2(11)-0.2   xini2(12)-0.2   tini2(1)-100    tini2(2)-100    tini2(3)-100    pi        0              -100      0   0.0055  0.0055  x0ini2(1)-1000  x0ini2(2)-500];
    ub2=    [0.001      0.001   0.001   0.001   0 0 0 0 xini2(9)+0.2    xini2(10)+0.2   xini2(11)+0.2   xini2(12)+0.2   tini2(1)+100    tini2(2)+100    tini2(3)+100    2*pi      deg2rad(6)     -35       0   0.0055  0.0055  x0ini2(1)+1000  x0ini2(2)+500];
elseif strcmp(optionsCal.LensType,'Entocentric')
    xini2=  [   0       0       0       0       0 0 0 0 qvini2(1)       qvini2(2)       qvini2(3)       qvini2(4)       tini2(1)        tini2(2)        tini2(3)        3/2*pi  degtorad(3) fini2   -50  pxsxini pxsyini x0ini2(1)       x0ini2(2)];
    lb2=    [  -0.001  -0.001  -1      -1       0 0 0 0 xini2(9)-0.2    xini2(10)-0.2   xini2(11)-0.2   xini2(12)-0.2   tini2(1)-200    tini2(2)-200    tini2(3)-200    pi      deg2rad(0) -100     -inf 0.0055  0.0055  x0ini2(1)-1000  x0ini2(2)-500];
    ub2=    [   0.001   0.001   0.001   0.001   0 0 0 0 xini2(9)+0.2    xini2(10)+0.2   xini2(11)+0.2   xini2(12)+0.2   tini2(1)+200    tini2(2)+200    tini2(3)+200    2*pi    deg2rad(6) -35       inf 0.0055  0.0055  x0ini2(1)+1000  x0ini2(2)+500];
end

problem2 = createOptimProblem('lsqnonlin','objective',...
                                f2,'x0',xini2,'lb',lb2,'ub',ub2,'options',optionsLsq);
  
tic
[x2,f2,exitflag2,output2,solutions2] = run(ms,problem2,NumOfStarts);
toc

q2=quaternion([x2(9) x2(10) x2(11) x2(12)]);
t2=[x2(13);x2(14);x2(15)];
c2=-RotateVector(inverse(q2),t2);

% Calculate reprojection error
if strcmp(optionsCal.LensType,'Telecentric')
[ReprojError2]=FindReprojError4StegerTelc(CalPhys2(baseMask2,:)',...
                                          x2,CorrCalCamMean2(baseMask2,:)');
elseif strcmp(optionsCal.LensType,'Entocentric')
[ReprojError2]=FindReprojError4StegerEntc(CalPhys2(baseMask2,:)',...
                                          x2,CorrCalCamMean2(baseMask2,:)');                               
end

%% Diagnostic bar plot of the solution
figure
bar((x1-(ub1+lb1)/2)./((ub1-lb1)/2))
hold on
bar((x2-(ub2+lb2)/2)./((ub2-lb2)/2))

%% Generate data structures with the camera calibration parameters
CamParam1=struct('Camera_ID',1,'Position',c1,'Direct_Empirical_Mapping_x',PolyCoefx11,'Direct_Empirical_Mapping_y',PolyCoefy11,...
                'Inverse_Empirical_Mapping_x',invPolyCoefx11,'Inverse_Empirical_Mapping_y',invPolyCoefy11,...
                'Mean_Reprojection_Error_px',mean(ReprojError1));
             
CamParam2=struct('Camera_ID',2,'Position',c2,'Direct_Empirical_Mapping_x',PolyCoefx21,'Direct_Empirical_Mapping_y',PolyCoefy21,...
                'Inverse_Empirical_Mapping_x',invPolyCoefx21,'Inverse_Empirical_Mapping_y',invPolyCoefy21,...
                'Mean_Reprojection_Error_px',mean(ReprojError2));
