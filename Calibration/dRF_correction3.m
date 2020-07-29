% Calculate the optimal position of the interface minimizing the
% discrepancy between the calibration and the experimental reference frames
% (dRF).
%
%
% Input:
% invPolyCoef+'i'+'+'j'+'k'  -  <vertical vec> empirical calibration
%                               coefficients. i - x,y; j - camera ID; k -
%                               calibration position number (always =1).
% 'Order'                    -  (Name-Value pair; Default='linear') Order of correction
% 
% Output:
% z-corr                     -  <scalar> Calculated optimal position of the
%                               interface.
% 
% Headnotes:
%           [1] :   Quadratic interpolation will produce two estimates. Typically those will
%                 	be quite different and it will be clear which one to use. The code will
%                   display and assign to the output the estimate which is closer to the
%                   sensor readings, input by the user. The second root will also be
%                   displayed separately just in case.
% 
% Headnotes:
%   [1] :   This method relies on the assumption that the calibration
%           target (and hence RF) is aligned with the cameras, ie the
%           cameras lie in or very closely in the yz of the calibration RF.
%           The more misaligned in that way the calibration target is, the
%           larger the allowable yLimit has to be set.
% 
% Changelog:
% version 2, Jan, 2018
% the script has been simplified to handle only a single cluster of
% particles and not a few as in the first version.
% 
% version 3, Jun 13, 2018
% Added option of 2nd order extrapolation. Use Name-Value pair
% 'Order','quadratic'
% 
% by Peyruz Gasimov, Jan 2018

function [z_corr]=dRF_correction3(invPolyCoefx11,invPolyCoefy11,invPolyCoefx21,invPolyCoefy21,varargin)

p=inputParser;

addRequired(p,'invPolyCoefx11',@(x)(isnumeric(x)&&isvector(x)));
addRequired(p,'invPolyCoefy11',@(x)(isnumeric(x)&&isvector(x)));
addRequired(p,'invPolyCoefx21',@(x)(isnumeric(x)&&isvector(x)));
addRequired(p,'invPolyCoefy21',@(x)(isnumeric(x)&&isvector(x)));
addParameter(p,'Order','linear',@(x)(any(validatestring(x,{'linear','quadratic'}))));

parse(p,invPolyCoefx11,invPolyCoefy11,invPolyCoefx21,invPolyCoefy21,varargin{:});

invPolyCoefx11=p.Results.invPolyCoefx11;
invPolyCoefy11=p.Results.invPolyCoefy11;
invPolyCoefx21=p.Results.invPolyCoefx21;
invPolyCoefy21=p.Results.invPolyCoefy21;


%% Image upload

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

dRFimPath11=cell(1,imNum11);

for ii=1:imNum11
dRFimPath11{ii}=strcat(PathName,fileName{ii});
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

dRFimPath12=cell(1,imNum12);

for ii=1:imNum12
dRFimPath12{ii}=strcat(PathName,fileName{ii});
end

if strcmp(p.Results.Order,'quadratic')
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
    
    dRFimPath13=cell(1,imNum13);
    
    for ii=1:imNum13
        dRFimPath13{ii}=strcat(PathName,fileName{ii});
    end
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

dRFimPath21=cell(1,imNum21);

for ii=1:imNum21
dRFimPath21{ii}=strcat(PathName,fileName{ii});
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

dRFimPath22=cell(1,imNum22);

for ii=1:imNum22
dRFimPath22{ii}=strcat(PathName,fileName{ii});
end

if strcmp(p.Results.Order,'quadratic')
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
    
    dRFimPath23=cell(1,imNum23);
    
    for ii=1:imNum23
        dRFimPath23{ii}=strcat(PathName,fileName{ii});
    end
end

%% Image Processing
batchName='dRFcorr'; % Batch name won't be used

optionsIP=struct( 'batchName',batchName, ...
                'startNum',1,...
                'detectCenters',true,...
                'ThreshMethod','scalar', ...
                'RSG_xPixelOverflowed',false, ...
                'minBlobDia',[], ...
                'maxBlobDia', [], ...
                'SubtractBackg',false, ...
                'Mask',true, ...
                'CenterFindingAlg','RadialSym/Gauss', ... % Possible options: 'Centroid','WeightedCentroid','RadialSym/Gauss'  
                'Smoothing',2,...
                'Average',true,...
                'StrictnessFactor',0.8,...
                'detectContours',false, ...
                'write2file',false, ...
                'saveDir',[]);
            
% Particle detection
[~,~,ImagePos11]=detectParticles5_4(dRFimPath11, optionsIP);
[~,~,ImagePos21]=detectParticles5_4(dRFimPath21, optionsIP);
[~,~,ImagePos12]=detectParticles5_4(dRFimPath12, optionsIP);
[~,~,ImagePos22]=detectParticles5_4(dRFimPath22, optionsIP);
if strcmp(p.Results.Order,'quadratic')
    [~,~,ImagePos13]=detectParticles5_4(dRFimPath13, optionsIP);
    [~,~,ImagePos23]=detectParticles5_4(dRFimPath23, optionsIP);
end

%% Dewarp the images
[invPolyCalSystem11]=invPolyCal2d7ver2(ImagePos11(:,1),ImagePos11(:,2));
sz11(:,1)=invPolyCalSystem11*invPolyCoefx11;
sz11(:,2)=invPolyCalSystem11*invPolyCoefy11;

[invPolyCalSystem21]=invPolyCal2d7ver2(ImagePos21(:,1),ImagePos21(:,2));
sz21(:,1)=invPolyCalSystem21*invPolyCoefx21;
sz21(:,2)=invPolyCalSystem21*invPolyCoefy21;

[invPolyCalSystem12]=invPolyCal2d7ver2(ImagePos12(:,1),ImagePos12(:,2));
sz12(:,1)=invPolyCalSystem12*invPolyCoefx11;
sz12(:,2)=invPolyCalSystem12*invPolyCoefy11;

[invPolyCalSystem22]=invPolyCal2d7ver2(ImagePos22(:,1),ImagePos22(:,2));
sz22(:,1)=invPolyCalSystem22*invPolyCoefx21;
sz22(:,2)=invPolyCalSystem22*invPolyCoefy21;

if strcmp(p.Results.Order,'quadratic')
    [invPolyCalSystem13]=invPolyCal2d7ver2(ImagePos13(:,1),ImagePos13(:,2));
    sz13(:,1)=invPolyCalSystem13*invPolyCoefx11;
    sz13(:,2)=invPolyCalSystem13*invPolyCoefy11;
    
    [invPolyCalSystem23]=invPolyCal2d7ver2(ImagePos23(:,1),ImagePos23(:,2));
    sz23(:,1)=invPolyCalSystem23*invPolyCoefx21;
    sz23(:,2)=invPolyCalSystem23*invPolyCoefy21;
end

%% Particle matching
% Position 1
IWsize=10;
Rnx=3;
Rny=3;

[NeighborMatrix1]=BuildNeighborMatrix4(sz11, IWsize,Rnx,Rny,1);

Rs=1;
[MatchMatrix1, NoMatchP1, Rs]=BuildInitialMatchMatrix3([sz11,zeros(size(sz11,1),1)],...
                                                       [sz21,zeros(size(sz21,1),1)], Rs);

Rq=0.005;    % We can be strict here, since we expect an almost parallel "motion" here.
[MatchMatrix1, NoMatchP1]=ParticleMatching(NeighborMatrix1,MatchMatrix1,NoMatchP1,...
                                            [sz11,zeros(size(sz11,1),1)],...
                                            [sz21,zeros(size(sz21,1),1)],Rq);

[velocityVecs1]=CalcVectors3([sz11,zeros(size(sz11,1),1)],...
                                                   	[sz21,zeros(size(sz21,1),1)],...
                                                   	MatchMatrix1,NeighborMatrix1,NoMatchP1);

% Position 2
[NeighborMatrix2]=BuildNeighborMatrix4(sz12, IWsize,Rnx,Rny,1);

[MatchMatrix2, NoMatchP2]=BuildInitialMatchMatrix3([sz12,zeros(size(sz12,1),1)],...
                                                       [sz22,zeros(size(sz22,1),1)], Rs);
                                                    
[MatchMatrix2, NoMatchP2]=ParticleMatching(NeighborMatrix2,MatchMatrix2,NoMatchP2,...
                                            [sz12,zeros(size(sz12,1),1)],...
                                            [sz22,zeros(size(sz22,1),1)],Rq);
                                        
[velocityVecs2]=CalcVectors3([sz12,zeros(size(sz12,1),1)],...
                                                  	[sz22,zeros(size(sz22,1),1)],...
                                                    MatchMatrix2,NeighborMatrix2,NoMatchP2);
                                                
if strcmp(p.Results.Order,'quadratic')                                                
    % Position 3
    [NeighborMatrix3]=BuildNeighborMatrix4(sz13, IWsize,Rnx,Rny,1);

    [MatchMatrix3, NoMatchP3]=BuildInitialMatchMatrix3([sz13,zeros(size(sz13,1),1)],...
                                                       [sz23,zeros(size(sz23,1),1)], Rs);
                                                    
    [MatchMatrix3, NoMatchP3]=ParticleMatching(NeighborMatrix3,MatchMatrix3,NoMatchP3,...
                                            [sz13,zeros(size(sz13,1),1)],...
                                            [sz23,zeros(size(sz23,1),1)],Rq);
                                        
    [velocityVecs3]=CalcVectors3([sz13,zeros(size(sz13,1),1)],...
                                                  	[sz23,zeros(size(sz23,1),1)],...
                                                    MatchMatrix3,NeighborMatrix3,NoMatchP3);                                                
end

% Filter out inprobable data based on limiting the y-deviation (see
% headnote 1)
yLimit1=0.05;
yLimit2=0.05;
if strcmp(p.Results.Order,'quadratic')
    yLimit3=0.05;
end

mask1=all([abs(velocityVecs1(:,2))<yLimit1, velocityVecs1(:,4)~=-1],2) ;
mask2=all([abs(velocityVecs2(:,2))<yLimit2, velocityVecs2(:,4)~=-1],2) ;
if strcmp(p.Results.Order,'quadratic')
mask3=all([abs(velocityVecs3(:,2))<yLimit3, velocityVecs3(:,4)~=-1],2) ;
end

dx1=median(velocityVecs1(mask1,1));
dx2=median(velocityVecs2(mask2,1));
if strcmp(p.Results.Order,'quadratic')
    dx3=median(velocityVecs3(mask3,1));
end

%% Calculate the optimal position
% Sensor readings for the above and below positions (mm)
z1_s=input('Input the sensor reading at position 1 (mm): ');
z2_s=input('Input the sensor reading at position 2 (mm): ');
if strcmp(p.Results.Order,'quadratic')
    z3_s=input('Input the sensor reading at position 3 (mm): ');
end

if strcmp(p.Results.Order,'linear')
    z_corr=z2_s+(z1_s-z2_s)/(dx1-dx2)*(-dx2);
elseif strcmp(p.Results.Order,'quadratic')
    polyCoef=[z1_s^2 z1_s 1; z2_s^2 z2_s 1; z3_s^2 z3_s 1 ]\[dx1; dx2; dx3];
    quadRoots=roots(polyCoef);
end

% For the quadratic order case, find the more plausable estimate 
% (Headnote [1])
if strcmp(p.Results.Order,'quadratic')
   [~,ix]=min(abs(quadRoots-z1_s));
   z_corr=quadRoots(ix);
end

% Display the result
fprintf('The optimal position is: %f mm\n', z_corr);
if strcmp(p.Results.Order,'quadratic')
   fprintf('Second Root: %f mm\n', quadRoots([1 2]~=ix)); 
end

figure
scatter([z1_s z2_s z3_s],[dx1 dx2 dx3]);
hold on
af=polyfit([z1_s z2_s z3_s],[dx1 dx2 dx3],1);
[yf]=polyval(af,[8,9]);
line([8,9],yf)
