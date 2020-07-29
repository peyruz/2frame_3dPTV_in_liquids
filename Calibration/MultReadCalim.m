% Read multiple images of the same position and find averaged positions of
% the dots. Additionally, as an option, output mean standard errors for
% each particle accross the set of images
% 
% by Peyruz Gasimov, Aug 2017
% 
% Note: the difference between ReadCalim and MultReadCalim is that the
% output CalCamMean contains only the valid, detected particles, whereas
% all calibration dots (even guesses for undetected ones) up until the last
% detected dot are output by ReadCalim
% 
% Usage Examples
% [CalCamMean, CalCamMaskAll, ixCol, RowNum]=MultReadCalim(imNum,CalimPath,DotPathName,colx,coly);
% [CalCamMean, CalCamMaskAll, ixCol, RowNum,xmat,ymat]=MultReadCalim(imNum,CalimPath,DotPathName,colx,coly);
% [CalCamMean, CalCamMaskAll, ixCol,RowNum,xmat,ymat,xerr,yerr]=MultReadCalim(imNum,CalimPath,DotPathName,colx,coly);
% 
% Input
% imNum -           number of images
% CalimPath -       directory to the representative image for the set
% DotPathName -     cell array of directories to the particle files
% colx, coly -      columns of the par file containing the x and y coordinates
%                   respectively
% 
% Output
% CalCamMean -      mean dot coordinates
% CalCamMaskAll -   Mask exposing the good dots of the pattern
% ixCol -           cell array with size (1,imNum). ixOCol contains
%                   information on which contours correspond to which dots.
%                   This array is used by 'arrangeContours' function which
%                   produces such a cell array of contours (contSort) where
%                   every contSort{:}{kk} contains contours of dot
%                   CalCamMean(kk,:).
% RowNum -          Number of rows in the calibration pattern
% xmat, ymat -      MxN matrix of dot coordinates from all images. M -
%                   number of dots, N - number of images
% xerr, yerr -      standard centroid finding error for each dot accross 
%                   the image set

function [CalCamMean, CalCamMaskAll, ixCol, varargout]=MultReadCalim(CalimPath,ParFilePath,colx,coly)

% Determine the number of images
imNum=max(size(ParFilePath));

% Preallocate
ixCol=cell(imNum,1);
x=cell(imNum,1);

% Upload and read the dot files
[CalCam, CalCamMask(:,1),ixCol{1},RowNum,principalDots]=ReadCalim_v1_4(CalimPath,ParFilePath{1},colx,coly,'Visualize','off');
x{1}=CalCam;

% The selected points are are used for the remaining dot files. This speeds
% up the process
for ii=2:imNum
    [CalCam, CalCamMask(:,ii),ixCol{ii}]=ReadCalim_v1_4(CalimPath,ParFilePath{ii},colx,coly,principalDots);
    x{ii}=CalCam;
end

% Mask in the dots appearing in all images
CalCamMaskAll=all(CalCamMask,2);

% Preaalocate xmat, ymat
xmat=zeros(sum(CalCamMaskAll),imNum);
ymat=zeros(sum(CalCamMaskAll),imNum);

% Collect x and y into a matrix
for ii=1:imNum
    xmat(:,ii)=x{ii}(CalCamMaskAll,1);  % the undetected dots are removed
    ymat(:,ii)=x{ii}(CalCamMaskAll,2);
end

% Mean coordinates
CalCamMean(:,1)=mean(xmat,2);
CalCamMean(:,2)=mean(ymat,2);

if nargout>=5
    varargout{2}=xmat;
    varargout{3}=ymat;
end

if nargout>=4
    varargout{1}=RowNum;
end

% Report errors and output if prompted
if nargout==7
    varargout{4}=std(xmat,[],2);
    varargout{5}=std(ymat,[],2);
    fprintf('Mean standard error in x: %f px \n',mean(varargout{2}))
    fprintf('Mean standard error in y: %f px \n',mean(varargout{3}))
end




