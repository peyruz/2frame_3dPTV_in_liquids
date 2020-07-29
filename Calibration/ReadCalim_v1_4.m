%% Analyze the calibration image.
% The script rearranges the output the list of detected particles in order.
% Part of software to operate 3dPTV setup.
%
% Input
% CalimPath - path to the calibration image
% ParFilePath - path to the particle file
% colx - number of column in par file containing x-coordinates
% coly - number of column in par file containing y-coordinates
% varargin{1} = principalDots - vector of principal points (the user can input
%               these instead of prompting providing them graphically)
% varargin{2} - Name-Value argument to request visualization of the result.
%               eg : (____,'Visualize','on')
%               Allowed values: 'on','off'.
% 
% Output
% CalDot - ordered array of calibration dots
% CalCamMask - logical mask exposing the dots detected by the processor
% varargout{1} = RowNum - number of rows in the calibration pattern
% varargout{2} = principalDots - vector of principal points
%
% Changelog
% version 1.1
% Ability to input and output set of principal points. This feat is useful
% for processing multiple images. The user is prompted to provide graphical
% input (choice pattern) only once for the whole image set instead of for each image.
%
% version 1.2
% Added sorting tracker ixCol. Given index of a dot (x,:) in CalCam, we can
% find the correspondent (identical) dot in the raw array of detected
% points (output from image processing, here: 'com') as com(ixCol(x),:).
% This tracker can be used to sort e.g. contours of the particles in
% accordance with CalCam.
%
% version 1.3, May 27, 2018
% - Improved robustness of the code. Now the code keeps track of not only
% vertical data gaps but also horizontal for those cases where relevant.
% - ixCol is now preallocated for better speed and memory management.
% 
% version 1.4, May 29, 2018
% - added robustness through parsing and validation of inputs
% - added option of visualizing the result (see inputs)
%%
function [CalCam, CalCamMask,ixCol,varargout]=ReadCalim_v1_4(CalimPath,ParFilePath,colx,coly,varargin)

% Input processing
p=inputParser;

addRequired(p,'CalimPath',@(x)(isstring(x)||ischar(x)));
addRequired(p,'ParFilePath',@(x)(isstring(x)||ischar(x)));
addRequired(p,'colx',@(x)(x>=0&&isnumeric(x)&&(rem(x,1)==0)));
addRequired(p,'coly',@(x)(x>=0&&isnumeric(x)&&(rem(x,1)==0)));
addOptional(p,'principalDots',[],@(x)(all(size(x)==[5 2])));
addParameter(p,'Visualize','off',@(x)any(validatestring(x,[{'on'},{'off'}])));
p.CaseSensitive=false;

parse(p,CalimPath,ParFilePath,colx,coly,varargin{:});

CalimPath=p.Results.CalimPath;
ParFilePath=p.Results.ParFilePath;
colx=p.Results.colx;
coly=p.Results.coly;

% 10k is an arbitrary upper bound for a number of calibration dots.
% If the limit is exceeded, the code will stop execution. This number
% should be increased for larger calibration targets.
maxDotsAllowed=10000;

ixCol=zeros(maxDotsAllowed,1);
% Needs to be increased for larger calibration targets.

% Load dot position file
com=dlmreadPG(ParFilePath,'\t',[2 colx inf coly]);
% com=com+[1 1];  % Correct for the shift (for Insight data)

% If principal dots are provided in input - proceeed to reading the image
% (useful in MultReadCalim). Otherwise request a GUI input from the user.
if ~isempty(p.Results.principalDots)
    xp=varargin{1}(:,1);
    yp=varargin{1}(:,2);
    distThr=10; % px, radius of search for the principal point
    
    for ii=1:5
        [dist,ix]=min(normMat2d((com-[xp(ii),yp(ii)])'));
        if dist<distThr
            xp(ii)=com(ix,1);
            yp(ii)=com(ix,2);
            
            if ii==1
                ixCol(1)=ix;
            end
            
            if ii==3
                ixCol(2)=ix;
            end
            
        else
            error('The provided principal dots are not accurate enough.');
        end
    end
    
else
    % Load image
    Calim=imread(CalimPath);
    
    % Preallocate
    xp=zeros(5,1);
    yp=zeros(5,1);
    
    % Ask user to manually mark some principal dots
    fh=choicePatternFig(imread('choicePattern.jpg'), Calim, com(:,1), com(:,2));
    hd=helpdlg('Please choose the dots in the calibration image in the order indicated in the diagram on the left.','Choose');
    waitfor(hd);
    
    for i=1:5
        
        [xp(i),yp(i)]=ginput(1);
        [~,ix]=min(normMat2d((com-[xp(i),yp(i)])'));
        xp(i)=com(ix,1);
        yp(i)=com(ix,2);
        
        % Record the sorting indices
        if i==1
            ixCol(1)=ix;
        end
        
        if i==3
            ixCol(2)=ix;
        end
        
        scatter(xp(i),yp(i),[],'yellow','filled');
        
    end
    if nargout==5
        varargout{2}=[xp,yp];
    end
end

% Initialize the intervals
dx=[xp(2)-xp(1),yp(2)-yp(1)];
dy=[xp(3)-xp(1),yp(3)-yp(1)];

% Initialize Mask
CalCamMask=false(5776,1);

% 5776 is the maximum number of dots possible for this target. Changes might be needed for other targets.
CalCamMask(1)=1;
CalCamMask(2)=1;

% Dot search
% Guess accuracy threshold (px)
gThresh=25;

% Vertical Gap size limit (number of dots). Headnote [1]
kkLim=4;

% We start with probe point '1'
CalCam(1,1)=xp(1);
CalCam(1,2)=yp(1);
CalCam(2,1)=xp(3);
CalCam(2,2)=yp(3);

% Initialize the while loop
ii=2;
lastDotNotReached=1;
kkv=0;   % Vertical gap size counter
Ny=100000; % ie some big number clearly larger than the number of rows in the target
top=0;
NewCol=0;

while lastDotNotReached
    
    if ii+1==maxDotsAllowed
        error('Execution terminated. Maximum allowable number of calibration dots reached.');
    end
    
    if mod(ii+1,Ny)==1 || top==1
        NewCol=1;   % New column flag
        
        % Check the size of the horizontal gap in the current row
        kkh=0;
        while kkh<kkLim
            if (ii-(kkh+1)*Ny+1)>0 && CalCamMask(ii-(kkh+1)*Ny+1)==0
                kkh=kkh+1;
            else
                gapTooLarge=false;
                break;
            end
            % If the neighbor dot is out of bounds of the calibration pattern
            if (ii-(kkh+1)*Ny+1)<0 
                gapTooLarge=false;
                break;
            end
            
            if kkh==kkLim
                gapTooLarge=true;
            end
        end
        
        if gapTooLarge
            fprintf('%0.0f dots analyzed. \n',ii+1);
            error(strcat('There appears to be a large horizontal gap present in the detected calibration pattern. ',...
                         'For better results use another image or a different thresholding strategy.'));
        end
        
        if floor((ii+1)/Ny)~=1
            dx=CalCam(ii-Ny+1,:)-CalCam(ii-2*Ny+1,:);
            % else ii+1 is our second principal point
        end
        
        CalCam(ii+1,:)=CalCam(ii-Ny+1,:)+dx;
        [diff,ix]=min(normMat2d((com-CalCam(ii+1,:))'));
        
        if diff<gThresh
            CalCam(ii+1,:)=com(ix,:);
            CalCamMask(ii+1)=1;
            top=0;
            ixCol(ii+1)=ix;
        else
            CalCamMask(ii+1)=0;     % We retain our guess, however mark it as such in the mask
            top=1;
            ixCol(ii+1)=NaN;
        end
        
        ii=ii+1;
        continue
    end
    
    % intermediate dots
    if NewCol==1
        % Check the size of the horizontal gap in the current row
        kkh=0;
        while kkh<kkLim
            if (ii-(kkh+1)*Ny+1)>0 && CalCamMask(ii-(kkh+1)*Ny+1)==0
                kkh=kkh+1;
            else
                gapTooLarge=false;
                break;
            end
            if (ii-(kkh+1)*Ny+1)<0
                gapTooLarge=false;
                break;
            end
            
            if kkh==kkLim
                gapTooLarge=true;
            end
        end
        
        if gapTooLarge
            fprintf('%0.0f dots analyzed. \n',ii+1);
            error(strcat('There appears to be a large horizontal gap present in the detected calibration pattern. ',...
                         'For better results use another image or a different thresholding strategy.'));
        end
        
        
        dy=CalCam(ii-Ny+1,:)-CalCam(ii-Ny,:);
        NewCol=0;
    end
    
    CalCam(ii+1,:)=CalCam(ii,:)+dy;
    [diff,ix]=min(normMat2d((com-CalCam(ii+1,:))'));
    if diff<gThresh
        CalCam(ii+1,:)=com(ix,:);
        CalCamMask(ii+1)=1;
        dy=CalCam(ii+1,:)-CalCam(ii,:);
        kkv=0;
        ixCol(ii+1)=ix;
    else
        CalCamMask(ii+1)=0; % We retain our guess, however mark it as such in the mask
        kkv=kkv+1;
        ixCol(ii+1)=0;
        if kkv==kkLim
            fprintf('%0.0f dots analyzed. \n',ii+1);
            error(strcat('There appears to be a large vertical gap present in the detected calibration pattern. ',...
                'For better results use another image or a different thresholding strategy.'));
        end
    end
    
    if  CalCam(ii+1,:)==[xp(4) yp(4)]
        Ny=ii+1;  % column size found
    end
    
    % End criterion
    if CalCam(ii+1,:)==[xp(5) yp(5)]
        lastDotNotReached=0;
    end
    
    ii=ii+1;
end

if nargout>3
    varargout{1}=Ny;
end


if strcmp(p.Results.Visualize,'on')
    %% Visual check
    % The dots are plotted in order, check if the ordering is correct.
    CalDotValid=CalCam(CalCamMask,:);
    
    % We create a group not to generate a large number of handles and entries in the plot browser
    hg=hggroup;
    
    axis([0 6600 0 4400])
    hold on
    for ll =1:length(CalDotValid)
        h(ll)=scatter(CalDotValid(ll,1),CalDotValid(ll,2),'MarkerFaceColor','green','MarkerFaceAlpha',1);
        set(h(ll),'Parent',hg);
        drawnow limitrate;
    end
end

if exist('fh','var')
    close(fh)
end
end