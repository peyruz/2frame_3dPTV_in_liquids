% Detect particles in the given image(s)
%
% Brief code summary
% The input images is first masked (if needed), then the background is
% subtracted (if needed). Next we binarize the image using either scalar or
% adaptive thresholding. Further we extract connected objects, filter them
% according to size, extract the countours of connected blobs and find
% their centers if requested. If needed, the data is either output into the
% parent workspace or into a file.
%
% Usage Warnings:
% Contour extraction feature is not recommended for non-calibration images.
% In principle if the option is activated, the contours will be extracted
% however due to the possibility of overlapping particles the arrays of the
% determined centers and contours in general will not correspond to each
% other.
%
% Headnotes:
%   [1] :   Overflow here is defined by when a pixel attains a maximum
%           value. First of all, if we find such a pixel, there is a good
%           chance that its neighbors are overflowed as well. Secondly, an
%           overflowed values (4095 here) is not a true luminescence value.
%           These two facts make it impossible to correctly locate a
%           particle center. For that reason, the particle is dismissed as
%           soon as we find that it has an overflowed pixel. Note: This
%           applied only to particles whose center is detected by
%           'RadialSym/Gauss' algorithm and not the 'Centroid' algorithm
%   [2] :   For the case of multiple peak connected object , we should
%           check for whether the two peaks are sufficiently  close in
%           intensities. This helps to avoid first, dubious  peaks coming
%           about due to random noise (these will be dimmer than a true
%           particle peak), second cases where a single peak is formed by
%           an overlap of particles. As a rule, such a multiple peak will
%           be brighter than the peaks adjacent to it. We want to avoid
%           dealing with such peaks because the found center location is
%           most likely false. This is the rationale for excluding such
%           particles from analysis.
%   [3] :   A question may arise why not set the masked region of the
%           image to zero? By trial I have found that setting it to zero
%           actually prevents in some cases (where there the background is
%           much illuminated) to detect particles at the boundary of the
%           unmasked region. Thresholds in those regions tend to be very
%           low as a lot of zero-valued pixels contribute to it. To keep it
%           neutral, I decided to instead assign a mean value of the
%           image to the masked pixels.
%
% Input
% imPath - cell array {1,imNum} storing the path to the images
% oprions - struct array with the following fields:
%
%          'batchName'          -   string; name of the image batch. Used as a base
%                                   name for the output files
%          'startNum'           -   int; starting number of the image batch.
%                                   Used for the output files
%          'detectCenters'      -   logical; detect particle centers?
%          'ThreshMethod'       -   Method for thresolding. 'scalar' and
%                                   'adaptive' are supported.
%          'Threshold'          -   <scalar> binarization threshold.
%                                   Signifies sensitivity for adaptive
%                                   thresholding (see imbinarize help) and
%                                   scalar threshold for scalar
%                                   thresholding.
%       'RSG_xPixelOverflowed'  -   logical; exlude from analysis particles
%                                   contaning overflowed pixels. (only for
%                                   RadialSym/Gauss center detection
%                                   algorithm).
%          'minBlobDia'         -   px; min allowable blog diameter to be
%                                   considered a valid object
%          'maxBlobDia'         -   px; max allowable blog diameter to be
%                                   considered a valid object
%          'SubtractBackg'      -   logical; subtract background from the image
%                                   batch (n/a for a single-image batch)
%          'Mask'               -   logical; Use mask for the batch?
%          'CenterFindingAlg'   -   string; algorithm for particle center
%                                   localization. Currently 'RadialSym/Gauss'
%                                   'Centroid' and 'WeightedCentroid' are supported
%          'Filter'             -   string; Supported: 'Gaussian', 'Median'
%       'NumOfFilterIterations' -   [0,1,2,..] - number of smoothing
%                                   iterations to apply to the particle images 
%                                   (effective for large particles).
%          'Average'            -   average the detected particle positions
%                                   over the image batch
%          'StrictnessFactor'   -   Strictness of averaging (see averagePositions help)
%          'detectContours'     -   logical; detect the edges of the
%                                   elliptical dot images (used for perspective
%                                   and distortion bias correction)
%          'write2file'         -   logical; write the center and contours'
%                                   information to a file?
%          'saveDir'            -   directory for the saved files. If not
%                                   indicated, will be requested after execution
%
% Output:
% centers - cell {images}(particles) ; cell array containing particle centers
%           centers{..}(:,1:2) - x,y of the detected particle centers.
%           centers{..}(:,3) - flags of detection. 1 - solitary peak; 2 -
%           part of a multiple-peak cluster.
% contours - cell {images}(particles); cell array containing particle contours
% mask - matrix; binary matrix containing a mask applied to the whole image batch
%
% Example usage:
% Output only particle centers
% [centers]=detectParticles(imPath, options)
% Generate and output a mask and average centers
% [centers,contours, avgCenters, mask]=detectParticles(imPath, options)
% input a ready mask
% [centers, contours]=detectParticles(imPath, options, mask)
%
% -------- Changelog
%  version 5_1, May 25, 2018
% - The masked region of image is now assigned a mean intensity value of the
% image and not zero as before. This allows for more efficient detection of
% calibration dots (see Headnote 3).
% - Some minor bugs fixed.
% - Adaptive thresholding is now more interactive (feedback loop with the
% user).
% - Processing time messages updated.
%
% version 5_2, June 2, 2018
% - centerDetection is now an option
% - some restructuring of the code
% - advanced contourDetection algorithm with automatically adapting
% threshold
% - improved robustness of Gaussian fit
% 
% version 5_3, June 14, 2018
% - option of smoothing added
% 
% version 5_4, June 25, 2018
% - option of averaging particle positions over the image batch added. This
% is useful for static images, such as coming from calibration or dRF
% correction.
% 
% by Peyruz Gasimov, June 25, 2018

function [varargout]=detectParticles5_5(imPath, options, varargin)

imNum=max(size(imPath));

%% Argument check
if ~options.detectCenters && ~options.detectContours
    error('Neither center nor contour detection has been requested in the options. Execution terminated.');
end

if nargout<1 && options.write2file==0
    error('No output method selected for center information. Execution stopped.');
end

if nargout<2 && options.detectContours && ~options.write2file
    warning('No output method selected for contour information. options.detectContours changed to "false".');
    options.detectContours=false;
end

if nargout>=2 && ~options.detectContours
    warning('Invalid output request. Contours will not be detected as per the input options. Contours = NaN.');
    varargout{2}=NaN;
end

if nargout>=3 && ~options.Average
    warning('Invalid output request. Average positions will not be generated as per the input options. avgCenters = NaN.');
    varargout{3}=NaN;
end

if nargout>=4 && ~options.Mask
    warning('Invalid output request. Mask will not be generated as per the input options. Mask = NaN.');
    varargout{4}=NaN;
end

if nargin>=3 && ~options.Mask
    warning('The input mask will be applied to the batch even though the input options.Mask=0.');
    options.Mask=1;
end

%% Preallocate
imArray=cell(1,1,imNum);
centers=cell(1,imNum);

if options.detectContours
    % Preallocate memory
    contours=cell(1,imNum);
end

%% Load Images
for ii=1:imNum
    imArray{1,1,ii}=imread(imPath{ii});
end

imSizeX=size(imArray{1,1,1},2);
imSizeY=size(imArray{1,1,1},1);

%% Mask
if options.Mask
    if nargin>=3
        batchMask=varargin{1};
        if nargout>=4
            warning('Since the mask has been provided at input and required in output, output mask = input mask.');
            varargout{4}=batchMask;
        end
    else
        batchMask=CreateMask(imPath{1});
        if nargout>=4
            varargout{4}=batchMask;
        end
    end
elseif nargout>=4
    varargout{4}=NaN;
    warning('No mask created to be output since options.Mask=false. Assigned NaN.');
end

%% Background subtraction (subtraction of batch minimum for each pixel)
% Image size
imSizex=size(imArray{1,1,1},2);
imSizey=size(imArray{1,1,1},1);

if options.SubtractBackg
    if imNum==1
        warning('Cannot subtract background from a single image. Ignored.')
    else
        batchMin=min(cell2mat(imArray),[],3);
        imArray=cellfun(@minus,imArray,repmat(mat2cell(batchMin,imSizey,imSizex),1,1,imNum),'Un',0);
    end
end

%% Particle center and contour detection
for ii=1:imNum
    %% --------------- Preprocessing --------------- %%
    
    % Apply Mask
    if options.Mask
        imArray{1,1,ii}(~batchMask)=mean(mean(imArray{1,1,ii})); % see Headnote 3
    end
    
    % Image Filter
    if exist('options.Filter','var') && exist('options.NumOfFilterIterations','var')
        if ~isempty(options.Filter) && ~isempty(options.NumOfFilterIterations)
            if options.NumOfFilterIterations>0
                switch options.Filter
                    case 'Gaussian'
                        for ff=1:options.NumOfFilterIterations
                            imArray{1,1,ii}=imgaussfilt(imArray{1,1,ii},0.5);
                        end
                    case 'Median'
                        for ff=1:options.NumOfFilterIterations
                            imArray{1,1,ii}=medfilt2(imArray{1,1,ii});
                        end
                end
            end
        end
    end
    
    % Thresholding
    if strcmp(options.ThreshMethod,'adaptive')
        
        % User feedback loop to find optimum thresholding sensitivity
        if ii==1
            if nargin<4
                
                satis=0;
                while satis==0
                    threshSens=input('Please input the thresholding sensitivity [0,1]: ');
                    BW=imbinarize(imArray{1,1,ii},'adaptive','ForegroundPolarity','bright','Sensitivity',threshSens);
                    fh=figure;
                    imThreshTrial=imArray{1,1,ii};
                    imThreshTrial(BW)=4095;
                    imagesc(imThreshTrial);
                    axis image
                    notClear=true;
                    satis=input('Looks good? (y/n): ','s');
                    while notClear
                        if satis=='y'
                            satis=true;
                            notClear=false;
                            if ishandle(fh)
                                close(fh)
                            end
                        elseif satis=='n'
                            satis=false;
                            notClear=false;
                            if ishandle(fh)
                                close(fh)
                            end
                        else
                            satis=input('Try again. Type "y" or "n": ','s');
                        end
                    end
                end
                if nargout>=5
                    varargout{5}=threshSens;
                end
            else
                threshSens=varargin{2};
                BW=imbinarize(imArray{1,1,ii},'adaptive','ForegroundPolarity','bright','Sensitivity',threshSens);
            end
        else
            BW=imbinarize(imArray{1,1,ii},'adaptive','ForegroundPolarity','bright','Sensitivity',threshSens);
        end
        
        % Scalar Binary thresholding
    elseif strcmp(options.ThreshMethod,'scalar')
        % User feedback loop
        if ii==1
            if nargin<4
                satis=0;
                while satis==0
                    imThresh=input('Please input an integer intensity threshold: ');
                    BW=imbinarize(imArray{1,1,ii},imThresh);
                    fh=figure;
                    imThreshTrial=imArray{1,1,ii};
                    imThreshTrial(BW)=4095;
                    imagesc(imThreshTrial);
                    axis image
                    notClear=true;
                    satis=input('Looks good? (y/n): ','s');
                    while notClear
                        if satis=='y'
                            satis=true;
                            notClear=false;
                            if ishandle(fh)
                                close(fh)
                            end
                        elseif satis=='n'
                            satis=false;
                            notClear=false;
                            if ishandle(fh)
                                close(fh)
                            end
                        else
                            satis=input('Try again. Type "y" or "n": ','s');
                        end
                    end
                end
                if nargout>=5
                    varargout{5}=imThresh;
                end
            else
                imThresh=varargin{2};
                BW=imbinarize(imArray{1,1,ii},imThresh);
            end
        else
            BW=imbinarize(imArray{1,1,ii},imThresh);
        end
    else
        error('Invalid Thresholding Method. options.ThreshMethod supports only "adaptive" and "scalar" methods. See help.');
    end
    
    % Calculate center of particle and bounding box
    if options.detectCenters
        if strcmp(options.CenterFindingAlg,'WeightedCentroid')
            imStats=regionprops('table',BW,imArray{1,1,ii},'WeightedCentroid','BoundingBox','EquivDiameter');
        elseif strcmp(options.CenterFindingAlg,'Centroid')
            imStats=regionprops('table',BW,imArray{1,1,ii},'Centroid','BoundingBox','EquivDiameter');
        elseif strcmp(options.CenterFindingAlg,'RadialSym/Gauss')
            imStats=regionprops('table',BW,imArray{1,1,ii},'BoundingBox','EquivDiameter','PixelList');
        elseif isempty(options.CenterFindingAlg)
            % Only extract info needed for contour detection
            warning('Center finding method not chosen. WeightedCentroid algorithm is used by default.');
            options.CenterFindingAlg='WeightedCentroid';
            imStats=regionprops('table',BW,imArray{1,1,ii},'WeightedCentroid','BoundingBox','EquivDiameter');
        else
            error('Unknown Algorithm for particle center localization. Accepted options are "Centroid", "WeightedCentroid" or "RadialSym/Gauss".');
        end
    else
        imStats=regionprops('table',BW,imArray{1,1,ii},'BoundingBox','EquivDiameter');
    end
    
    % Filter out the blobs according to size
    % if min/max blob size not provided, request user input.
    if isempty(options.minBlobDia) || isempty(options.maxBlobDia)
        fh2=figure;
        histogram(imStats.EquivDiameter);
        title('Blob size histogram')
        options.minBlobDia=input('Please input min blob diameter (px): ');
        options.maxBlobDia=input('Please input max blob diameter (px): ');
        if ishandle(fh2)
            close(fh2)
        end
    end
    
    sizeMask=imStats.EquivDiameter>options.minBlobDia & imStats.EquivDiameter<options.maxBlobDia;
    
    %% --------------- Center Finding --------------- %%
    
    tic
    if options.detectCenters

        if strcmp(options.CenterFindingAlg,'WeightedCentroid')
            imStats=table(imStats.BoundingBox(sizeMask,:), ...
                imStats.EquivDiameter(sizeMask,:), ...
                imStats.WeightedCentroid(sizeMask,:), ...
                'VariableNames',{'BoundingBox','EquivDiameter','WeightedCentroid'});
            
            centers{ii}(:,1:2)=imStats.WeightedCentroid;
            centers{ii}(:,3)=0;
            
        elseif strcmp(options.CenterFindingAlg,'Centroid')
            imStats=table(imStats.BoundingBox(sizeMask,:), ...
                imStats.EquivDiameter(sizeMask,:), ...
                imStats.Centroid(sizeMask,:), ...
                'VariableNames',{'BoundingBox','EquivDiameter','Centroid'});
            
            centers{ii}(:,1:2)=imStats.Centroid;
            centers{ii}(:,3)=0;
            
        elseif strcmp(options.CenterFindingAlg,'RadialSym/Gauss')
            imStats=table(imStats.BoundingBox(sizeMask,:), ...
                imStats.EquivDiameter(sizeMask,:),imStats.PixelList(sizeMask,:), ...
                'VariableNames',{'BoundingBox','EquivDiameter','PixelList'});
            
            % Find the center via radial symmetry algorithm
            % Introduce a lag counter for the cases of multiple overlapping
            % paricles being detected.
            lagCounter=0;
            for kk=1:size(imStats,1)
                
                % Crop the region subimage. Consider only the pixels belonging to the
                % connected object.
                dotIm=imArray{1,1,ii}(ceil(imStats.BoundingBox(kk,2)):floor(imStats.BoundingBox(kk,2)+imStats.BoundingBox(kk,4)),...
                    ceil(imStats.BoundingBox(kk,1)):floor(imStats.BoundingBox(kk,1)+imStats.BoundingBox(kk,3)));
                
                actlix=uint16(sub2ind_2d(size(dotIm,1),imStats.PixelList{kk}(:,2)-floor(imStats.BoundingBox(kk,2)),...
                    imStats.PixelList{kk}(:,1)-floor(imStats.BoundingBox(kk,1))));
                
                mask1=false(size(dotIm));
                mask1(actlix)=true;
                dotIm(~mask1)=0;
                
                % Find the number of peaks
                conn=8;
                peakMask = imregionalmax(dotIm,conn);
                
                peakNum=sum(peakMask(:));
                
                % Check the feasibility of the peaks (remove intensity outliers by median test)
                if peakNum>1
                    
                    % Check if overflow has occurred. (see headnote 1)
                    PixelMax=4095;
                    if options.RSG_xPixelOverflowed
                        if any( dotIm(peakMask) == PixelMax )
                            warning('Pixel overflow. Image #: %i, Particle #: %i. Excluded.', ii, kk);
                            lagCounter=lagCounter-1;
                            continue
                        end
                    end
                    [ix,iy]=find(peakMask);
                    maxPeakDiff=500;    % Maximum allowable deviation from median for a peak intensity. Always >=0. (see headnote 2)
                    medPeak=median(dotIm(peakMask));
                    
                    for mm=1:size(ix,1)
                        peakDiff=abs(dotIm(ix(mm),iy(mm))-medPeak);
                        if peakDiff>maxPeakDiff
                            peakMask(ix(mm),iy(mm))=false;
                        end
                    end
                    peakNum=sum(peakMask(:));
                end
                
                if peakNum>1
                    
                    try
                        [xc,yc,fitError,peakNum]=GaussFit2d(dotIm);
                    catch
                        warning('Problem fitting Gaussian to the blob. Image #: %i, Particle #: %i. Excluded.', ii, kk);
                        lagCounter=lagCounter-1;
                        continue
                    end
                    
                    % Fit error threshold for including the position into the output
                    fitError=sqrt(fitError/numel(dotIm));
                    errThresh=100;
                    
                    
                    if fitError<errThresh
                        centers{ii}(kk+lagCounter:kk+lagCounter+peakNum-1,1)=xc+floor(imStats.BoundingBox(kk,1));
                        centers{ii}(kk+lagCounter:kk+lagCounter+peakNum-1,2)=yc+floor(imStats.BoundingBox(kk,2));
                        centers{ii}(kk+lagCounter:kk+lagCounter+peakNum-1,3)=2;
                        lagCounter=lagCounter+peakNum-1;
                    else
                        lagCounter=lagCounter-1;
                    end
                    
                elseif peakNum==1
                    [cenTemp(1),cenTemp(2)]=radialcenter(double(dotIm));
                    if any(isnan(cenTemp))
                        lagCounter=lagCounter-1;
                        warning('Error finding center via symmetry algorithm. Image #: %i, Particle #: %i. Excluded.', ii, kk);
                        continue
                    else
                        centers{ii}(kk+lagCounter,1)=cenTemp(1);
                        centers{ii}(kk+lagCounter,2)=cenTemp(2);
                        centers{ii}(kk+lagCounter,1:2)=centers{ii}(kk+lagCounter,1:2)+[floor(imStats.BoundingBox(kk,1)), floor(imStats.BoundingBox(kk,2))];
                        centers{ii}(kk+lagCounter,3)=1;
                    end
                else
                    warning('No feasible maximum found in subImage. Image #: %i, Particle #: %i. Excluded.', ii, kk);
                    lagCounter=lagCounter-1;
                end
            end
            
        end
        
        % Display center localization elapsed time
        toc1=toc;
        ImageCounterMsg=['Image (',num2str(ii),'/',num2str(imNum),'): '];
        CLocTimeMsg=['Center localization time: ',num2str(toc1),' sec. \n'];
        fprintf([ImageCounterMsg, CLocTimeMsg]);        
    else
        imStats=table(imStats.BoundingBox(sizeMask,:), 'VariableNames',{'BoundingBox','Equivdiameter'});
    end
    
    %% --------------- Contour Detection --------------- %%
    
    
    if ( (options.detectContours && nargout>=2) || (options.detectContours && options.write2file) )
        
        % For contour detection we will extend the bounding box by 'frame' pixels
        frame=5;
        
        % Set boundary threshold (see subpixelEdges help)
        bndThresh=0.3;
        
        % Extract the ellipse boundaries and
        for kk=1:size(imStats,1)
            ytop=floor(imStats.BoundingBox(kk,2)-frame);
            ybtm=ceil(imStats.BoundingBox(kk,2)+imStats.BoundingBox(kk,4)+2*frame);
            xleft=floor(imStats.BoundingBox(kk,1)-frame);
            xright=ceil(imStats.BoundingBox(kk,1)+imStats.BoundingBox(kk,3)+2*frame);
            
            if ytop<=0
                ytop=1;
            end
            
            if ybtm>imSizeY
                ybtm=imSizeY;
            end
            
            if xleft<=0
                xleft=1;
            end
            
            if xright>=imSizeX
                xright=imSizeX;
            end
            
            dotIm=imArray{1,1,ii}(ytop:ybtm,xleft:xright);
            
            [contours{ii}{kk}] = subpixelEdges(dotIm, bndThresh);
            
            % Check the detected contours
            minDotNum=0.5*pi*imStats.EquivDiameter(kk)/2;
            maxDotNum=2.2*pi*imStats.EquivDiameter(kk)/2;
            
            notDetected=true;
            maxTries=6;
            tries=0;
            bndThresh_0=bndThresh;
            while notDetected
                if isempty(contours{ii}{kk}.x) || numel(contours{ii}{kk}.x)<minDotNum
                    
                    bndThresh=bndThresh*0.7;
                    tries=tries+1;
                    
                    [contours{ii}{kk}] = subpixelEdges(dotIm, bndThresh);
                    if numel(contours{ii}{kk}.x)>=minDotNum && numel(contours{ii}{kk}.x)<maxDotNum
                        notDetected=false;
                    end
                    if tries>=maxTries && notDetected==true
                        contours{ii}{kk}=NaN;
                        break
                    end
                elseif numel(contours{ii}{kk}.x)>=maxDotNum
                    
                    bndThresh=bndThresh*1.3;
                    tries=tries+1;
                    
                    [contours{ii}{kk}] = subpixelEdges(dotIm, bndThresh);
                    if numel(contours{ii}{kk}.x)<maxDotNum && numel(contours{ii}{kk}.x)>=minDotNum 
                        notDetected=false;
                    end
                    if tries>=maxTries && notDetected==true
                        contours{ii}{kk}=NaN;
                        break
                    end
                else
                    notDetected=false;
                end
            end
            % Revert the boundary threshold back to the original value
            bndThresh=bndThresh_0;
            if ~notDetected
                [contours{ii}{kk}] = [contours{ii}{kk}.x contours{ii}{kk}.y]+[floor(imStats.BoundingBox(kk,1))-frame-1, floor(imStats.BoundingBox(kk,2))-frame-1];
            end
        end
        
        % Print Contour detection elapsed time
        if options.detectCenters
            fprintf([repmat(' ',1,length(ImageCounterMsg)+2), 'Contour Detection time: ',num2str(toc-toc1),' sec.\n']);
%                            ^-Alignment
        else
            fprintf('Image (%i/%i): Contour Detection time: %0.1f sec.\n',ii,imNum,toc-toc1);
        end
    end
end

%% --------------- Average positions over the image batch --------------- %%
if options.Average 
    if imNum>1
        fprintf('\n ---- Averaging positions over all images. \n');
        avgCenters=averagePositions2(centers,'StrictnessFactor',options.StrictnessFactor);
    elseif imNum==1
        avgCenters=centers{1};
    end
end
    

%% --------------- Assign Outputs --------------- %%

if nargout>=1
    if options.detectCenters
    varargout{1}=centers;
    else
    warning('Positions output not assigned, since options.detectCenters=false. See detectParticles help.');    
    end
end

if nargout>=2
    if options.detectContours
    varargout{2}=contours;
    else
    warning('Contours output not assigned, since options.detectContours=false. See detectParticles help.');    
    end
end

if nargout>=3 
    if options.Average
    varargout{3}=avgCenters;
    else
        warning('Average positions output not assigned, since options.Average=false. See detectParticles help.');
    end
end

if options.write2file
    % Check the save directory
    if isempty(options.saveDir)
        [saveFolder]=uigetdir(fileparts(imPath{1}),'Save directory.');
        if ~exist(strcat(saveFolder,'\AnalysisPG'),'dir')
            mkdir(saveFolder,'AnalysisPG')
        end
        saveDir=strcat(saveFolder,'\AnalysisPG');
    else
        if ~exist(strcat(options.saveDir,'AnalysisPG'),'dir')
            mkdir(options.saveDir,'AnalysisPG')
        end
        saveDir=strcat(options.saveDir,'\AnalysisPG');
    end
    
    % Write contours to file
    if options.detectContours
        save(strcat(saveDir,'\contours_',options.batchName), 'contours');
    end
    
    % Write mask to file
    if options.Mask
        save(strcat(saveDir,'\mask_',options.batchName), 'batchMask');
    end
    
    % Write particle centers to file
    if isempty(options.startNum)
        options.startNum=0;     % Default value
    end
    
    for ii=1:imNum
        fid=fopen(strcat(saveDir,'\centers_',options.batchName,'_',num2str(options.startNum+ii-1),'.3dp'),'w');
        fprintf(fid, ['ImageID: ',options.batchName,'_',num2str(options.startNum+ii-1), ...
            '   ', datestr(datetime('today')),'   ', 'CenterFindingAlgorithm: ', options.CenterFindingAlg,...
            '   ', 'Mask: ', num2str(options.Mask),'\r\n']);
        fprintf(fid, 'X \t\t Y \t\t Flag\r\n');
        fprintf(fid,'%f\t%f\t%i\r\n',centers{ii}');
        fclose(fid);
    end
    
    if options.Average
        fid=fopen(strcat(saveDir,'\centersAVG_',options.batchName,'.3dp'),'w');
        fprintf(fid, ['Batch Name: ',options.batchName,'   ', 'Averaged over ', num2str(imNum), ' images' ...
            '   ', datestr(datetime('today')),'   ', 'CenterFindingAlgorithm: ', options.CenterFindingAlg,...
            '   ', 'Mask: ', num2str(options.Mask),'\r\n']);
        fprintf(fid, 'X \t\t Y \r\n');
        fprintf(fid,'%f\t%f\r\n',avgCenters');
        fclose(fid);
    end
    
end

end




