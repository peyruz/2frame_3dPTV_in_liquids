% Find average positions of particles from a few static images.
%
% Input:
%   posCellArray        -   <cell>{N x 2} (N - number of positions in a
%                           given cell); cell array with each cell
%                           containing detected positions from one of the
%                           images.
%   StrictnessFactor    -   <Name-Value pair>, default=0; not all particles
%                           will show up in all images. StrictnessFactor
%                           controls in how many images the position needs
%                           to show up in order to be considered valid.
%                           StrictnessFactor=0 - a particle will be
%                           considered valid in any case, even if it shows
%                           up only in one image.
%                           StrictnessFactor=1 - a particle needs to show
%                           up in all images to be considered valid.
%                           StrictnessFactor=0.5 - a particle needs to show
%                           up in at least half of the images to be
%                           considered valid.
% 
% Output:
%   posAvgArray         -   <matrix> - array of valid averaged positions
% 
% ---- Changelog
%  version 2, Jun 26, 2018
% - code refactored, now using cluster analysis. Significant speedup (order of 100).
% 
% by Peyruz Gasimov, Jun 26, 2018

function posAvgArray=averagePositions2(posCellArray,varargin)

p=inputParser;
addRequired(p,'posCellArray',@(x)(iscell(x)&&isvector(x)));
addParameter(p,'StrictnessFactor',0,@(x)(isscalar(x)&&(x>=0 && x<=1)));

parse(p,posCellArray,varargin{:});

posCellArray=p.Results.posCellArray;
StrictnessFactor=p.Results.StrictnessFactor;

% Collect all the positions into one array
imNum=numel(posCellArray);
if isrow(posCellArray)
    posCellArray=posCellArray';
end

posArray=cell2mat(posCellArray);

% Trim the position array to 2d
posArray=posArray(:,1:2);

% Cluster Analysis
% Number of clusters is maximum number of detected particles among all
% images
nClust=max(cellfun(@(x)(length(x)),posCellArray));
tic
fprintf('Cluster Analysis..');
[idx,C]=kmeans(posArray,nClust);
fprintf('Done, %f sec \n',toc);

% Generate the output
% Calculate the minimum number of image the particle needs to be detected
% to be considered valid
imEnough=StrictnessFactor*imNum;

posAvgArray=zeros(nClust,2);
for ii=1:nClust
    if sum(idx==ii)>=imEnough
        posAvgArray(ii,:)=C(ii,:);
    else
        posAvgArray(ii,:)=NaN;
    end
end

posAvgArray=posAvgArray(~isnan(posAvgArray(:,1)),:);

% fprintf('RMS= %f px\n',std();
end