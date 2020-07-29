% Build Cost matrix for Cross Camera matching, via e.g. Hungarian
% Algorithm. Parallelized algorithm.

% Cross Camera Matching is part of a stereo 3dPTV processing. After the
% particles are localized in the images, we need to correspond the same
% particles in Camera 1 and 2 before we can find their 3d position. We
% "dewarp" the detected particle positions and project the them onto the
% z=0 plane (interface) in the physical reference frame. We call this
% projections interface images. Further for each particle in the 1st Camera
% we need to define a region where we could expect the correspondent
% interface image of the 2nd Camera to be. This region is a curve, we call
% "Match curve" (function genMatchCurve). The Match Curve is built using
% Snell's law. After we have the match curve for the particle we find the
% 2nd camera interface images closest to the match curve (function
% dist2QuadCurve) and calculate the distances of each to the curve. These
% distances are input into the cost Matrix. The images which are too far
% from the curve and thus have little chance of being the correct match are
% assigned inf cost.
%
% Input:
% sz01, sz02 -      interface images of particle from Camera 1 and 2
%                   respectively.
% min/maxDepth -    minimum and maximum depths expected in the experiment.
%                   Note that 0>minDepth>maxDepth
% c1,c2 -           positions of camera projections centers in the world reference
%                   frame
% n1,n2 -           refractive indices of the media. For my experiment, 1 - air, 2 -
%                   polyurethane
%
% Output:
% HungCostMatrix -  Cost Matrix inputtable e.g. into Hungarian algorithm
%                   for assignment.

% by Peyruz Gasimov, Sep 2017

function [HungCostMatrix]=Par_buildHungCostMat(sz01,sz02,minDepth,maxDepth,c1,c2,n1,n2)

% Settings
maxAccrete=0;
maxDist=0.15;       % mm
meanDistTh=0.25;    % mm

% Preallocate memory
HungCostMatrix=inf(size(sz01,2),size(sz02,2));

% Triangulation of Cam 2 interface image
tri2=delaunayTriangulation(sz02');

% Find Match Curves and potential match points
% parfor ii=1:size(sz01,2)
parfor ii=1:size(sz01,2)
    % Find match curve
    [PolyFitCoef{ii},x1(ii),x2(ii),MatchCurvePoints{ii}]=genMatchCurve(sz01(:,ii),c1,c2,n1,n2,minDepth,maxDepth);
end

for ii=1:size(sz01,2)
    
    % Check the quality of the points
    if all(isnan(MatchCurvePoints{ii}(:)))
        continue
    end
    
    % Find the Cam2 interface image points close to the Match Curve
    % Initialize
    vi=unique(mat2vec(tri2.ConnectivityList(cell2vec(vertexAttachments(tri2,mat2vec(tri2.ConnectivityList(nonnan(unique(pointLocation(tri2,MatchCurvePoints{ii}(1:2,:)'))),:)))),:)));
    % Accrete
    if maxAccrete>0
        for jj=1:maxAccrete
            vi=unique(mat2vec(tri2.ConnectivityList(cell2vec(vertexAttachments(tri2,vi)),:)));
            dist = dist2QuadCurve(sz02(:,vi)',PolyFitCoef{ii}(1),PolyFitCoef{ii}(2),PolyFitCoef{ii}(3),x1(ii),x2(ii));
            % Stop accretion if mean distance is larger than a threshold
            if mean(dist)>meanDistTh
                break
            end
        end
    else
        dist = dist2QuadCurve(sz02(:,vi)',PolyFitCoef{ii}(1),PolyFitCoef{ii}(2),PolyFitCoef{ii}(3),x1(ii),x2(ii));
    end
    % Filter out the overly long distances
    dist(dist>maxDist)=inf;
    % Fill in a row in the Cost Matrix
    HungCostMatrix(ii,vi)=dist';
end



