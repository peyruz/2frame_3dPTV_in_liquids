%% Function to calculate vectors of the matched particles
% This function is part of the 3d PTV relaxation algorithm
% 
% Input:
% pos1  -   stores positions of particles in snapshot one. size = [N1,3],
%           where N1 - number of particles detected in snapshot one
% pos2  -   stores positions of particles in snapshot two. size = [N2,3],
%           where N2 - number of particles detected in snapshot two
% MatchMatrix - contains information on who are potential Match Candidates
%           of each particle in image one. size = [N1,N2]
% NeighborMatrix - contains information on who are Neighbors of each
%           particle. size = [N1,N1]
% 
% Output:
% velocityVecs - collection of calculated vectors
% MatchRatio
% InvalidMask - logical mask showing vectors recovered as bad during
%               validation
% ValidMatchedMask - a mask showing only vectors passed validation
% 
% Headnote [1]:
% When CheckFor is set to 'Magnitude+Deviation', the code does not
% distinguish between them, ie a vector marked as invalid may be
% marked as such due to discrepancy in either magnitude, or deviation or
% combination of both. 'Deviation' mode of the checker will rely only on the
% angle between a given vector and the local median vector. This is helpful
% in high-gradient settings, since then large discrepancies from local median
% occur labeling vectors wrongfully invalid.
% 
%% by Peyruz Gasimov, April, 2017
%%
function [velocityVecs,MatchRatio, InvalidMask,ValidMatchedMask]=CalcVectors3(pos1,pos2,MatchMatrix,NeighborMatrix,NoMatchP)

N1=size(pos1,1);    % Number of particles in snapshot one

[MatchProb, MatchedParticlesIx]=max(MatchMatrix,[],2);  % A match candidate with max probability is selected

MatchThresh=0.5;       % Probability Threshold over which a match candidate is considered the match (see next line)
MatchedMask=~(NoMatchP>0.5) &  MatchProb>MatchThresh;    % Mask of matched particles         

velocityVecs=zeros(N1,4);
velocityVecs(~MatchedMask,4)=-1; % Flag Unmatched particles

velocityVecs(MatchedMask,1:3)=pos2(MatchedParticlesIx(MatchedMask),:)-pos1(MatchedMask,1:3);  % Calculate velocity vectors for the matched particles

MatchRatio=nnz(MatchedMask)/N1;

%% Vector Validation

% Vector validation checks how different each vector is from its Neighbors
% and if it is too so, marks it invalid (the vector is assigned -2 flag).
% Median mode is recommended since it is less susceptable to outliers which
% is a very useful property for Vector Validation.

% Choose setting s.a. 'mean' or 'median' (recommended)
ValidationMethod='median';

% Check for angular deviation or both deviation and magnitude as compared
% to the median vector? Options: 'Magnitude+Deviation', 'Deviation'. See
% [1].
% CheckFor='Deviation';
CheckFor='Magnitude+Deviation';

% Choose Maximum Allowable Deviation (MaxDev)
% The table below contains angle, deg (left) corresponding to
% MaxDev (right) for a case when both vectors are of the same magnitude.
% e.g. a MaxDev of 0.5 allows ~30 degrees of deviation from Neighbor Median
% direction. MaxDev=SQRT(2-2*COS(deg2rad(Angle)))
% _Angle___MaxDev___
% | 0	|   0   	|
% | 15	|   0.26 	|
% | 30	|   0.51 	|
% | 45	|   0.76 	|
% | 60	|   1    	|
% ------------------|
MaxDev=20;

% CheckFor 'Magnitude+Deviation'
if strcmp(CheckFor,'Magnitude+Deviation')
switch ValidationMethod
    case 'mean'
        for i=1:N1
            if velocityVecs(i,4)~=-1
                VecVariation=1/norm(mean(velocityVecs(NeighborMatrix(i,:),1:3),1))*(norm(mean(velocityVecs(NeighborMatrix(i,:),1:3),1)-velocityVecs(i,1:3)));
                if VecVariation>MaxDev
                    velocityVecs(i,4)=-2;
                end
            end
        end
    case 'median'
        for i=1:N1
            if velocityVecs(i,4)~=-1
                VecVariation=1/norm(median(velocityVecs(NeighborMatrix(i,:),1:3),1))*(norm(median(velocityVecs(NeighborMatrix(i,:),1:3),1)-velocityVecs(i,1:3)));
                if VecVariation>MaxDev
                    velocityVecs(i,4)=-2;
                end
            end
        end
end

% CheckFor 'Deviation'
elseif strcmp(CheckFor,'Deviation')
    switch ValidationMethod
    case 'mean'
        for i=1:N1
            if velocityVecs(i,4)~=-1
                VecVariation=rad2deg(acos((dot(mean(velocityVecs(NeighborMatrix(i,:),1:3),1),velocityVecs(i,1:3)))...
                            /(norm(mean(velocityVecs(NeighborMatrix(i,:),1:3),1))*norm(velocityVecs(i,1:3)))));
                if VecVariation>MaxDev
                    velocityVecs(i,4)=-2;
                end
            end
        end
    case 'median'
        for i=1:N1
            if velocityVecs(i,4)~=-1
                VecVariation=rad2deg(acos((dot(median(velocityVecs(NeighborMatrix(i,:),1:3),1),velocityVecs(i,1:3)))...
                            /(norm(median(velocityVecs(NeighborMatrix(i,:),1:3),1))*norm(velocityVecs(i,1:3)))));
                if VecVariation>MaxDev
                    velocityVecs(i,4)=-2;
                end
            end
        end
    end
else
    error('Unknown CheckFor mode. Valid options are: ''Magnitude+Deviation'' and ''Deviation''.');
end

InvalidMask=velocityVecs(:,4)==-2; 
ValidMatchedMask=all([MatchedMask';~InvalidMask']);

% %% Plot
% figure
% quiver3(pos1(ValidMatchedMask,1), pos1(ValidMatchedMask,2), pos1(ValidMatchedMask,3), velocityVecs(ValidMatchedMask,1),...
%     velocityVecs(ValidMatchedMask,2), velocityVecs(ValidMatchedMask,3),'AutoScale','off');    % Plot velocity vectors
% hold on
% quiver3(pos1(InvalidMask,1), pos1(InvalidMask,2), pos1(InvalidMask,3), velocityVecs(InvalidMask,1),...
%     velocityVecs(InvalidMask,2), velocityVecs(InvalidMask,3),'r','AutoScale','off');    % Plot velocity vectors
% 
% view(0,90)
% axis equal