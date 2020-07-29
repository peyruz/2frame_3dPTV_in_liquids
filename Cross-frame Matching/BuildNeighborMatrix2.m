%% Construct Neighbor Matrix

% We start by discretizing the space of interest into subregions or
% Interrogation Windows (IWs). Then we find the IW "address" of every
% particle. Given a particle in IW [i,j], any particle contained in IWs
% [i-1,j-1], [i,j-1], [i+1,j-1], [i,j-1], [i,j], [i1,j+1], [i+1,j-1],
% [i+1,j] and finally [i+1,j+1] will be considered its neighbors. This sort
% of preprocessing will speed up the process of actual Neighbor Search,
% where these preliminary candidates will further filtered based on
% distance to the particle. For finding a match to a particle it is
% critical to find its neighbors. Otherwise a no-match status will be
% assigned to it. For that reason we need to come up with a large enough
% neighbor search radius so that each particle has at least one neighbor.
% This stil does not ensure a match however since the neighbor particles
% may truly themselves, or the particle of matter itself may lack a true
% match.

% Input :
% pos1 - postions of particles in the first snapshot; Matrix <double>; size: [N1,3]
% Rn - Radius of neighborhood

% Output :
% NeighborMatrix; Matrix <boolean>; size: [N1,N1]

%% by Peyruz Gasimov, April, 2017 
%%
function [NeighborMatrix, MeanDist, Rn]=BuildNeighborMatrix2(pos1, Rn)
fprintf('Building Neighbor Matrix..')
tic;
%% Input data analysis

N1=size(pos1,1); % Number of particles in snapshot one

%% Preallocate memory
NeighborMatrix=boolean(zeros(N1));
NeighborParticles=boolean(zeros(1,N1));
MeanDist=zeros(1,N1);
NumOfNbrs=zeros(1,N1);

%% Initialize variables
frame=0.01; %, mm

%%
[IWaddress]=ceil(pos1(:,1:2)./repmat(Rn,N1,2));    % Find IW address of each particle

for i=1:N1
    NeighborRefined=boolean(zeros(1,N1));           % Initialize "carrier" vector
    
    IWaddressSearch=IWaddress-repmat(IWaddress(i,1:2),N1,1);    % Search for neighbors (see eligibility criteria at the header of this code)
    NeighborParticles= sum(IWaddressSearch==0 | IWaddressSearch==-1 | IWaddressSearch==1, 2)==2;
    
    size1=uint32(nnz(NeighborParticles));
    
    NeighborRefined(NeighborParticles)=normMat([pos1(NeighborParticles,:)-repmat(pos1(i,:),size1,1)]')<Rn; % Refine the Matrix, restricting neighbors to neighborhood radius Rn
    NeighborMatrix(i,:)=NeighborRefined;    % Fill in the Neighbor Matrix row
    NumOfNbrs(i)=sum(NeighborMatrix(i,:),2);
end

crit=min(NumOfNbrs)<2;  % The NeighborMatrix will be recalculated with higher Rn until each particle has got at least one neighbor.
                        % We ought to put 2 and not 1 since at this point
                        % the code considers that particle i itself if its
                        % neighbor. This is corrected at the end by
                        % subtracting identity matrix from the
                        % NeighborMatrix
if crit
    Rn=Rn*1.1;
end

while crit
    NeighborMatrix=boolean(zeros(N1));  % Since previous build was a failure, we clear the NeighborMatrix and start from the beginning
    [IWaddress]=ceil(pos1(:,1:2)./repmat(Rn,N1,2));    % Find IW address of each particle
    for i=1:N1      
        NeighborRefined=boolean(zeros(1,N1));           % Initialize "carrier" vector
        
        IWaddressSearch=IWaddress-repmat(IWaddress(i,1:2),N1,1);    % Search for neighbors (see eligibility criteria at the header of this code)
        NeighborParticles= sum(IWaddressSearch==0 | IWaddressSearch==-1 | IWaddressSearch==1, 2)==2;
        
        size1=uint32(size(nonzeros(NeighborParticles),1));
        
        NeighborRefined(NeighborParticles)=normMat([pos1(NeighborParticles,:)-repmat(pos1(i,:),size1,1)]')<Rn; % Refine the Matrix, restricting neighbors to neighborhood radius Rn
        NeighborMatrix(i,:)=NeighborRefined;    % Fill in the Neighbor Matrix row
        NumOfNbrs(i)=sum(NeighborMatrix(i,:),2);
    end
    crit=min(NumOfNbrs)<3;
    if crit
        Rn=1.1*Rn;      % Since previous build was a failure, we increase the neighbor search radius and try again
    end
end

NeighborMatrix=NeighborMatrix-eye(N1);  % The code counts among a given particle's neighbors the particle itself. We correct it at this line.
NeighborMatrix=sparse(boolean(NeighborMatrix));

%% Estimate Mean Distance between particles (probabalistic approach)
dx=max(pos1(:,1))-min(pos1(:,1));
dy=max(pos1(:,2))-min(pos1(:,2));
dz=max(pos1(:,3))-min(pos1(:,3));

if dz==0
    MeanDist=sqrt((dx*dy)/pi/N1);
else
    MeanDist=(3/4*(dx*dy*dz)/pi/N1)^(1/3);
end
ElapsedTime=toc;
fprintf('Done (%1.2f sec)\n', ElapsedTime);