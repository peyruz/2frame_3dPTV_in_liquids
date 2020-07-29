%% Construct Neighbor Matrix
% 
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
% 
% Input :
% pos1 - postions of particles in the first snapshot; Matrix <double>; size: [N1,3]
% Rn - Radius of neighborhood
% 
% Output :
% NeighborMatrix; Matrix <boolean>; size: [N1,N1]
% 
% Changelog:
% version 3
% Neighborhood radius will be varied accoring to the local aperture, which
% is calculated based on the minimum and maximum z-coordinate of particles
% contained within the box. Instead of inputting a fixed neighborhood
% radius R_n as in previous version, the function now requires an input of
% an integer factor k, which will be used to calculate R_n with the
% following expression: R_n=b_ij/k, where b_ij is the local aperture
% estimated
% 
%% by Peyruz Gasimov, April, 2017
%%
function [NeighborMatrix, MeanDist, Rn]=BuildNeighborMatrix3(pos1, IWsize, k)
fprintf('Building Neighbor Matrix..')
tic;
%% Input data analysis

N1=size(pos1,1); % Number of particles in snapshot one

%% Preallocate memory
NeighborMatrix=false(N1);
NumOfNbrs=zeros(1,N1);
Rn=zeros(1,N1);

%% Code settings
% Minimum number of Neighbors required
minNeiNum=2;

%%
[IWaddress]=ceil((pos1(:,1:2)-min(pos1(:,1:2),[],1)+0.001)./repmat(IWsize,N1,2));    % Find IW address of each particle

%% Calculate the initial neighborhood radius for each particle
for ii=1:N1
    Rn(ii)=(max(pos1(all(IWaddress==IWaddress(ii,:),2),3))-min(pos1(all(IWaddress==IWaddress(ii,:),2),3)))/k;
    if Rn(ii)==0
        Rn(ii)=Rn(ii-1);
    end
end
NeighborParticles=cell(1,N1);
size1=zeros(1,N1);
for ii=1:N1
    NeighborRefined=boolean(zeros(1,N1));           % Initialize "carrier" vector
    
    IWaddressSearch=IWaddress-repmat(IWaddress(ii,1:2),N1,1);    % Search for neighbors (see eligibility criteria at the header of this code)
    NeighborParticles{ii}= sum(IWaddressSearch==0 | IWaddressSearch==-1 | IWaddressSearch==1, 2)==2;
    
    size1(ii)=sum(NeighborParticles{ii});
    
    NeighborRefined(NeighborParticles{ii})=normMat([pos1(NeighborParticles{ii},:)-repmat(pos1(ii,:),size1(ii),1)]')<Rn(ii); % Refine the Matrix, restricting neighbors to neighborhood radius Rn
    NeighborMatrix(ii,:)=NeighborRefined;    % Fill in the Neighbor Matrix row
    NumOfNbrs(ii)=sum(NeighborMatrix(ii,:),2)-1;
end


crit=NumOfNbrs<minNeiNum;   % The NeighborMatrix will be recalculated with higher Rn until each particle has got at least one neighbor.
                            % We ought to put 2 and not 1 since at this point
                            % the code considers that particle i itself if its
                            % neighbor. This is corrected at the end by
                            % subtracting identity matrix from the
                            % NeighborMatrix
if any(crit)
    Rn(crit)=Rn(crit)*1.1;
end

counter=linspace(1,N1,N1);

while any(crit)

    for ii=counter(crit)
        NeighborRefined=false(1,N1);           % Initialize "carrier" vector
        
        NeighborRefined(NeighborParticles{ii})=normMat([pos1(NeighborParticles{ii},:)-repmat(pos1(ii,:),size1(ii),1)]')<Rn(ii); % Refine the Matrix, restricting neighbors to neighborhood radius Rn
        NumOfNbrs(ii)=sum(NeighborRefined)-1;
        if NumOfNbrs(ii)>=minNeiNum
            NeighborMatrix(ii,:)=NeighborRefined;    % Fill in the Neighbor Matrix row
        end
    end
    crit=NumOfNbrs<minNeiNum;
    if any(crit)
        Rn(crit)=Rn(crit)*1.1;  % Since previous build was a failure, we increase the neighbor search radius and try again
        if any(Rn(crit)>2*IWsize)
            error('The Neighbor Radius exceeded the maximum allowable value of 2*IWsize. Consider increasing IWsize or providing more position measurements.')
        end
    end
end

% The code counts among a given particle's neighbors the particle itself. We correct it in the next line.
NeighborMatrix=NeighborMatrix-eye(N1);  
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