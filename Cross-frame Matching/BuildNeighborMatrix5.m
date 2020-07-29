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
% pos1 -        positions of particles in the first snapshot; Matrix <double>; size: [N1,3]
% IWsize -      <double> Size of interrogation window used to narrow down neighbor search
%               domain. Used to set an upper bound for Rnx,Rny,Rnz.
% Rnx, Rny -    <double> dimensions of neighborhood prism (see changelog version 4
%               for details)
% k -           <int> Approximate number of prisms discretizing the domain
%               depthwise. The Rnz is changed to accomodate for variable
%               aperture
% dimension -   dimensions of the problem. Accepted values: 2,3
% 
%
% Output :
% NeighborMatrix; Matrix <boolean>; size: [N1,N1]
%
% Changelog:
% version 3, Dec, 2017
% - Neighborhood radius will be varied accoring to the local aperture, which
% is calculated based on the minimum and maximum z-coordinate of particles
% contained within the box. Instead of inputting a fixed neighborhood
% radius R_n as in previous version, the function now requires an input of
% an integer factor k, which will be used to calculate R_n with the
% following expression: R_n=b_ij/k, where b_ij is the local aperture
% estimated
%
% version 4, Jan, 2018
% - Search radius (sphere) is replaced by search prism which is
% characterized by dimensions [Rnx,Rny,Rnz]. This setting should be more
% helpful in settings where gradients in one direction are stronger than in
% the other, such as for example in Hele-Shaw or fracture flow. Rnz is
% fixed and is based on the estimate of local aperture. Rnx and Rny are
% enlarged until a sufficient number of neighbors is found.
% + particles with no neighbors found are not causing error, but only a
% warning, thus making the code more robust.
% 
% version 5, Dec, 2018
% Now the particles for which no neighbors have been found are removed.
% 
%% by Peyruz Gasimov, April, 2017
%%
function [NeighborMatrix, pos1, MeanDist, Rnx, Rny, Rnz]=BuildNeighborMatrix5(pos1, IWsize,Rnx,Rny, k)
fprintf('Building Neighbor Matrix..')
tic;

%% Input data analysis

N1=size(pos1,1); % Number of particles in snapshot one
Rnx=repmat(Rnx,1,N1);
Rny=repmat(Rny,1,N1);

% Determine the dimension of the problem
dimension=size(pos1,2);

if dimension~=2 && dimension~=3
    error('Wrong dimension of the position array. <.,2> and <.,3> are accepted dimensions. Try transposing the input position array.')
end

if dimension==2
    pos1=[pos1, zeros(N1,1)];
end

%% Preallocate memory
NeighborMatrix=false(N1);
singlesMask=false(N1,1);
neiNum1=zeros(1,N1);
Rnz=zeros(1,N1);

%% Code settings
% Minimum number of Neighbors required
minNeiNum=2;

%%
[IWaddress]=ceil((pos1(:,1:2)-min(pos1(:,1:2),[],1)+0.001)./repmat(IWsize,N1,2));    % Find IW address of each particle

%% Calculate the initial neighborhood radius for each particle
if dimension==3
    for ii=1:N1
        Rnz(ii)=(max(pos1(all(IWaddress==IWaddress(ii,:),2),3))-min(pos1(all(IWaddress==IWaddress(ii,:),2),3)))/k;
        if Rnz(ii)==0
            Rnz(ii)=Rnz(ii-1);
        end
    end
else
    Rnz=ones(1,N1);
end
NeighborParticles=cell(1,N1);
numNei0=zeros(1,N1);
for ii=1:N1
    NeighborRefined=false(1,N1);           % Initialize "carrier" vector
    
    IWaddressSearch=IWaddress-IWaddress(ii,1:2);    % Search for neighbors (see eligibility criteria at the header of this code)
    NeighborParticles{ii}= sum(IWaddressSearch==0 | IWaddressSearch==-1 | IWaddressSearch==1, 2)==2;
    
    numNei0(ii)=sum(NeighborParticles{ii});
    temp=abs([pos1(NeighborParticles{ii},:)-repmat(pos1(ii,:),numNei0(ii),1)]');
    NeighborRefined(NeighborParticles{ii})=all([temp(1,:)<Rnx(ii); temp(2,:)<Rny(ii); temp(3,:)<Rnz(ii)],1);

    NeighborMatrix(ii,:)=NeighborRefined;
    neiNum1(ii)=sum(NeighborMatrix(ii,:),2)-1;
end

%% Search Volume Expansion
% Increase the search volume for those particles for which a minimum number
% of neighbors has not been found

crit=neiNum1<minNeiNum;   
                            
if any(crit)
    Rnx(crit)=Rnx(crit)*1.1;
    Rny(crit)=Rny(crit)*1.1;
end

% counter=1:N1;

while any(crit)
    
    critIx=find(crit);
    
    for ii=1:length(critIx)
        NeighborRefined=false(1,N1);           % Initialize "carrier" vector
        temp=abs( (pos1(NeighborParticles{critIx(ii)},:)-pos1(critIx(ii),:))' );

        NeighborRefined(NeighborParticles{critIx(ii)})=all([temp(1,:)<Rnx(critIx(ii)); temp(2,:)<Rny(critIx(ii)); temp(3,:)<Rnz(critIx(ii))],1);
        neiNum1(critIx(ii))=sum(NeighborRefined)-1;

        if neiNum1(critIx(ii))>=minNeiNum
            NeighborMatrix(critIx(ii),:)=NeighborRefined;    % Fill in the Neighbor Matrix row
            crit(critIx(ii))=false;
        elseif Rnx(critIx(ii))*1.1>2*IWsize
            singlesMask(critIx(ii))=true;
            warning('The Neighbor Search dimension exceeded the maximum allowable value of 2*IWsize. Consider increasing IWsize or providing more position measurements.')
            crit(critIx(ii))=false;
        else
            Rnx(critIx(ii))=Rnx(critIx(ii))*1.1;
            Rny(critIx(ii))=Rny(critIx(ii))*1.1;
        end
    end
end

% Remove the particles with no found neighbors
NeighborMatrix(singlesMask,:)=[];
NeighborMatrix(:,singlesMask)=[];
pos1(singlesMask,:)=[];
N1=size(NeighborMatrix,1);


% The code counts among a given particle's neighbors the particle itself. We correct it in the next line.
NeighborMatrix=NeighborMatrix-eye(N1);

% Make the Meighbor Matrix symmetric
% NeighborMatrix=(NeighborMatrix+NeighborMatrix')>0;

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