% Find neighbors of each particle  in a particle swarm
% 
% Input:
%   positions       -   < N x 2|3 matrix >; matrix containing row position
%                       vectors. 
%   searchRadius    -   <scalar>; radius of the neighbor search
% 
% by Peyruz Gasimov, Jun 26, 2018

function [NeighborMatrix]=findNeighbors(positions, searchRadius)

p=inputParser;
addRequired(p,'positions',@(x)(ismatrix(x)&&isnumeric(x)&&isreal(x)));
addRequired(p,'SearchRadius',@(x)(isscalar(x)&&isnumeric(x)&&isreal(x)&&(x>0)));

parse(p,positions,searchRadius);

positions=p.Results.positions;
searchRadius=p.Results.SearchRadius;

N1=size(positions,1); % Number of particles

% Determine the dimension of the problem
dimension=size(positions,2);

if dimension~=2 && dimension~=3
    error('Wrong dimension of the position array. <.,2> and <.,3> are accepted dimensions. Try transposing the input position array.')
end

% The code further treats only the 3d case, so:
if dimension==2
    positions=[positions, zeros(N1,1)];
end

%% Preallocate memory
NeighborMatrix=false(N1);

%% Find IW address of each particle
[IWaddress]=ceil((positions(:,1:2)-min(positions(:,1:2),[],1)+0.001)./repmat(searchRadius,N1,2));    

%% Calculate the initial neighborhood radius for each particle
for ii=1:N1
    IWaddressSearch=IWaddress-repmat(IWaddress(ii,1:2),N1,1);    % Search for neighbors (see eligibility criteria at the header of this code)
    NeighborParticlesMask= sum(IWaddressSearch==0 | IWaddressSearch==-1 | IWaddressSearch==1, 2)==2;
    
    temp=abs(positions(NeighborParticlesMask,:)-positions(ii,:));
    NeighborMatrix(ii,NeighborParticlesMask)=(normMat2(temp,'VectorShape','row')<searchRadius)';           % Fill in the Neighbor Matrix row
end

% As of now, particles are listed as their own neighbors in the matrix. Fix that and convert to boolean
NeighborMatrix=boolean(NeighborMatrix-eye(size(NeighborMatrix,1)));