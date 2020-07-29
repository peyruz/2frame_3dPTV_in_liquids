%% Construct Initial Match Matrix
% 
% We are looking for neighbors of a given particle of snapshot one among
% particles of snapshot two. This matrix will thus contain candidates for
% particle-matching. Match matrix will be constructed in the same way as
% the Neighbor Matrix (see header of BuildNeighborMatrix.m).
% 
% Input :
% pos1 - postions of particles in the first snapshot; Matrix <double>; size: [N1,3]
% pos2 - postions of particles in the first snapshot; Matrix <double>; size: [N2,3]
% Rs - Radius of search
% 
% Output :
% MatchMatrix   :   Matrix <double>; size: [N1,N1]
% NoMatchP      :   Vector <double>; size: [N1,1]
% Rs            :   The final search radius
% 
% About Match Matrix.
% Match Matrix contains match probabilities between particles of the first
% and the second snapshots. An (i,j) entry of the matrix contains the
% probability of match between particles i of the first snapshot and
% particle j of the second snapshot.
% 
% Changelog
% version 4 
% Allow for particles with no match. Assign no-match probability of these
% 100%. This helps dealing with the outliers. Also, the particles in frames
% A and B for which no match was found are removed.
% 
%% by Peyruz Gasimov, April, 2017
%%
function [IniMatchMatrix,NeighborMatrix, IniNoMatchP, pos1, pos2, Rs]=BuildInitialMatchMatrix4(pos1, pos2, Rs,NeighborMatrix)
fprintf('Building Initial Match Matrix..')
tic;

% Settings
minCandNum=1;

% Analyze inputs
N1=size(pos1,1); % Number of particles in snapshot one
N2=size(pos2,1); % Number of particles in snapshot two

%% Preallocate memory
IniMatchMatrix=double(zeros(N1,N2));
IniNoMatchP=double(zeros(N1,1));


%% Initialize variables
frame=0.01; %, mm
sIW=2*Rs+frame;

%%
[IWaddress1]=ceil(pos1(:,1:2)./repmat(sIW,N1,2));    % Find IW address of each particle, snaphot 1
[IWaddress2]=ceil(pos2(:,1:2)./repmat(sIW,N2,2));    % Find IW address of each particle, snaphot 2

for ii=1:N1
    % Search for neighbors (see eligibility criteria at the header of this code)
    IWaddressSearch=IWaddress2-repmat(IWaddress1(ii,1:2),N2,1);    
    NeighborParticles= sum(IWaddressSearch==0 | IWaddressSearch==-1 | IWaddressSearch==1, 2)==2;
    neiNum=sum(NeighborParticles);
    if neiNum==0
        warning('No match candidates found for particle #%i. No-match probability set to 100%',ii);
        IniNoMatchP(ii)=1;
        continue
    end
    
    
    % Refine the Matrix, restricting match candidates to neighborhood radius Rs
    NeighborParticles(NeighborParticles)=normMat((pos2(NeighborParticles,:)-pos1(ii,:))')<Rs; 
    neiNum=sum(NeighborParticles);
    
    % If no match candidates are found, increase Rs
    Rs_new=Rs;
    fail=false;
    while neiNum<minCandNum
        Rs_new=Rs_new*1.1;
        
        NeighborParticles(NeighborParticles)=normMat((pos2(NeighborParticles,:)-pos1(ii,:))')<Rs;
        neiNum=sum(NeighborParticles);
       
        % Terminate when the search radius becomes larger than 1.5 * initial
        % IW size (i.e. 3*Rs). There is no use in the search beyond that
        % radius.
        if Rs_new>3*Rs
            warning('No match candidates found for particle #%i. No-match probability set to 100%',ii);
            IniNoMatchP(ii)=1;
            fail=true;
            break
        end    
    end
    
    if fail
        continue
    end
    
    IniMatchMatrix(ii,NeighborParticles)=1/(neiNum+1);
    IniNoMatchP(ii)=1/(neiNum+1);
    
end

% Clean up the data
pos1(IniNoMatchP==1,:)=[];
pos2(sum(IniMatchMatrix,1)==0,:)=[];
IniMatchMatrix(IniNoMatchP==1,:)=[];
NeighborMatrix(IniNoMatchP==1,:)=[];
NeighborMatrix(:,IniNoMatchP==1)=[];
IniMatchMatrix(:,sum(IniMatchMatrix,1)==0)=[];
IniNoMatchP(IniNoMatchP==1)=[];



% Convert into sparse matrix
IniMatchMatrix=sparse(IniMatchMatrix);

fprintf('Done (%1.2f sec)\n', toc);





