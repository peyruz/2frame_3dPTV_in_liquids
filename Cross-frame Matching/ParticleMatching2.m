% Function ParticleMatching is part of 3dPTV relaxation algorithm.
% 
% The function takes as an imput the constructed Neighbor and Match
% Matrices to further refine the match between the particles in the first
% and second snapshots. This done using the relaxation algorithm (see eg
% Pereira et al., 2005) wherein the MatchMatrix is successively iterated to
% refined the match probabilities. As the iteration progreses the
% candidates more likely to be a match will gain more probability while the rest
% of the candidates' probabilities will decrease.
% 
% User may vary the number of iterations which usually helps to decrease
% the number of "borderline" cases of match, however does not fully protect
% from thee.
% 
% Changelog
% version 2
% Simply added treatment of particles with no possible match or with no neighbors.
% 
% by Peyruz Gasimov, April, 2017

function [MatchMatrix, NoMatchP]=ParticleMatching2(NeighborMatrix,MatchMatrix,NoMatchP,pos1,pos2,Rq)
tic;

%% Initialize variables
iterNum=15;
A=0.3;
B=10;
N1=size(MatchMatrix,1);
N2=size(MatchMatrix,2);

%% Preallocate memory
Temp1=cell(N1,1);
NeighborIx=cell(N1,1);
NeighborNum=zeros(N1,1);
MatchCandIx=cell(N1,1);
MatchCandNum=zeros(N1,1);

%% Precalculate broadcast vars
for tt=1:N1
    NeighborNum(tt)=nnz(NeighborMatrix(tt,:));   % Number of neighbors of all particles
    NeighborIx{tt}=find(NeighborMatrix(tt,:));   % Indices of neighbors of all particles
    
    MatchCandNum(tt)=nnz(MatchMatrix(tt,:));     % Number of match candidates for all particles i
    MatchCandIx{tt}=find(MatchMatrix(tt,:));     % Indices of match candidates of all particles i    
end

for ttt=1:N1
    Temp1{ttt}=MatchMatrix(NeighborIx{ttt},:);
end
%% Iterations
fprintf('\n----Particle Matching \n');
p = ProgressBar(iterNum);   % Progress bar
for ii=1:iterNum  
    for i=1:N1    
        
        if  NoMatchP(i)==1
            continue
        end
        
        MatchCandNum_i=MatchCandNum(i);    % Number of match candidates for particle i
        MatchCandIx_i=MatchCandIx{i};     % Indices of match candidates of particle i
        
        NeighborNum_i=NeighborNum(i);   % Number of neighbors of particle i
        NeighborIx_i=NeighborIx{i};   % Indices of neighbors of particle i
        
        if  NeighborNum_i==0
            continue
        end
        
        size1_i=[NeighborNum_i,N2];
        
        PosTempx=zeros(size1_i);  % Preallocate memory
        PosTempy=zeros(size1_i);
        PosTempz=zeros(size1_i);
        
        [lin, col]=find(Temp1{i});
        lind=sub2ind_2d(size1_i(1),lin,col);
        
        PosTempx(lind)=pos2(col,1);     % Stores positions of match-cadidates of neighbor particles
        PosTempy(lind)=pos2(col,2);     % I use three separate matrices for x,y and z to avoid using 3d arrrays
        PosTempz(lind)=pos2(col,3);
        
        Mask=Temp1{i}~=0;                  
        Sumask=any(Mask,1);             % We "sum" the masks of the neighbor rows into one [1,N2] vector mask
        
        Temp2=Temp1{i}(:,Sumask);        
        
        PosTempx=PosTempx(:,Sumask);
        PosTempy=PosTempy(:,Sumask);
        PosTempz=PosTempz(:,Sumask);
        
        size2=size(PosTempx,2);
        
        Probabilities=MatchMatrix(i,MatchCandIx_i);
        for j=1:MatchCandNum_i          % Loop on match candidates
            diff=pos2(MatchCandIx_i(j),1:3)-pos1(i,1:3);
            
            MdispNx=pos1(NeighborIx_i,1)+repmat(diff(1),NeighborNum_i,1);   % The code performs a number of calculations which wont be used. 
            MdispNy=pos1(NeighborIx_i,2)+repmat(diff(2),NeighborNum_i,1);   % This was sacrificed in favor of vectorization and readability
            MdispNz=pos1(NeighborIx_i,3)+repmat(diff(3),NeighborNum_i,1);
            
            Qx=repmat(MdispNx,1,size2)-PosTempx;
            Qy=repmat(MdispNy,1,size2)-PosTempy;
            Qz=repmat(MdispNz,1,size2)-PosTempz;
                    
            Q=sqrt(Qx.^2+Qy.^2+Qz.^2)<Rq;       % Calculation of weigh factors Q
            
            Probabilities(j)=Probabilities(j)*(A+B*sum(sum((Temp2.*Q))));
        end
        MatchMatrix(i,MatchCandIx_i)=Probabilities;
        
        sumProb_i=sum(MatchMatrix(i,MatchCandIx_i))+NoMatchP(i);    % Sum of the the received probabilities
        MatchMatrix(i,MatchCandIx_i)=MatchMatrix(i,MatchCandIx_i)/sumProb_i;    % Normalize
        NoMatchP(i)=NoMatchP(i)/sumProb_i;
    end
    p.progress; % Progress bar
end
p.stop; % Progress bar
TimeElapsed=toc;
fprintf('----Done (%1.2f sec)\n',TimeElapsed);