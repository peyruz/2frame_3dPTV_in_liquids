% Reveal invalid vectors in a field using Duncan et al 2010 [1] method.
% 
% Input:
% positions                    	-   <Nx2, Nx3, 3xN, 2xN>; positions of the vectors (2d or 3d)
% vecs                        	-   <size(vecs)=size(positions)> matrixvectors (2d or 3d)
% Rn0                           -   <scalar>; initial neighborhood radius,
%                                   these may grow if not enough neighbors
%                                   are found within the radius
% threshRes                     -   <scalar> - threshold normalized residual value
% tol (optional)                -   <scalar>; tolerance (see [1] for
%                                   details)
% 
% Output:
% validMask                     -   mask of valid vectors;
%                                   size(validMask)=size(vecs)
% varargout{1}                  -   threshold normalized residual value
% 
% Usage example:
% [validMask]=validateVelField_Duncan(positions,vecs,Rn,maxDev,'Visualize','on');
% 
% by Peyruz Gasimov, Jul 6, 2018

%%% Reference:
% [1] J Duncan et al 2010 Meas. Sci. Technol. 21 057002


function [validMask,varargout]=validateVelField_Duncan(positions,vecs,Rn0,varargin)

%% Input preprocessing
vecNum=length(vecs);

if vecNum<5
    error('Minimum number of vectors required: 5');
end

% Convert into a vertical array if needed
if size(vecs,2)==vecNum
    vecs=vecs';
end

if size(positions,2)==length(positions)
    positions=positions';
end

% Dimension
dim=size(vecs,2);

if dim==2
    vecs=[vecs,zeros(vecNum,1)];
    positions=[positions,zeros(vecNum,1)];
end

p=inputParser;
addRequired(p,'positions',@(x)(isnumeric(x)&&(size(x,2)==2 || size(x,2)==3)));
addRequired(p,'vecs',@(x)(isnumeric(x)&&(size(x,2)==2 || size(x,2)==3)&&all((size(positions)==size(x)))));
addRequired(p,'Rn0',@(x)(isnumeric(x)&&(x>0)));
addOptional(p,'threshRes',[],@(x)(isnumeric(x)&&(x>0)));
addOptional(p,'tol',0.1,@(x)(isnumeric(x)&&(x>0)));
addOptional(p,'minNbrsNum',4,@(x)(isnumeric(x)&&(x>0)));
p.CaseSensitive=false;
parse(p,positions,vecs,Rn0,varargin{:});

Rn0=p.Results.Rn0;
vecs=p.Results.vecs;
positions=p.Results.positions;
threshRes=p.Results.threshRes;
tol=p.Results.tol;
minNbrsNum=p.Results.minNbrsNum;

%% Identify the neighbor vectors
% Number of iterations in the search of neighbors (neighbors of neighbors
% of nei...)
neibrSearchOrder=2;

if dim==3
    TR=delaunayTriangulation(positions);
else
    TR=delaunayTriangulation(positions(:,1:2));
end

% Preallocate
NDel=cell(vecNum,1);
N=cell(vecNum,1);
D=N;
for ii=1:vecNum
    NDel{ii}=ii;
end

for kk=1:neibrSearchOrder
    for ii=1:vecNum
        NDel{ii}=unique(TR.ConnectivityList(cell2vec(vertexAttachments(TR,NDel{ii})),:));
        if ~iscolumn(NDel{ii})
            NDel{ii}=NDel{ii}';
        end
        if kk==neibrSearchOrder
            NDel{ii}=NDel{ii}(NDel{ii}~=ii);
        end
    end
end


% Filter out neighbors farther than Rn
Rn=Rn0;
for ii=1:vecNum
    notEnoughNbrs=true;
    while notEnoughNbrs
        dist=positions(NDel{ii},:)-positions(ii,:);
        dist=sqrt(dist(:,1).^2+dist(:,2).^2+dist(:,3).^2);
        N{ii}=NDel{ii}(dist<Rn);
        D{ii}=dist(dist<Rn);
        % If not enough neighbors found, increase the neighbor radius and
        % rerun the loop
        if length(D{ii})<minNbrsNum+1
            Rn=Rn*1.5;
            continue;
        else
            Rn=Rn0;
            notEnoughNbrs=false;
        end
    end
end


% Calculate normalized residuals of the vectors
validMask=false(vecNum,1);
nRes=zeros(vecNum,1);

for ii=1:vecNum
    nRes(ii)=norm( vecs(ii,:)/(median(D{ii})+tol) - median(vecs(N{ii},:)./(D{ii}+tol)) ) /...
        ( median( norm(vecs(N{ii},:)./(D{ii}+tol)-median(vecs(N{ii},:)./(D{ii}+tol))) )+tol );
end

    
%% Interactive thresholding if requested
% if strcmp(p.Results.Interactive,'on')
if isempty(threshRes)
    fhh=figure;
    histogram(nRes);
    
    satis=0;
    while satis==0
        threshRes=input('Please input the threshold value for normalized residual: ');
        
        fh=figure;
        view([0 90]);
        axis image
        hold on
        
        if ishandle(fhh)
            close(fhh);
        end
        
        validMask(nRes<threshRes)=true;
        hq1=quiver3(positions(validMask,1),positions(validMask,2),positions(validMask,3),vecs(validMask,1),vecs(validMask,2),vecs(validMask,3),'blue','AutoScale','off');
        hq2=quiver3(positions(~validMask,1),positions(~validMask,2),positions(~validMask,3),vecs(~validMask,1),vecs(~validMask,2),vecs(~validMask,3),'red','AutoScale','off');
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
                    delete(hq1);
                    delete(hq2);
                end
            else
                satis=input('Try again. Type "y" or "n": ','s');
            end
        end
    end
    if nargout>=2
        varargout{1}=threshRes;
    end
elseif nargout>=2
    varargout{1}=threshRes;
end

