% Reveal invalid vectors in a field by comparing it to its neigbors.
% 
% Input:
% positions                    	-   <Nx2, Nx3, 3xN, 2xN>; positions of the vectors (2d or 3d)
% vecs                        	-   <size(vecs)=size(positions)> matrixvectors (2d or 3d)
% Rn                           	-   <scalar>; neighborhood radius
% maxDev (optional)             -   <scalar>; max allowed deviation to be considered valid
% 'Visualize'
% (optional name-value pair)  	-   'on'/'off' (def), Plot the vector field
%                                   at the end
% 
% Output:
% validMask                     -   mask of valid vectors;
%                                   size(validMask)=size(vecs)
% 
% Usage example:
% [validMask]=validateVelField(positions,vecs,Rn,maxDev,'Visualize','on');
% 
% by Peyruz Gasimov, Jun 5, 2018

function [validMask]=validateVelField(positions,vecs,Rn,varargin)

%% Input preprocessing
vecNum=length(vecs);

if vecNum<4
    error('Minimum number of vectors required: 4');
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
addRequired(p,"positions",@(x)(isnumeric(x)&&(size(x,2)==2 || size(x,2)==3)));
addRequired(p,"vecs",@(x)(isnumeric(x)&&(size(x,2)==2 || size(x,2)==3)&&all((size(positions)==size(x)))));
addRequired(p,"Rn",@(x)(isnumeric(x)&&(x>0)));
addOptional(p,"maxDev",0.5,@(x)(isnumeric(x)&&(x>0)));
addParameter(p,'TypeOfTest','median',@(x)(any(validatestring(x,[{'mean'},{'median'}]))));
addParameter(p,'Visualize','off',@(x)(any(validatestring(x,[{'on'},{'off'}]))));
p.CaseSensitive=false;
parse(p,positions,vecs,Rn,varargin{:});

%% Identify the neighbor vectors
% Number of iterations in the search of neighbors (neighbors of neighbors
% of nei...)
neibrSearchOrder=3;

if dim==3
    TR=delaunayTriangulation(positions);
else
    TR=delaunayTriangulation(positions(:,1:2));
end

N=cell(vecNum,1);
for ii=1:vecNum
    N{ii}=ii;
end

for kk=1:neibrSearchOrder
    for ii=1:vecNum
        N{ii}=unique(TR.ConnectivityList(cell2vec(vertexAttachments(TR,N{ii})),:));
        if ~iscolumn(N{ii})
            N{ii}=N{ii}';
        end
        % After the last search iteration, perform filtering based on the
        % input neighbor radius Rn
        if kk==neibrSearchOrder
            dist=positions(N{ii},:)-positions(ii,:);
            dist=sqrt(dist(:,1).^2+dist(:,2).^2+dist(:,3).^2);
            neibrMask=dist<Rn;
            N{ii}=N{ii}(neibrMask);
        end
    end
end

%% Test the vectors
validMask=false(vecNum,1);

if strcmp(p.Results.TypeOfTest,'mean')
    for ii=1:vecNum
        VecVariation=1/norm(mean(vecs(N{ii},:),1))*(norm(mean(vecs(N{ii},:),1)-vecs(ii,:)));
        if VecVariation>p.Results.maxDev
            validMask(ii)=false;
        else
            validMask(ii)=true;
        end
    end
elseif strcmp(p.Results.TypeOfTest,'median')
    for ii=1:vecNum
        VecVariation=1/norm(median(vecs(N{ii},:),1))*(norm(median(vecs(N{ii},:),1)-vecs(ii,:)));
        if VecVariation>p.Results.maxDev
            validMask(ii)=false;
        else
            validMask(ii)=true;
        end
    end
end

% Visualize if requested
if strcmp(p.Results.Visualize,'on')
    figure
    quiver3(positions(validMask,1),positions(validMask,2),positions(validMask,3),vecs(validMask,1),vecs(validMask,2),vecs(validMask,3),'blue','AutoScale','on','AutoScaleFactor',5);
    hold on
    quiver3(positions(~validMask,1),positions(~validMask,2),positions(~validMask,3),vecs(~validMask,1),vecs(~validMask,2),vecs(~validMask,3),'red','AutoScale','on','AutoScaleFactor',5);
    view([0 90]);
    axis image
    axis([-100 6700 -100 4500]);
end

