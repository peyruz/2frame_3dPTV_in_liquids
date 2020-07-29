%% Parsed assignment
% Parse the Cost matrix into several overlapping submatrices and run the hungarian
% algorithm on each. This will not work for a general cost matrix, but is
% written for handling matrices with non-inf values focused around the diagonal.
% 
% Input: 
%   costMat             - <matrix>; cost matrix
%   subDim              - <scalar>; 1st dimension of the submatrices
%   overlap             - <scalar>; number of overlaping rows among the submatrices. For larger       
% Output:
%   finalAssignment     - <vector>; contains the assignment indices  
%
%  by Peyruz Gasimov, Dec 2018

function [finalAssignment]=parsedHungAlg(costMat,subDim,overlap)

%% Analyze the cost matrix and fit a polynomial
% Extract active rows and columns
[r,c]=find(costMat~=inf);
rmin=min(r);
cmax=max(c);

% Polynomial fit
[af,bf,cf]=polyfit(r,c,4);

%% Generate column parseIndices (column boundaries of the submatrices)
parseIndicesRow=[1,subDim];

ii=2;
while max(parseIndicesRow(:))<size(costMat,1)
    parseIndicesRow(ii,:)=parseIndicesRow(ii-1,:)+subDim-overlap;
    ii=ii+1;
end


% Truncate the row parseIndices if exceeded the size of the matrix
parseIndicesRow(parseIndicesRow(:,2)>size(costMat,1),2)=size(costMat,1);
parseIndicesRow(parseIndicesRow(:,1)>=parseIndicesRow(:,2),:)=[];

% Find the corresponding column boundaries
parseIndicesCol=floor(polyval(af,parseIndicesRow,bf,cf));

if cmax+overlap>size(costMat,2)
    parseIndicesCol(parseIndicesCol(:,2)>size(costMat,2),2)=size(costMat,2);
    parseIndicesCol(end,2)=size(costMat,2);
else
    parseIndicesCol(parseIndicesCol(:,2)>size(costMat,2),2)=cmax+overlap;
    parseIndicesCol(end,2)=cmax+overlap;
end

parseIndicesCol(1)=parseIndicesCol(1)-overlap;
parseIndicesCol(parseIndicesCol<1)=1;

parseIndicesRow(parseIndicesCol(:,1)>=parseIndicesCol(:,2),:)=[];
parseIndicesCol(parseIndicesCol(:,1)>=parseIndicesCol(:,2),:)=[];

startInd=sum(parseIndicesCol(:,1)==1);
endInd=size(parseIndicesCol,1)-sum(parseIndicesCol(:,2)==parseIndicesCol(end,2)+1);

parseIndicesCol=parseIndicesCol(startInd:endInd,:);
parseIndicesRow=parseIndicesRow(startInd:endInd,:);


%% Assignment via Hungarian Algorithm
% Create a waitbar
wb=waitbar(0,'Cross-camera particle matching...');

% Preallocate
finalAssignment=zeros(size(costMat,1),1);

for ii=1:size(parseIndicesRow,1)
    
    subMatrix=costMat(parseIndicesRow(ii,1):parseIndicesRow(ii,2),parseIndicesCol(ii,1):parseIndicesCol(ii,2));
    assignment=assignmentoptimal(subMatrix)+parseIndicesCol(ii,1)-1;
    assignment(assignment==(parseIndicesCol(ii,1)-1))=0;
    if ~isempty(assignment)
        if ii==1
            finalAssignment(parseIndicesRow(ii,1):parseIndicesRow(ii,2)-overlap/2)=assignment(1:end-overlap/2);
        else
            finalAssignment(parseIndicesRow(ii,1)+overlap/2:parseIndicesRow(ii,2)-overlap/2)=assignment(overlap/2+1:end-overlap/2);
        end
    else
        continue
    end
    
    waitbar(ii/size(parseIndicesRow,1),wb);
    
end
close(wb);

end
