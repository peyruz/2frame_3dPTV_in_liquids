% cell array to vertical vector

function [vertVec]=cell2vec(cellArray)

vertVec=0;  % initialize

for ii=1:numel(cellArray)
    
    vertVec=[vertVec;cellArray{ii}(:)];
    
end

vertVec=vertVec(2:end);