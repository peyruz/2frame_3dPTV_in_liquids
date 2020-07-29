% Script for loading files

dir='Test Data\';           % Indicate the directory of the files
fileBaseName='image';       % Base name to which the id will be attached (eg 'image' in image000, image001 etc)
fileNum=[0:5];              % File counter. Can be either a range (eg [1:10]) or a set of chosen numbers (eg [1 4 5 2])

% Create a text file name, and read the file.

ParticlePositions=cell(length(fileNum),1);    % Preallocate Memory

for k=1:length(fileNum);
    FileName = [dir,'image' num2str(fileNum(k),'%03u') '.txt']; % Construct a filename
    if exist(FileName, 'file')
        FileHandle = fopen(FileName, 'rt');     % Open with read permission
        ParticlePositions{k}=textscan(FileHandle, '%f %f %f');  % Read the contents of the file
        fclose(FileHandle);
    else
        fprintf('File %s was not found.\n', FileName);
    end
end

%-----------------------