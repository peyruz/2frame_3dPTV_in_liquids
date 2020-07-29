% Create a polygon Image Mask interactively
% 
% Example usage:
% [maskName]=CreateMask; - will ask user to select the image file using gui
% [maskName]=CreateMask(image) - will use 'image'
% 
% by Peyruz Gasimov, Oct 2017
function [maskName]=CreateMask(varargin)


if nargin==0
    %% Upload the Image
    [fileName,PathName,~] = uigetfile({'*.tif*';'*.tiff'},'Load the image to be masked.');
    if fileName==0
        error('Image not chosen. Try again.');
    end
    imPath=strcat(PathName,fileName);
    theImage=imread(imPath);
elseif nargin==1
    theImage=imread(varargin{1});
else
    error('Too many input arguments. I just need a grayscale image');
end

%% Graphical input of the mask vertices
fh=createFigMasking(theImage);
hold on

% Max number of vertices of the mask
maxNum=666;

for ii=1:maxNum
    try
        [x(ii,1),y(ii,1)]=ginput(1);
    catch
        line('XData',[x(1); x(ii-1)],'YData',[y(1); y(ii-1)],'Color', 'yellow');
        break
    end
    
    scatter(x(ii,1),y(ii,1),'yellow','x');
    if ii>1
        line('XData',[x(ii-1); x(ii)],'YData',[y(ii-1); y(ii)],'Color', 'yellow');
        if ii==maxNum
            line('XData',[x(1); x(ii)],'YData',[y(1); y(ii)],'Color', 'yellow');
        end
    end
end

%% Create the binary Mask
maskName=poly2mask(x,y,size(theImage,1),size(theImage,2));

% Close the figure
close(fh);

end