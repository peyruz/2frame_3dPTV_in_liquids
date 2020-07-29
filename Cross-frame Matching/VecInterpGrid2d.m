function [xq,yq,v1q,v2q]=VecInterpGrid2d(pos1,velocityVecs,ValidMatchedMask)

% Interpolate 3dPTV vectors onto a uniform grid

% Construct quiery grid
maxX=max(pos1(:,1));
maxY=max(pos1(:,2));
minY=min(pos1(:,2));
minX=min(pos1(:,1));

XSpan=maxX-minX
YSpan=maxY-minY

nodeY=linspace(0,1,50);
nodeX=linspace(0,1,50);

nodeX=minX+nodeX*XSpan;
nodeY=minY+nodeY*YSpan;

[xq,yq]=meshgrid(nodeX,nodeY);

% Interpolate the vectors
v1q=griddata(pos1(ValidMatchedMask,1),pos1(ValidMatchedMask,2),velocityVecs(ValidMatchedMask,1),xq,yq);
v2q=griddata(pos1(ValidMatchedMask,1),pos1(ValidMatchedMask,2),velocityVecs(ValidMatchedMask,2),xq,yq);

% Plot interpolated velocity field
quiver(xq,yq,v1q,v2q,'AutoScale','off');

% Heatmap of velocity magnitudes
% heatmap(sqrt(v1q.^2+v2q.^2)); 